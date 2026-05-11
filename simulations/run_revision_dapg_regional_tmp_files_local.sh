#$ -S /bin/bash
#$ -e /net/home/nconnally/comparative_colocalization/fastenloc_simulations/eo
#$ -o /net/home/nconnally/comparative_colocalization/fastenloc_simulations/eo
#$ -l mf=48G
#$ -V
#$ -N dapg
#$ -pe serial 16
#$ -binding linear:16

eval "$(/net/apps/conda/miniconda3/bin/conda shell.bash hook)"
conda activate fact

# category should be cgtex, cgwas, hgtex, or hgwas
dir="$1"
category="$2"
startline="$3"
endline="$4"
ld_ctrl="$5"

cat_abbv=$(echo $category | sed -r 's/_.*//g')
cd $dir/$cat_abbv
mkdir -p outputs logs

# ---------------------------------------------------------------------------
# Thread caps -- keep dap-g inside its allocation without pinning CPU IDs.
# OMP_THREAD_LIMIT is the hard ceiling that overrides dap-g's internal
# omp_set_num_threads() call (which is why -t alone wasn't enough).
# ---------------------------------------------------------------------------
TOTAL_CORES=${NSLOTS:-16}
N_PARALLEL=4
THREADS_PER=$(( TOTAL_CORES / N_PARALLEL ))

export OMP_NUM_THREADS=$THREADS_PER
export OMP_THREAD_LIMIT=$THREADS_PER
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export MKL_NUM_THREADS=$THREADS_PER
export OPENBLAS_NUM_THREADS=$THREADS_PER
export BLIS_NUM_THREADS=$THREADS_PER
export VECLIB_MAXIMUM_THREADS=$THREADS_PER
export NUMEXPR_NUM_THREADS=$THREADS_PER

# ---------------------------------------------------------------------------
# Node-local scratch on /local/var/tmp (RAID0 ext4, NOT shared between
# nodes, NOT redundant). Everything we put in here we are responsible for
# cleaning up. Each job gets its own subdirectory keyed on JOB_ID + PID
# so the cleanup at the end can never touch another job's files.
# ---------------------------------------------------------------------------
LOCAL_ROOT="/local/var/tmp"
LOCAL_SCRATCH="${LOCAL_ROOT}/${USER}_dapg_${JOB_ID:-nojid}_$$"
mkdir -p "$LOCAL_SCRATCH/tmp_sbams" "$LOCAL_SCRATCH/sort_tmp" \
         "$LOCAL_SCRATCH/outputs"  "$LOCAL_SCRATCH/logs"

# Best-effort copy of any in-flight dap-g outputs back to NFS, then nuke
# the local scratch. Runs on normal exit and on SIGINT/SIGTERM.
salvage_and_cleanup() {
    # Don't let salvage interrupt itself
    trap '' INT TERM
    # Move any completed-but-not-yet-moved outputs back to NFS
    if [ -d "$LOCAL_SCRATCH/outputs" ]; then
        find "$LOCAL_SCRATCH/outputs" -type f -name '*.dapg.out' \
            -exec mv -n {} outputs/ \; 2>/dev/null || true
    fi
    if [ -d "$LOCAL_SCRATCH/logs" ]; then
        find "$LOCAL_SCRATCH/logs"    -type f -name '*.dapg.out' \
            -exec mv -n {} logs/    \; 2>/dev/null || true
    fi
    rm -rf "$LOCAL_SCRATCH"
}
trap 'kill $(jobs -p) 2>/dev/null; salvage_and_cleanup' INT TERM
trap salvage_and_cleanup EXIT

# ---------------------------------------------------------------------------
# Stage (or build) the tabix-indexed genotype file on node-local disk.
#
# Two cases:
#   (a) NFS already has a built ${category}_geno.sbams.gz + .tbi -- copy
#       both to local scratch once at job start.
#   (b) NFS doesn't have it yet -- build it on local scratch, then publish
#       to NFS so other jobs can reuse it instead of rebuilding.
# Either way, GENO_GZ ends up pointing at the local copy and every tabix
# query during the run hits local disk only.
# ---------------------------------------------------------------------------
GENO_RAW="${category}_geno.sbams"
NFS_GENO_GZ="${category}_geno.sbams.gz"
LOCAL_GENO_GZ="$LOCAL_SCRATCH/${category}_geno.sbams.gz"

if [ -f "$NFS_GENO_GZ" ] && [ -f "${NFS_GENO_GZ}.tbi" ]; then
    echo "Staging existing tabix index from NFS to $LOCAL_SCRATCH ..."
    cp "$NFS_GENO_GZ"        "$LOCAL_GENO_GZ"
    cp "${NFS_GENO_GZ}.tbi"  "${LOCAL_GENO_GZ}.tbi"
else
    echo "Building tabix index for $GENO_RAW on local scratch ..."
    awk 'BEGIN{OFS="\t"} {pos=$2; sub(/^snp/,"",pos); print "1", pos, $0}' \
        "$GENO_RAW" \
      | sort -k2,2n --parallel=8 -S 4G -T "$LOCAL_SCRATCH/sort_tmp" \
      | bgzip -@ 8 > "${LOCAL_GENO_GZ}.tmp"
    mv "${LOCAL_GENO_GZ}.tmp" "$LOCAL_GENO_GZ"
    tabix -s 1 -b 2 -e 2 "$LOCAL_GENO_GZ"

    echo "Publishing built index back to NFS for reuse by other jobs ..."
    cp "$LOCAL_GENO_GZ"       "${NFS_GENO_GZ}.tmp"        \
        && mv "${NFS_GENO_GZ}.tmp"        "$NFS_GENO_GZ"
    cp "${LOCAL_GENO_GZ}.tbi" "${NFS_GENO_GZ}.tbi.tmp"    \
        && mv "${NFS_GENO_GZ}.tbi.tmp"    "${NFS_GENO_GZ}.tbi"
fi

GENO_GZ="$LOCAL_GENO_GZ"   # all queries below hit local disk

# Stage the phenotype sbams to local disk too -- it gets awk-scanned once
# per trait, and on a 16-trait batch even a small NFS file is worth caching.
LOCAL_PHENO="$LOCAL_SCRATCH/${category}_pheno.sbams"
cp "${category}_pheno.sbams" "$LOCAL_PHENO"

# Same for the PCA file (concatenated into every per-trait sbams).
LOCAL_PCA="$LOCAL_SCRATCH/${cat_abbv}.pca"
cp "../pca/${cat_abbv}.pca" "$LOCAL_PCA"

# ---------------------------------------------------------------------------
# Worker: fine-map a single trait. Reads everything from local disk, writes
# dap-g's output and log to local disk, then promotes them to NFS.
# ---------------------------------------------------------------------------
process_trait() {
    local ph=$1
    local nfs_out="outputs/${category}_${ph}.dapg.out"
    local nfs_log="logs/${category}_${ph}.dapg.out"

    if [ -f "$nfs_out" ]; then
        echo "skip $ph (output exists)"
        return
    fi

    local position=${ph#tr}
    local start=$((position - 1000000))
    [ $start -lt 1 ] && start=1
    local end=$((position + 1000000))

    local tmp="$LOCAL_SCRATCH/tmp_sbams/${category}_${ph}.sbams"
    local local_out="$LOCAL_SCRATCH/outputs/${category}_${ph}.dapg.out"
    local local_log="$LOCAL_SCRATCH/logs/${category}_${ph}.dapg.out"

    awk -v trait="$ph" '$2 == trait' "$LOCAL_PHENO" > "$tmp"
    tabix "$GENO_GZ" "1:${start}-${end}" | cut -f3- >> "$tmp"
    cat "$LOCAL_PCA" >> "$tmp"

    echo "dap-g: $ph (pos $position, ${THREADS_PER} threads)"
    dap-g -d "$tmp" --output_all -t $THREADS_PER -msize 5 \
          -o "$local_out" -l "$local_log" \
          -ld_control "$ld_ctrl" \
        || echo "dap-g exited non-zero for $ph (partial output kept)"

    # Promote results to NFS as soon as they're done, so a killed job
    # still leaves completed traits on shared storage.
    [ -f "$local_out" ] && mv "$local_out" "$nfs_out"
    [ -f "$local_log" ] && mv "$local_log" "$nfs_log"

    rm -f "$tmp"
}

# ---------------------------------------------------------------------------
# Dispatch traits N_PARALLEL at a time.
# ---------------------------------------------------------------------------
slot=0
pids=()
for ph in $(awk -v s=$startline -v e=$endline \
            'NR>s && NR<=e {print $2}' "$LOCAL_PHENO"); do

    process_trait "$ph" &
    pids[$slot]=$!
    slot=$((slot + 1))

    if [ $slot -ge $N_PARALLEL ]; then
        for p in "${pids[@]}"; do wait "$p"; done
        pids=()
        slot=0
    fi
done

for p in "${pids[@]}"; do wait "$p"; done
