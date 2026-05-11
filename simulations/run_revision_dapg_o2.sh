#!/bin/bash
#SBATCH -J dapg_array
# Resource directives (-p, -c, --mem, -t, --array, -o, -e) are supplied by the
# submitter at sbatch time so the same script handles small/medium/large
# profiles without duplication.
#
# Per-array-task worker: processes CHUNK_SIZE traits sequentially, each via the
# same dap-g call as run_revision_dapg_regional_tmp_files_local.sh.
#
# Required env (passed via sbatch --export):
#   STAGE_DIR        /n/scratch/.../<sim>/<cat_abbv>
#   CATEGORY         e.g. cgwas_scaling_35
#   LD_CTRL          dap-g -ld_control argument
#   CHUNK_SIZE       traits per array task
#   PUBLISH_DIR      /n/data2/.../<sim>/<cat_abbv>  (for "skip if already done")

set -euo pipefail

TABIX=/n/data2/hms/dbmi/sunyaev/lab/nconnally/bin/htslib-1.23.1/tabix
DAPG=/n/data2/hms/dbmi/sunyaev/lab/nconnally/bin/dap/dap_src/dap-g
HTSLIB_LIB=/n/data2/hms/dbmi/sunyaev/lab/nconnally/bin/htslib-1.23.1
export LD_LIBRARY_PATH="${HTSLIB_LIB}:${LD_LIBRARY_PATH:-}"

: "${STAGE_DIR:?STAGE_DIR not set}"
: "${CATEGORY:?CATEGORY not set}"
: "${LD_CTRL:?LD_CTRL not set}"
: "${CHUNK_SIZE:?CHUNK_SIZE not set}"
: "${PUBLISH_DIR:?PUBLISH_DIR not set}"

CAT_ABBV=$(echo "$CATEGORY" | sed -r 's/_.*//g')

PHENO="$STAGE_DIR/${CATEGORY}_pheno.sbams"
GENO_GZ="$STAGE_DIR/${CATEGORY}_geno.sbams.gz"
PCA="$STAGE_DIR/pca/${CAT_ABBV}.pca"
MANIFEST="$STAGE_DIR/manifest.txt"
OUT_DIR="$STAGE_DIR/outputs"
LOG_DIR="$STAGE_DIR/logs"

for f in "$PHENO" "$GENO_GZ" "${GENO_GZ}.tbi" "$PCA" "$MANIFEST"; do
    [ -f "$f" ] || { echo "missing input: $f" >&2; exit 1; }
done
mkdir -p "$OUT_DIR" "$LOG_DIR"

# Thread caps -- keep dap-g inside its allocation. Reused verbatim from the
# muhee local script: OMP_THREAD_LIMIT is the hard ceiling that overrides
# dap-g's internal omp_set_num_threads() call.
THREADS_PER=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$THREADS_PER
export OMP_THREAD_LIMIT=$THREADS_PER
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export MKL_NUM_THREADS=$THREADS_PER
export OPENBLAS_NUM_THREADS=$THREADS_PER
export BLIS_NUM_THREADS=$THREADS_PER
export VECLIB_MAXIMUM_THREADS=$THREADS_PER
export NUMEXPR_NUM_THREADS=$THREADS_PER

# SLURM gives us a per-job TMPDIR on node-local disk on O2; fall back to /tmp.
TASK_TMP="${TMPDIR:-/tmp}/dapg_${SLURM_JOB_ID:-nojid}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TASK_TMP"
trap 'rm -rf "$TASK_TMP"' EXIT

N_TRAITS=$(wc -l < "$MANIFEST")
START_IDX=$(( ${SLURM_ARRAY_TASK_ID:-0} * CHUNK_SIZE + 1 ))
END_IDX=$(( START_IDX + CHUNK_SIZE - 1 ))
[ $END_IDX -gt $N_TRAITS ] && END_IDX=$N_TRAITS

echo "[$(date -Is)] task ${SLURM_ARRAY_TASK_ID:-0}: traits $START_IDX-$END_IDX of $N_TRAITS, ${THREADS_PER} threads"

process_trait() {
    local ph=$1
    local out="$OUT_DIR/${CATEGORY}_${ph}.dapg.out"
    local log="$LOG_DIR/${CATEGORY}_${ph}.dapg.out"
    local published_out="$PUBLISH_DIR/outputs/${CATEGORY}_${ph}.dapg.out"

    if [ -f "$out" ] || [ -f "$published_out" ]; then
        echo "skip $ph (output exists)"
        return
    fi

    local position=${ph#tr}
    local start=$(( position - 1000000 ))
    [ $start -lt 1 ] && start=1
    local end=$(( position + 1000000 ))

    local tmp="$TASK_TMP/${CATEGORY}_${ph}.sbams"

    awk -v t="$ph" '$2 == t' "$PHENO" > "$tmp"
    "$TABIX" "$GENO_GZ" "1:${start}-${end}" | cut -f3- >> "$tmp"
    cat "$PCA" >> "$tmp"

    echo "dap-g: $ph (pos $position, ${THREADS_PER} threads)"
    "$DAPG" -d "$tmp" --output_all -t "$THREADS_PER" -msize 5 \
            -o "${out}.tmp" -l "${log}.tmp" \
            -ld_control "$LD_CTRL" \
        || echo "dap-g exited non-zero for $ph (partial output kept)"

    [ -f "${out}.tmp" ] && mv -f "${out}.tmp" "$out"
    [ -f "${log}.tmp" ] && mv -f "${log}.tmp" "$log"
    rm -f "$tmp"
}

while IFS= read -r ph; do
    [ -z "$ph" ] && continue
    process_trait "$ph"
done < <(sed -n "${START_IDX},${END_IDX}p" "$MANIFEST")

echo "[$(date -Is)] task ${SLURM_ARRAY_TASK_ID:-0}: done"
