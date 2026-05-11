#!/bin/bash
#
# Submitter: schedule dap-g regional fine-mapping for one (sim, category) on
# Harvard's O2 (SLURM). Spawns up to four jobs with dependencies:
#
#   stage   (cp inputs /n/data2 -> /n/scratch)   [skipped if already staged]
#     |
#   index   (build tabix.gz for the geno file)   [skipped if already built]
#     |
#   array   (run_revision_dapg_o2.sh, one task per CHUNK_SIZE traits)
#     |
#   finalize (rsync outputs back to /n/data2)
#
# Run from an O2 login node (does not require SLURM allocation itself).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER="$SCRIPT_DIR/run_revision_dapg_o2.sh"
FINALIZER="$SCRIPT_DIR/finalize_revision_dapg_o2.sh"

TABIX=/n/data2/hms/dbmi/sunyaev/lab/nconnally/bin/htslib-1.23.1/tabix
BGZIP=/n/data2/hms/dbmi/sunyaev/lab/nconnally/bin/htslib-1.23.1/bgzip
HTSLIB_LIB=/n/data2/hms/dbmi/sunyaev/lab/nconnally/bin/htslib-1.23.1

SIM_DIR=""
CATEGORY=""
LD_CTRL=""
CHUNK_SIZE=""
ARRAY_CONCURRENCY=50
PROFILE=""
SCRATCH_BASE="/n/scratch/users/n/njc12/fine-mapping"
DRY_RUN=0

usage() {
    cat <<EOF
Usage: $0 --sim-dir DIR --category CAT --ld-ctrl FLOAT [options]

Required:
  --sim-dir DIR              e.g. /n/data2/hms/dbmi/sunyaev/lab/nconnally/simulation_fine_mapping/positive_selection_simulation
  --category CAT             e.g. cgwas_scaling_35 (cgwas|cgtex|hgwas|hgtex prefix)
  --ld-ctrl FLOAT            r2 cutoff passed to dap-g -ld_control

Options:
  --chunk-size N             traits per array task (default: profile-dependent)
  --array-concurrency M      simultaneous array tasks (default: $ARRAY_CONCURRENCY)
  --profile P                small|medium|large (default: auto-detect from sample count)
  --scratch-base PATH        scratch root (default: $SCRATCH_BASE)
  --dry-run                  print sbatch commands without submitting
  -h, --help                 this message
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --sim-dir)            SIM_DIR="$2"; shift 2;;
        --category)           CATEGORY="$2"; shift 2;;
        --ld-ctrl)            LD_CTRL="$2"; shift 2;;
        --chunk-size)         CHUNK_SIZE="$2"; shift 2;;
        --array-concurrency)  ARRAY_CONCURRENCY="$2"; shift 2;;
        --profile)            PROFILE="$2"; shift 2;;
        --scratch-base)       SCRATCH_BASE="$2"; shift 2;;
        --dry-run)            DRY_RUN=1; shift;;
        -h|--help)            usage; exit 0;;
        *) echo "Unknown arg: $1" >&2; usage; exit 1;;
    esac
done

[ -z "$SIM_DIR"  ] && { echo "--sim-dir required"  >&2; exit 1; }
[ -z "$CATEGORY" ] && { echo "--category required" >&2; exit 1; }
[ -z "$LD_CTRL"  ] && { echo "--ld-ctrl required"  >&2; exit 1; }
[ -x "$WORKER"   ] || chmod +x "$WORKER"   2>/dev/null || true
[ -x "$FINALIZER" ] || chmod +x "$FINALIZER" 2>/dev/null || true

CAT_ABBV=$(echo "$CATEGORY" | sed -r 's/_.*//g')
SIM_NAME=$(basename "$SIM_DIR")
PUBLISH_DIR="$SIM_DIR/$CAT_ABBV"
PCA_SRC="$SIM_DIR/pca/${CAT_ABBV}.pca"
PHENO_SRC="$PUBLISH_DIR/${CATEGORY}_pheno.sbams"
GENO_SRC="$PUBLISH_DIR/${CATEGORY}_geno.sbams"

for f in "$PHENO_SRC" "$GENO_SRC" "$PCA_SRC"; do
    [ -f "$f" ] || { echo "missing source file: $f" >&2; exit 1; }
done

STAGE_DIR="$SCRATCH_BASE/$SIM_NAME/$CAT_ABBV"
mkdir -p "$STAGE_DIR/pca" "$STAGE_DIR/outputs" "$STAGE_DIR/logs" "$STAGE_DIR/slurm_logs" "$STAGE_DIR/sort_tmp"

# ---------------------------------------------------------------------------
# Manifest of trait ids (one per line). Matches local script's selector
# (run_revision_dapg_regional_tmp_files_local.sh:166): traits are field 2 of
# every pheno line. Order-stable: same line number across resubmissions.
# ---------------------------------------------------------------------------
MANIFEST="$STAGE_DIR/manifest.txt"
awk '{print $2}' "$PHENO_SRC" > "$MANIFEST"
N_TRAITS=$(wc -l < "$MANIFEST")
[ "$N_TRAITS" -gt 0 ] || { echo "no traits found in $PHENO_SRC" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Auto-detect profile from sample count = NF - 2 of the geno file's first
# row (cols are: tag, snpID, then one column per individual).
# ---------------------------------------------------------------------------
N_SAMPLES=$(awk 'NR==1 {print NF - 2; exit}' "$GENO_SRC")
if [ -z "$PROFILE" ]; then
    if   [ "$N_SAMPLES" -le 16000  ]; then PROFILE=small
    elif [ "$N_SAMPLES" -le 80000  ]; then PROFILE=medium
    else                                   PROFILE=large
    fi
fi

case "$PROFILE" in
    small)  PARTITION=short;  CPUS=4;  MEM=16G;  TIME=0-04:00; DEFAULT_CHUNK=8 ;;
    medium) PARTITION=short;  CPUS=8;  MEM=64G;  TIME=0-12:00; DEFAULT_CHUNK=4 ;;
    large)  PARTITION=medium; CPUS=16; MEM=200G; TIME=2-00:00; DEFAULT_CHUNK=1 ;;
    *) echo "Unknown profile: $PROFILE" >&2; exit 1;;
esac
[ -z "$CHUNK_SIZE" ] && CHUNK_SIZE=$DEFAULT_CHUNK

N_ARRAY_TASKS=$(( (N_TRAITS + CHUNK_SIZE - 1) / CHUNK_SIZE ))
ARRAY_SPEC="0-$((N_ARRAY_TASKS - 1))%$ARRAY_CONCURRENCY"

echo "Sim:        $SIM_NAME"
echo "Category:   $CATEGORY  (cat_abbv=$CAT_ABBV)"
echo "Source:     $PUBLISH_DIR"
echo "Stage:      $STAGE_DIR"
echo "Samples:    $N_SAMPLES   Profile: $PROFILE"
echo "Traits:     $N_TRAITS    Chunk: $CHUNK_SIZE   Array tasks: $N_ARRAY_TASKS"
echo "Resources:  -p $PARTITION -c $CPUS --mem=$MEM -t $TIME"
echo "Array spec: --array=$ARRAY_SPEC"

run_or_print() {
    # First arg: a short label used as the synthetic job id in --dry-run.
    local tag="$1"; shift
    if [ "$DRY_RUN" -eq 1 ]; then
        # Print the sbatch command to stderr so it stays human-visible without
        # contaminating the captured "job id" on stdout.
        { printf '  $ '; printf '%q ' "$@"; printf '\n'; } >&2
        echo "DRYJID_${tag}"
    else
        local out
        # sbatch with --parsable prints just "<jobid>" (or "<jobid>;<cluster>").
        out=$("$@")
        echo "$out" | awk -F';' '{print $1}'
    fi
}

# ---------------------------------------------------------------------------
# 1) stage: cp inputs from /n/data2 -> /n/scratch  (skip if already staged)
# ---------------------------------------------------------------------------
need_stage=0
for src in "$PHENO_SRC" "$PCA_SRC"; do
    base=$(basename "$src")
    if [ "$src" = "$PCA_SRC" ]; then dst="$STAGE_DIR/pca/$base"
    else                              dst="$STAGE_DIR/$base"
    fi
    if [ ! -f "$dst" ] || [ "$(stat -c %s "$src" 2>/dev/null || stat -f %z "$src")" != "$(stat -c %s "$dst" 2>/dev/null || stat -f %z "$dst")" ]; then
        need_stage=1
    fi
done
GENO_GZ_DST="$STAGE_DIR/${CATEGORY}_geno.sbams.gz"
GENO_RAW_DST="$STAGE_DIR/${CATEGORY}_geno.sbams"
if [ ! -f "$GENO_GZ_DST" ] && [ ! -f "$GENO_RAW_DST" ]; then
    need_stage=1
fi

JID_STAGE=""
if [ "$need_stage" -eq 1 ]; then
    STAGE_CMD="set -euo pipefail; \
        cp -u '$PHENO_SRC' '$STAGE_DIR/'; \
        cp -u '$PCA_SRC'   '$STAGE_DIR/pca/'; \
        if [ ! -f '$GENO_GZ_DST' ]; then cp -u '$GENO_SRC' '$STAGE_DIR/'; fi; \
        echo staged"
    JID_STAGE=$(run_or_print STAGE sbatch \
        -J dapg_stage_${CAT_ABBV} \
        -p short -t 0-02:00 -c 2 --mem=8G \
        -o "$STAGE_DIR/slurm_logs/stage_%j.out" \
        -e "$STAGE_DIR/slurm_logs/stage_%j.err" \
        --parsable \
        --wrap "$STAGE_CMD")
    echo "stage job:    $JID_STAGE"
else
    echo "stage job:    skipped (inputs already on scratch)"
fi

# ---------------------------------------------------------------------------
# 2) index: build tabix-bgzip geno on scratch  (skip if already built)
# Pipeline mirrors run_revision_dapg_regional_tmp_files_local.sh:94-99.
# ---------------------------------------------------------------------------
JID_INDEX=""
if [ ! -f "${GENO_GZ_DST}.tbi" ] || [ ! -f "$GENO_GZ_DST" ]; then
    INDEX_CPUS=8
    INDEX_MEM=32G
    [ "$PROFILE" = "large" ] && INDEX_CPUS=16 && INDEX_MEM=64G
    INDEX_CMD="set -euo pipefail; \
        export LD_LIBRARY_PATH='$HTSLIB_LIB:\${LD_LIBRARY_PATH:-}'; \
        cd '$STAGE_DIR'; \
        awk 'BEGIN{OFS=\"\\t\"} {pos=\$2; sub(/^snp/,\"\",pos); print \"1\", pos, \$0}' '${CATEGORY}_geno.sbams' \
          | sort -k2,2n --parallel=$INDEX_CPUS -S 4G -T sort_tmp \
          | '$BGZIP' -@ $INDEX_CPUS > '${CATEGORY}_geno.sbams.gz.tmp'; \
        mv '${CATEGORY}_geno.sbams.gz.tmp' '${CATEGORY}_geno.sbams.gz'; \
        '$TABIX' -s 1 -b 2 -e 2 '${CATEGORY}_geno.sbams.gz'; \
        rm -f '${CATEGORY}_geno.sbams'; \
        echo indexed"

    INDEX_DEP=()
    [ -n "$JID_STAGE" ] && INDEX_DEP=(--dependency=afterok:"$JID_STAGE")
    JID_INDEX=$(run_or_print INDEX sbatch \
        -J dapg_index_${CAT_ABBV} \
        -p short -t 0-04:00 -c $INDEX_CPUS --mem=$INDEX_MEM \
        -o "$STAGE_DIR/slurm_logs/index_%j.out" \
        -e "$STAGE_DIR/slurm_logs/index_%j.err" \
        "${INDEX_DEP[@]}" \
        --parsable \
        --wrap "$INDEX_CMD")
    echo "index job:    $JID_INDEX"
else
    echo "index job:    skipped (tabix index already exists)"
fi

# ---------------------------------------------------------------------------
# 3) array: run dap-g per trait
# ---------------------------------------------------------------------------
ARRAY_DEP=()
if   [ -n "$JID_INDEX" ]; then ARRAY_DEP=(--dependency=afterok:"$JID_INDEX")
elif [ -n "$JID_STAGE" ]; then ARRAY_DEP=(--dependency=afterok:"$JID_STAGE")
fi

JID_ARRAY=$(run_or_print ARRAY sbatch \
    -J "dapg_${CAT_ABBV}_${SIM_NAME}" \
    -p $PARTITION -t $TIME -c $CPUS --mem=$MEM \
    --array="$ARRAY_SPEC" \
    -o "$STAGE_DIR/slurm_logs/dapg_%A_%a.out" \
    -e "$STAGE_DIR/slurm_logs/dapg_%A_%a.err" \
    "${ARRAY_DEP[@]}" \
    --export=ALL,STAGE_DIR="$STAGE_DIR",CATEGORY="$CATEGORY",LD_CTRL="$LD_CTRL",CHUNK_SIZE="$CHUNK_SIZE",PUBLISH_DIR="$PUBLISH_DIR" \
    --parsable \
    "$WORKER")
echo "array job:    $JID_ARRAY"

# ---------------------------------------------------------------------------
# 4) finalize: rsync outputs/logs back to /n/data2 (afterany so partial
# results from a failed array still get published).
# ---------------------------------------------------------------------------
JID_FINALIZE=$(run_or_print FINALIZE sbatch \
    --dependency=afterany:"$JID_ARRAY" \
    -o "$STAGE_DIR/slurm_logs/finalize_%j.out" \
    -e "$STAGE_DIR/slurm_logs/finalize_%j.err" \
    --export=ALL,STAGE_DIR="$STAGE_DIR",PUBLISH_DIR="$PUBLISH_DIR" \
    --parsable \
    "$FINALIZER")
echo "finalize job: $JID_FINALIZE"

cat <<EOF

Submitted. Monitor with:
  squeue -u \$USER -j $JID_ARRAY
  tail -f $STAGE_DIR/slurm_logs/dapg_${JID_ARRAY}_0.out
  seff $JID_ARRAY   # after array completes
Outputs land in: $PUBLISH_DIR/outputs/
EOF
