#!/bin/bash
# Run dap-g for a single trait. Extracted from run_revision_dapg_o2.sh so the
# Snakemake stage3_dapg_locus rule can fan out one job per trait.
#
# Required args (all --flag value):
#   --pheno PATH       SBAMS phenotype file (one trait per line, field 2 = trait id)
#   --geno PATH        bgzipped + tabix-indexed geno file (chrom-prepended SBAMS)
#   --trait ID         trait id (e.g. tr1234567); position derived from the suffix
#   --window N         flanking window in bp (default 1_000_000)
#   --ld-ctrl FLOAT    dap-g -ld_control argument
#   --dapg PATH        path to dap-g binary
#   --tabix PATH       path to tabix binary
#   --htslib-lib PATH  htslib lib dir to LD_LIBRARY_PATH
#   --dapg-libpath PATH  extra dir prepended to LD_LIBRARY_PATH for dap-g
#                        (typically the conda env that has libgsl.so.28 etc.)
#   --threads N        threads for dap-g (also exported as OMP_NUM_THREADS, etc.)
#   --out PATH         dap-g output file
#   --log PATH         dap-g log file
# Optional:
#   --pca PATH         SBAMS "controlled" covariate rows to append. Left in place
#                      for ad-hoc use but NOT passed by rules/common.smk any more.
#                      Stage 2 used to supply 20 PCs from ts.pca(), which only
#                      supports mode="branch" -- eigenvectors of the local
#                      genealogy of the very region being fine-mapped. In the
#                      cattle arm (4*Ne*r*L of order 10-40 across 2 Mb) that
#                      spans most of the genotype space and absorbs the causal
#                      variant along with it. See scripts/run_plink_glm.sh for
#                      the measured effect.

set -euo pipefail

PHENO="" ; GENO="" ; PCA="" ; TRAIT=""
WINDOW=1000000 ; LD_CTRL="" ; DAPG="" ; TABIX="" ; HTSLIB_LIB="" ; DAPG_LIBPATH=""
THREADS="${SLURM_CPUS_PER_TASK:-1}"
OUT="" ; LOG=""

while [ $# -gt 0 ]; do
    case "$1" in
        --pheno)         PHENO="$2"; shift 2;;
        --geno)          GENO="$2"; shift 2;;
        --pca)           PCA="$2"; shift 2;;
        --trait)         TRAIT="$2"; shift 2;;
        --window)        WINDOW="$2"; shift 2;;
        --ld-ctrl)       LD_CTRL="$2"; shift 2;;
        --dapg)          DAPG="$2"; shift 2;;
        --tabix)         TABIX="$2"; shift 2;;
        --htslib-lib)    HTSLIB_LIB="$2"; shift 2;;
        --dapg-libpath)  DAPG_LIBPATH="$2"; shift 2;;
        --threads)       THREADS="$2"; shift 2;;
        --out)           OUT="$2"; shift 2;;
        --log)           LOG="$2"; shift 2;;
        *) echo "Unknown arg: $1" >&2; exit 1;;
    esac
done

for v in PHENO GENO TRAIT LD_CTRL DAPG TABIX OUT LOG; do
    [ -n "${!v}" ] || { echo "missing required --$(echo $v | tr A-Z a-z)" >&2; exit 1; }
done

# Prepend any supplied lib paths to LD_LIBRARY_PATH so dap-g finds libgsl
# (the dap-g binary on this cluster needs libgsl.so.28).
[ -n "$DAPG_LIBPATH" ] && export LD_LIBRARY_PATH="${DAPG_LIBPATH}:${LD_LIBRARY_PATH:-}"
[ -n "$HTSLIB_LIB" ]   && export LD_LIBRARY_PATH="${HTSLIB_LIB}:${LD_LIBRARY_PATH:-}"

export OMP_NUM_THREADS="$THREADS"
export OMP_THREAD_LIMIT="$THREADS"
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export MKL_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export BLIS_NUM_THREADS="$THREADS"
export VECLIB_MAXIMUM_THREADS="$THREADS"
export NUMEXPR_NUM_THREADS="$THREADS"

if [ -f "$OUT" ]; then
    echo "skip $TRAIT (output exists at $OUT)"
    exit 0
fi

mkdir -p "$(dirname "$OUT")" "$(dirname "$LOG")"

POSITION="${TRAIT#tr}"
START=$(( POSITION - WINDOW )); [ $START -lt 1 ] && START=1
END=$(( POSITION + WINDOW ))

TASK_TMP="${TMPDIR:-/tmp}/dapg_${SLURM_JOB_ID:-nojid}_${SLURM_ARRAY_TASK_ID:-0}_${TRAIT}"
mkdir -p "$TASK_TMP"
trap 'rm -rf "$TASK_TMP"' EXIT
TMP_SBAMS="$TASK_TMP/locus.sbams"

awk -v t="$TRAIT" '$2 == t' "$PHENO" > "$TMP_SBAMS"
"$TABIX" "$GENO" "1:${START}-${END}" | cut -f3- >> "$TMP_SBAMS"
# Covariates are optional now and the pipeline supplies none. An `if` block
# rather than `[ -n "$PCA" ] && cat ...` because `set -e` is on and a trailing
# false test would abort the script.
if [ -n "$PCA" ]; then
    cat "$PCA" >> "$TMP_SBAMS"
fi

echo "dap-g: $TRAIT (pos $POSITION, $THREADS threads)"
"$DAPG" -d "$TMP_SBAMS" --output_all -t "$THREADS" -msize 5 \
        -o "${OUT}.tmp" -l "${LOG}.tmp" \
        -ld_control "$LD_CTRL" \
    || echo "dap-g exited non-zero for $TRAIT (partial output kept)"

[ -f "${OUT}.tmp" ] && mv -f "${OUT}.tmp" "$OUT"
[ -f "${LOG}.tmp" ] && mv -f "${LOG}.tmp" "$LOG"
