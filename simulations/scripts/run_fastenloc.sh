#!/bin/bash
# Run fastEnloc on a (gwas, gtex) pair of DAP-G output directories.
#
# Required args:
#   --gwas-dir DIR      contains <outputs-subdir>/*.dapg.out for the GWAS category
#   --gtex-dir DIR      contains <outputs-subdir>/*.dapg.out for the GTEx category
#   --annot-gwas VCF.gz annotation VCF for the GWAS category
#   --annot-gtex VCF.gz annotation VCF for the GTEx category
#   --out-prefix PREFIX prefix for output files (.enloc.*.out)
# Optional:
#   --prior MODE        "coloc_default" (default) or "estimated".
#
#                       coloc_default passes fastEnloc's --coloc_default_prior,
#                       which fixes p1 = p2 = 1e-4 and p12 = 1e-5 (the coloc
#                       package's defaults, i.e. alpha0 = -9.21, alpha1 = 6.91)
#                       and SKIPS the EM enrichment estimation. That is what we
#                       want for a cross-species comparison: the estimated prior
#                       depends on -total_variants, which under "auto" differed
#                       ~39x between the arms (~66k SNPs for human vs ~1.7k for
#                       cattle), so the two arms' RCPs were not on a common
#                       scale. Expect RCP/LCP to sit higher than estimated-prior
#                       runs -- alpha1 = 6.91 is well above a typical fitted
#                       value of ~4.2.
#
#                       NOTE the flag takes TWO dashes. fastEnloc's own
#                       docs/options.md writes "-coloc_default_prior" with one,
#                       but main.cc only matches "--coloc_default_prior" and
#                       unknown options are fatal (exit 1), not ignored.
#
#                       In coloc_default mode fastEnloc writes only
#                       .enloc.{gene,sig,snp}.out -- .enrich.out and .mi.out
#                       come from the enrichment estimator that is being
#                       skipped. helpers/paths.py:stage4_output_kinds() mirrors
#                       this; keep the two in step.
#   --outputs-subdir SUBDIR  subdir of {gwas,gtex}-dir containing *.dapg.out
#                       files. Defaults to "outputs" so the original layout is
#                       unchanged; set to e.g. "outputs_r2_0_25" to consume a
#                       parallel r2-tagged result set.
#   --total-snps N      genome-wide SNP count for fastEnloc -total_variants.
#                       Ignored entirely in coloc_default mode (fastEnloc never
#                       reads it once p1/p2/p12 are set), so the expensive
#                       "auto" intersection is skipped there.
#   --dap2enloc PATH    summarize_dap2enloc.pl (default: looks on PATH)
#   --fastenloc PATH    fastenloc binary (default: looks on PATH)

set -euo pipefail

GWAS_DIR="" ; GTEX_DIR="" ; ANNOT_GWAS="" ; ANNOT_GTEX="" ; PREFIX=""
OUTPUTS_SUBDIR="outputs"
TOTAL="auto"
PRIOR="coloc_default"
DAP2ENLOC="summarize_dap2enloc.pl"
FASTENLOC="fastenloc"
FASTENLOC_LIBPATH=""

while [ $# -gt 0 ]; do
    case "$1" in
        --gwas-dir)         GWAS_DIR="$2"; shift 2;;
        --gtex-dir)         GTEX_DIR="$2"; shift 2;;
        --annot-gwas)       ANNOT_GWAS="$2"; shift 2;;
        --annot-gtex)       ANNOT_GTEX="$2"; shift 2;;
        --out-prefix)       PREFIX="$2"; shift 2;;
        --outputs-subdir)   OUTPUTS_SUBDIR="$2"; shift 2;;
        --total-snps)       TOTAL="$2"; shift 2;;
        --prior)            PRIOR="$2"; shift 2;;
        --dap2enloc)        DAP2ENLOC="$2"; shift 2;;
        --fastenloc)        FASTENLOC="$2"; shift 2;;
        --fastenloc-libpath) FASTENLOC_LIBPATH="$2"; shift 2;;
        *) echo "Unknown arg: $1" >&2; exit 1;;
    esac
done

# Prepend the fastenloc library path so the binary finds libgsl.so.28 etc.
if [ -n "$FASTENLOC_LIBPATH" ]; then
    export LD_LIBRARY_PATH="${FASTENLOC_LIBPATH}:${LD_LIBRARY_PATH:-}"
fi

for v in GWAS_DIR GTEX_DIR ANNOT_GWAS ANNOT_GTEX PREFIX; do
    [ -n "${!v}" ] || { echo "missing --$(echo $v | tr A-Z a-z | tr _ -)" >&2; exit 1; }
done

# Extract column-3 SNP IDs from a (possibly headerless) gzipped VCF.
# Header-flexible: skips lines starting with '#' (VCF ## meta and #CHROM) AND
# any line whose POS column ($2) is not purely numeric (catches a tab-separated
# header like "chrom\tpos\tid\t..." that has no leading '#'). `gzip -dc` is
# portable across mac (vs. gzcat) and linux (vs. zcat).
ids_in() {
    gzip -dc "$1" \
        | awk '$1 !~ /^#/ && $2 ~ /^[0-9]+$/ && NF >= 3 {print $3}' \
        | sort -u
}

case "$PRIOR" in
    coloc_default|estimated) ;;
    *) echo "ERROR: --prior must be 'coloc_default' or 'estimated', got '$PRIOR'" >&2; exit 1;;
esac

# Only the estimated-prior path reads -total_variants, and computing it "auto"
# means a full sort+comm over both gunzipped annotation VCFs.
if [ "$PRIOR" = "estimated" ] && { [ "$TOTAL" = "auto" ] || [ -z "$TOTAL" ]; }; then
    TOTAL=$(comm -12 <(ids_in "$ANNOT_GWAS") <(ids_in "$ANNOT_GTEX") | wc -l | tr -d ' ')
    [ "$TOTAL" -gt 0 ] || { echo "ERROR: intersection of annotation VCFs is empty" >&2; exit 1; }
    echo "fastenloc -total_variants (auto): $TOTAL SNPs (intersection of $ANNOT_GWAS and $ANNOT_GTEX)"
fi

OUT_DIR="$(dirname "$PREFIX")"
mkdir -p "$OUT_DIR"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

GWAS_DAP="$TMP/gwas.dap.gz"
GTEX_DAP="$TMP/gtex.dap.gz"

"$DAP2ENLOC" -dir "$GWAS_DIR/$OUTPUTS_SUBDIR" -vcf "$ANNOT_GWAS" | gzip > "$GWAS_DAP"
"$DAP2ENLOC" -dir "$GTEX_DIR/$OUTPUTS_SUBDIR" -vcf "$ANNOT_GTEX" | gzip > "$GTEX_DAP"

if [ "$PRIOR" = "coloc_default" ]; then
    echo "fastenloc prior: --coloc_default_prior (p1=p2=1e-4, p12=1e-5; enrichment estimation skipped)"
    "$FASTENLOC" \
        -eqtl "$GTEX_DAP" \
        -gwas "$GWAS_DAP" \
        --coloc_default_prior \
        -prefix "$PREFIX"
else
    echo "fastenloc prior: estimated from the data (-total_variants $TOTAL)"
    "$FASTENLOC" \
        -eqtl "$GTEX_DAP" \
        -gwas "$GWAS_DAP" \
        -total_variants "$TOTAL" \
        -prefix "$PREFIX"
fi
