#!/bin/bash
# Run fastEnloc on a (gwas, gtex) pair of DAP-G output directories.
#
# Required args:
#   --gwas-dir DIR      contains outputs/*.dapg.out for the GWAS category
#   --gtex-dir DIR      contains outputs/*.dapg.out for the GTEx category
#   --annot-gwas VCF.gz annotation VCF for the GWAS category
#   --annot-gtex VCF.gz annotation VCF for the GTEx category
#   --out-prefix PREFIX prefix for output files (.enloc.{enrich,gene,mi,sig,snp}.out)
# Optional:
#   --total-snps N      genome-wide SNP count for fastEnloc -t. Use "auto" (the
#                       default) to compute it as the size of the intersection
#                       of column-3 IDs in the two annotation VCFs.
#   --dap2enloc PATH    summarize_dap2enloc.pl (default: looks on PATH)
#   --fastenloc PATH    fastenloc binary (default: looks on PATH)

set -euo pipefail

GWAS_DIR="" ; GTEX_DIR="" ; ANNOT_GWAS="" ; ANNOT_GTEX="" ; PREFIX=""
TOTAL="auto"
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
        --total-snps)       TOTAL="$2"; shift 2;;
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

if [ "$TOTAL" = "auto" ] || [ -z "$TOTAL" ]; then
    TOTAL=$(comm -12 <(ids_in "$ANNOT_GWAS") <(ids_in "$ANNOT_GTEX") | wc -l | tr -d ' ')
    [ "$TOTAL" -gt 0 ] || { echo "ERROR: intersection of annotation VCFs is empty" >&2; exit 1; }
    echo "fastenloc -t (auto): $TOTAL SNPs (intersection of $ANNOT_GWAS and $ANNOT_GTEX)"
fi

OUT_DIR="$(dirname "$PREFIX")"
mkdir -p "$OUT_DIR"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

GWAS_DAP="$TMP/gwas.dap.gz"
GTEX_DAP="$TMP/gtex.dap.gz"

"$DAP2ENLOC" -dir "$GWAS_DIR/outputs" -vcf "$ANNOT_GWAS" | gzip > "$GWAS_DAP"
"$DAP2ENLOC" -dir "$GTEX_DIR/outputs" -vcf "$ANNOT_GTEX" | gzip > "$GTEX_DAP"

"$FASTENLOC" \
    -eqtl "$GTEX_DAP" \
    -gwas "$GWAS_DAP" \
    -annotation "$ANNOT_GWAS" \
    -t "$TOTAL" \
    -prefix "$PREFIX"
