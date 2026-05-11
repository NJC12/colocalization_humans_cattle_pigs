#!/bin/bash
# Run plink2 GLM for one category (GWAS or GTEx).
#
# Mirrors the recipe at the bottom of create_gwas_files_and_phenotypes.py:
#   1. plink2 --vcf ... --make-pgen
#   2. plink2 --pfile ... --maf --indep-pairwise -> ldpruned.prune.in
#   3. plink2 --pfile ... --extract ... --pca approx -> *_pca.eigenvec
#   4. plink2 --pfile ... --glm hide-covar --covar *.eigenvec --pheno ... -> *.glm.linear
#
# Required args:
#   --vcf PATH         (.vcf or .vcf.gz)
#   --pheno PATH       traits.scaling_*.pheno (FID, IID, trait columns)
#   --out-prefix PREF  prefix for outputs (creates PREF.tr*.glm.linear)
# Optional:
#   --maf FLOAT        min allele frequency (default 0.01)
#   --plink2 PATH      plink2 binary (default: plink2)

set -euo pipefail

VCF="" ; PHENO="" ; PREFIX="" ; MAF=0.01
PLINK2="plink2"

while [ $# -gt 0 ]; do
    case "$1" in
        --vcf)         VCF="$2"; shift 2;;
        --pheno)       PHENO="$2"; shift 2;;
        --out-prefix)  PREFIX="$2"; shift 2;;
        --maf)         MAF="$2"; shift 2;;
        --plink2)      PLINK2="$2"; shift 2;;
        *) echo "Unknown arg: $1" >&2; exit 1;;
    esac
done

for v in VCF PHENO PREFIX; do
    [ -n "${!v}" ] || { echo "missing --$(echo $v | tr A-Z a-z)" >&2; exit 1; }
done

OUT_DIR="$(dirname "$PREFIX")"
mkdir -p "$OUT_DIR"

WORK="${PREFIX}.work"
PFILE="${WORK}.pgen"
LD="${WORK}_ldpruned"
PCA="${WORK}_pca"

if [ ! -f "$PFILE" ]; then
    "$PLINK2" --vcf "$VCF" --make-pgen --out "$WORK"
fi
"$PLINK2" --pfile "$WORK" --maf "$MAF" --indep-pairwise 500 0.5 --out "$LD"
"$PLINK2" --pfile "$WORK" --extract "${LD}.prune.in" --pca approx --out "$PCA"
"$PLINK2" --pfile "$WORK" --glm hide-covar \
          --covar "${PCA}.eigenvec" \
          --pheno "$PHENO" \
          --out "$PREFIX"
