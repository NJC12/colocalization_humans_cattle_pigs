#!/bin/bash
# Run plink2 GLM for one category (GWAS or GTEx).
#
#   1. plink2 --vcf ... --make-pgen
#   2. plink2 --pfile ... --maf --glm hide-covar allow-no-covars --pheno ... -> *.glm.linear
#
# NO PCA COVARIATES, deliberately. This used to LD-prune, run --pca, and pass the
# eigenvectors as --covar. Both simulated arms are single panmictic populations
# (the human arm samples only EUR) and the environmental noise is drawn i.i.d.,
# so there is no confounding to correct -- the PCs could only remove signal, and
# because they were computed from inside the 2 Mb test region they removed a lot
# of it. Measured on round-2 output, at leads where the lead SNP IS the causal
# SNP, SE_observed / sqrt(1/(2pq*n)) was 2.05 (cattle GWAS) and 2.29 (cattle
# GTEx) against 1.10 / 1.11 for human -- i.e. ~76-81% of the causal variant's
# genotype variance absorbed in cattle versus ~18% in human, deflating cattle
# chi-square 4-5x. Cattle's 2 Mb region has 4*Ne*r*L of order 10-40 versus ~10^3
# in human EUR, so a handful of local PCs span most of its genotype space.
# If a structure control is ever needed, compute it from an UNLINKED region.
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

if [ ! -f "$PFILE" ]; then
    "$PLINK2" --vcf "$VCF" --make-pgen --out "$WORK"
fi
# allow-no-covars is mandatory, not cosmetic: plink2 refuses to run --glm with no
# covariate file unless you say so explicitly ("Since running --glm without at
# least e.g. principal component covariates is usually an analytical mistake").
# Here it is not a mistake -- see the header.
#
# --maf now applies to the GLM itself, not just to the (removed) pruning step.
# Without it the GLM tested every variant including singletons -- 208k instead of
# ~9k in the round-2 human sample -- and the lead-SNP awk below took its minimum
# p-value over that whole set.
"$PLINK2" --pfile "$WORK" --maf "$MAF" --glm hide-covar allow-no-covars \
          --pheno "$PHENO" \
          --out "$PREFIX"

# Collect the lead SNP (smallest p-value, column 15) from each per-trait
# *.glm.linear file. ENDFILE fires once per input file; the resulting TSV has
# one row per trait. NA / non-numeric p-values are skipped (string vs. number
# comparison fails). Trivial CPU/memory.
LEAD_OUT="${PREFIX}_lead_snps.tsv"
if compgen -G "${PREFIX}"*.glm.linear >/dev/null; then
    gawk 'BEGIN {OFS="\t"; out=""} FNR == 1 {min_p=1; out=""} $15 < min_p {min_p=$15; out=$0} ENDFILE {print FILENAME"\t"out}' \
        "${PREFIX}"*.glm.linear > "${LEAD_OUT}"
    echo "Lead SNPs: wrote $(wc -l < "${LEAD_OUT}") rows to ${LEAD_OUT}"
else
    echo "No ${PREFIX}*.glm.linear files found — skipping lead SNP collection" >&2
fi
