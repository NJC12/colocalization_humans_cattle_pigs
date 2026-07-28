#!/bin/bash
#SBATCH --partition=short
#SBATCH --time=6:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#
# FAST stage-5-only controller variant of controller_2Mb.sbatch. Used to re-run
# ONLY stage 5 (plink GLM, exact --pca) for replicates whose stages 1-4 are
# already complete on disk. Differs from controller_2Mb.sbatch in two ways, both
# to maximize backfill scheduling under depleted fairshare:
#   1. short / 6 h  (not long / 30 days): the stage-5 pipeline is 3 quick plink
#      jobs (~10 min compute), so a 30-day reservation was un-backfillable.
#   2. --set-resources trims the child stage5_plink_glm job from short/8h/32GB to
#      short/1h/4GB so those child jobs also fit backfill gaps (plink on 2 Mb is
#      seconds-minutes; the smoke test ran in ~2 s).
# Stages 1-4 are cached, so Snakemake only executes stage 5. If a replicate's
# stage 1-4 outputs were somehow missing this controller could time out at 6 h;
# that is not the case for the A/B/E/F re-run (1-4 present, only stage 5 to do).
#
# Parameterized via --export env vars (same as controller_2Mb.sbatch):
#   CONFIGFILE, SEED, WD, PD, and optional EXTRA_CONFIG.

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SNAKEMAKE=/home/njc12/miniconda3/envs/coloc_sims/bin/snakemake

: "${CONFIGFILE:?}" "${SEED:?}" "${WD:?}" "${PD:?}"
EXTRA_CONFIG="${EXTRA_CONFIG:-}"
mkdir -p "$WD"
cd "$WD"

# --unlock releases any lock from the cancelled long controller before re-running.
"$SNAKEMAKE" --snakefile "$REPO/Snakefile" --configfile "$REPO/$CONFIGFILE" --unlock \
    --config stage1_seed="$SEED" stage2_seed="$SEED" workdir="$WD" publishdir="$PD" $EXTRA_CONFIG || true

# Same --set-resources as the full controller (stage2 -> medium, human stage1 ->
# long) PLUS a tight stage5_plink_glm profile (60 min / 4 GB) for fast backfill.
"$SNAKEMAKE" --snakefile "$REPO/Snakefile" --configfile "$REPO/$CONFIGFILE" \
    --profile "$REPO/profiles/o2" --rerun-triggers mtime \
    --set-resources stage2_split_pheno:slurm_partition=medium stage1_human:slurm_partition=long \
                    stage5_plink_glm:runtime=60 stage5_plink_glm:mem_mb=4000 \
    --config stage1_seed="$SEED" stage2_seed="$SEED" workdir="$WD" publishdir="$PD" $EXTRA_CONFIG
