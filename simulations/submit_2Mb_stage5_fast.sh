#!/bin/bash
# FAST stage-5-only re-run for A, B, E, F (40 controllers), using the short/6h
# controller_2Mb_stage5.sbatch so they backfill quickly under depleted fairshare.
# Replaces the earlier long/30-day controllers (which were cancelled).
#
# Stages 1-4 are already complete on disk for all 40 (E/F/G cached; A/B had only
# their stage-5 dirs deleted), so each controller runs ONLY stage 5 (plink GLM,
# exact --pca) -- ~3 quick child jobs each.
#
# G is intentionally EXCLUDED: its controllers are still RUNNING the full
# analysis (stage-3 DAP-G) on `long` and self-heal at stage 5.
#
# Run once from a login node:  bash submit_2Mb_stage5_fast.sh

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb
PUBLISH=/n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake/simulations_round_2_2Mb
cd "$REPO"

NREP=10

# submit_cat PREFIX CONFIGFILE SEED_BASE [EXTRA_CONFIG]
submit_cat() {
    local PREFIX=$1 CFG=$2 BASE=$3 EXTRA="${4:-}"
    for N in $(seq 1 "$NREP"); do
        local SEED=$((BASE + N)) ID="${PREFIX}${N}"
        local WD="$SCRATCH/${ID}" PD="$PUBLISH/${ID}"
        mkdir -p "$WD"
        local JID
        JID=$(sbatch --parsable --job-name="${ID}" \
            --output="${WD}/controller_%j.out" --error="${WD}/controller_%j.err" \
            --export=ALL,CONFIGFILE="${CFG}",SEED="${SEED}",WD="${WD}",PD="${PD}",EXTRA_CONFIG="${EXTRA}" \
            controller_2Mb_stage5.sbatch)
        echo "  ${ID} (seed ${SEED}): controller $JID"
    done
}

echo "== A: human, directional-negative =="
submit_cat A config/human_2Mb.yaml 10
echo "== B: human, neutral trait variants =="
submit_cat B config/human_2Mb.yaml 20 "neutral_trait_vars=True"
echo "== E: cattle baseline =="
submit_cat E config/cattle_baseline_from_midpoint_2Mb.yaml 50
echo "== F: cattle + positive selection, bottlenecked =="
submit_cat F config/cattle_sel_bottlenecked_2Mb.yaml 60

echo "Submitted A+B+E+F (40) on short/6h fast controllers. G left running (full analysis)."
