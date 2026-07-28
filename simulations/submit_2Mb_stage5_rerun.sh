#!/bin/bash
# Re-run stage 5 (plink GLM) after fixing run_plink_glm.sh to use exact --pca
# (2 Mb cattle GTEx has only ~121 pruned variants; "--pca approx" needs >=220).
#
# Covers A, B (human) and E, F (cattle). For all four, stages 1-4 are already
# complete on disk, so each controller only re-runs the missing/deleted stage 5.
#   - E/F: stage5 never succeeded (sentinels absent) -> stage5 runs automatically.
#   - A/B: stage5 previously completed with approx PCA; their stage5 dirs must be
#     DELETED first (done by the companion cleanup step) so exact-PCA stage5 re-runs.
#
# G is intentionally EXCLUDED: its controllers are still RUNNING (in stage 3) and
# will invoke the freshly-deployed run_plink_glm.sh when they reach stage 5, so
# they self-heal. Re-submitting G would collide with the live .snakemake lock.
#
# Run once from a login node:  bash submit_2Mb_stage5_rerun.sh
# Nothing lands in priority.

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
            controller_2Mb.sbatch)
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

echo "Submitted A+B+E+F (40 controllers) for stage-5 re-run. G excluded (self-heals)."
