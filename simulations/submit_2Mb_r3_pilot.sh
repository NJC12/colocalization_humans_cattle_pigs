#!/bin/bash
# Round-3 PILOT: one human replicate (A1) + one cattle replicate (E1).
#
# Run this before submit_2Mb_r3.sh. The point is to exercise the full stage
# 1->5 path on the round-3 code once, cheaply, so that the changes with a
# file-contract component fail here rather than 50 times in parallel:
#   - PCA removal      (stage 3 dapg inputs, stage 5 glm covariates)
#   - fastenloc_prior  (stage 4 output contract drops from 5 files to 3)
#   - the Q=1 -> Q=0.01 cattle handoff (E1 only)
#
# Seeds follow submit_2Mb_rerun.sh's banding (category letter position x10 +
# replicate), so the pilot's A1/E1 are the same seeds the full grid will use and
# their outputs are reusable rather than throwaway:
#   A -> 11..20   B -> 21..30   E -> 51..60   F -> 61..70   G -> 71..80
#
# E1 requires the rescaled handoff to already exist at
#   $SCRATCH/cattle_baseline/stage1/farm_selection_Q_0.01.L_2000000.seed_20250303.ep7.ts
# (produced by helpers/rescale_checkpoint.py --from-Q 1 --to-Q 0.01). The check
# below is here because without it the failure surfaces only after the controller
# has queued, as a one-line stderr from inside the from-midpoint rule.
#
# Run once from a login node:  bash submit_2Mb_r3_pilot.sh

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_round_3_2Mb
PUBLISH=/n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake/simulations_round_3_2Mb
cd "$REPO"

HANDOFF="$SCRATCH/cattle_baseline/stage1/farm_selection_Q_0.01.L_2000000.seed_20250303.ep7.ts"
if [ ! -f "$HANDOFF" ]; then
    echo "ERROR: cattle handoff checkpoint missing: $HANDOFF" >&2
    echo "Run helpers/rescale_checkpoint.py on the Q=1 .ep7.ts first." >&2
    exit 1
fi

# submit_one PREFIX CONFIGFILE SEED [EXTRA_CONFIG]
submit_one() {
    local ID=$1 CFG=$2 SEED=$3 EXTRA="${4:-}"
    local WD="$SCRATCH/${ID}" PD="$PUBLISH/${ID}"
    mkdir -p "$WD"
    local JID
    JID=$(sbatch --parsable --job-name="${ID}" \
        --output="${WD}/controller_%j.out" --error="${WD}/controller_%j.err" \
        --export=ALL,CONFIGFILE="${CFG}",SEED="${SEED}",WD="${WD}",PD="${PD}",EXTRA_CONFIG="${EXTRA}" \
        controller_2Mb.sbatch)
    echo "  ${ID} (seed ${SEED}): controller $JID"
}

echo "== A1: human, directional-negative =="
submit_one A1 config/human_2Mb_r3.yaml 11
echo "== E1: cattle baseline (from Q=1 handoff) =="
submit_one E1 config/cattle_baseline_from_midpoint_2Mb_r3.yaml 51

echo
echo "Submitted 2 pilot controllers. Watch with:"
echo "  squeue -u \$USER -o '%.10i %.20j %.9P %.2t %.10M %R'"
echo "Then verify against the checklist before launching the full grid."
