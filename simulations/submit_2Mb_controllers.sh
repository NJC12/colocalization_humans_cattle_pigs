#!/bin/bash
# Submit ONLY the 16 replicate controllers for the 2 Mb run. The cattle baseline
# 2 Mb selection is already running and is NOT resubmitted (its job id is passed
# in / defaulted below). Run once from a login node:
#   bash submit_2Mb_controllers.sh           # uses the default baseline job id
#   bash submit_2Mb_controllers.sh <jobid>   # or pass the baseline job id
# Uses the fixed controller_2Mb.sbatch (partition=long; child rules stage2 ->
# medium, human stage1 -> long via --set-resources) so NOTHING here lands in
# priority -- priority is reserved for the cattle baseline selection only.

set -euo pipefail
JOB_BASE=${1:-45511685}   # already-running cattle baseline 2 Mb selection
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb
PUBLISH=/n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake/simulations_round_2_2Mb
cd "$REPO"
echo "Cattle replicate controllers will start after baseline job $JOB_BASE."

# 8 cattle from-midpoint replicate controllers (start after the running baseline).
for N in 1 2 3 4 5 6 7 8; do
    SEED=$((50 + N)); WD="$SCRATCH/cattle_E${N}"; PD="$PUBLISH/cattle_E${N}"; mkdir -p "$WD"
    JID=$(sbatch --parsable --job-name="cattle_E${N}" \
        --dependency=afterok:${JOB_BASE} \
        --output="${WD}/controller_%j.out" --error="${WD}/controller_%j.err" \
        --export=ALL,CONFIGFILE=config/cattle_baseline_from_midpoint_2Mb.yaml,SEED=${SEED},WD=${WD},PD=${PD} \
        controller_2Mb.sbatch)
    echo "  cattle E${N} (seed ${SEED}): controller $JID  [after $JOB_BASE]"
done

# 8 human replicate controllers (independent; fresh runs).
for N in 1 2 3 4 5 6 7 8; do
    SEED=$((60 + N)); WD="$SCRATCH/human_E${N}"; PD="$PUBLISH/human_E${N}"; mkdir -p "$WD"
    JID=$(sbatch --parsable --job-name="human_E${N}" \
        --output="${WD}/controller_%j.out" --error="${WD}/controller_%j.err" \
        --export=ALL,CONFIGFILE=config/human_2Mb.yaml,SEED=${SEED},WD=${WD},PD=${PD} \
        controller_2Mb.sbatch)
    echo "  human H${N} (seed ${SEED}): controller $JID"
done

echo "Submitted 8 cattle + 8 human replicate controllers (none in priority)."
