#!/bin/bash
# One-shot submission for the 2 Mb re-run. Run ONCE from an O2 login node:
#   cd /n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
#   bash submit_2Mb.sh
# This only fires sbatch jobs and exits -- it does NOT stay resident, so there
# is no login-node tmux/screen to be lost on a node restart. Every unit of work
# runs as a SLURM job.
#
# PRECONDITION: the extracted 2 Mb burn-in must already be in place at
#   .../simulations_round_2_2Mb/cattle_baseline/stage1/farm_burn_in.Q_0.01.L_2000000.seed_20250303.cycle_25000.ts
# (produced by helpers/extract_subregion.py from the 10 Mb cycle_25000.ts).
#
# Seeds: cattle replicates E1..E8 -> 51..58 ; human replicates H1..H8 -> 61..68.

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb
PUBLISH=/n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake/simulations_round_2_2Mb
cd "$REPO"

# 1) Shared cattle baseline selection (standalone, long; produces epoch_8.ts).
JOB_BASE=$(sbatch --parsable submit_cattle_baseline_2Mb_selection.sbatch)
echo "cattle baseline 2Mb selection: job $JOB_BASE"

# 2) 8 cattle from-midpoint replicate controllers (start after the baseline).
for N in 1 2 3 4 5 6 7 8; do
    SEED=$((50 + N))
    WD="$SCRATCH/cattle_E${N}"; PD="$PUBLISH/cattle_E${N}"; mkdir -p "$WD"
    JID=$(sbatch --parsable --job-name="cattle_E${N}" \
        --dependency=afterok:${JOB_BASE} \
        --output="${WD}/controller_%j.out" --error="${WD}/controller_%j.err" \
        --export=ALL,CONFIGFILE=config/cattle_baseline_from_midpoint_2Mb.yaml,SEED=${SEED},WD="${WD}",PD="${PD}" \
        controller_2Mb.sbatch)
    echo "  cattle E${N} (seed ${SEED}): controller job $JID  [after $JOB_BASE]"
done

# 3) 8 human replicate controllers (independent; fresh runs).
for N in 1 2 3 4 5 6 7 8; do
    SEED=$((60 + N))
    WD="$SCRATCH/human_E${N}"; PD="$PUBLISH/human_E${N}"; mkdir -p "$WD"
    JID=$(sbatch --parsable --job-name="human_E${N}" \
        --output="${WD}/controller_%j.out" --error="${WD}/controller_%j.err" \
        --export=ALL,CONFIGFILE=config/human_2Mb.yaml,SEED=${SEED},WD="${WD}",PD="${PD}" \
        controller_2Mb.sbatch)
    echo "  human H${N} (seed ${SEED}): controller job $JID"
done

echo "Submitted: 1 baseline selection + 8 cattle replicate controllers + 8 human replicate controllers."
