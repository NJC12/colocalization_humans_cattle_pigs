#!/bin/bash
# Launch the x20 phenotype-scaling variant of category E: 10 replicates.
#
#   E  cattle baseline (no artificial selection)  -> cattle_baseline_from_midpoint_2Mb_x20.yaml
#      (identical to cattle_baseline_from_midpoint_2Mb.yaml but
#       gwas_scaling=gtex_scaling=20, not 35)
#
# Seeds match the x35 E run exactly (E{N} -> seed 50+N -> 51..60) so each x20
# replicate REUSES the corresponding x35 stage-1 tree sequence
# (E{N}/stage1/farm_selection_from_ep8.*.seed_{seed}.full.ts) via
# stage1_search_dirs -- no genetic simulation reruns (the from-ep8 SLiM step is
# short-circuited to a symlink). Only stages 2-5 run.
#
# Outputs land under a separate ..._2Mb_x20 root so x35 and x20 sit side by side
# and never collide (stages 3-5 are named only by the stage-1 filename, which has
# no scaling in it, so a shared workdir would overwrite the x35 results).
#
# Nothing lands in `priority` or `long`: the controller is on `medium`, its child
# stage2 is forced to `medium`, and stages 3-5 (plus the cattle stage-1 symlink
# job) are `short` (see controller_2Mb_x10.sbatch).
#
# ASSUMES the x35 category-E run lives at $X35_SCRATCH/E{N}/ (the layout produced
# by submit_2Mb_rerun.sh). The per-replicate guard below verifies each reused
# stage-1 file exists before any job is submitted.
#
# Run once from a login node:  bash submit_2Mb_x20_E.sh

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake

# Source of the reused x35 stage-1 tree sequences (do NOT write here).
X35_SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb
# Destination roots for the x20 outputs.
XOUT_SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb_x20
XOUT_PUBLISH=/n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake/simulations_round_2_2Mb_x20
cd "$REPO"

NREP=10
BASE=50                 # E{N} -> seed BASE+N -> 51..60
CFG=config/cattle_baseline_from_midpoint_2Mb_x20.yaml

# Constants that appear in the reused stage-1 filename (must match the x35 config
# values -- see helpers/paths.py::stage1_cattle_baseline_from_midpoint_full).
Q=0.01
L=2000000
CB=20250303             # cattle_baseline_seed

echo "== E (x20): cattle baseline, gwas_scaling=gtex_scaling=20 =="
for N in $(seq 1 "$NREP"); do
    SEED=$((BASE + N))
    ID="E${N}"
    STAGE1_SRC="$X35_SCRATCH/${ID}/stage1"
    WD="$XOUT_SCRATCH/${ID}"
    PD="$XOUT_PUBLISH/${ID}"

    # Pre-submit guard: the reused x35 stage-1 tree sequence must already exist.
    # Without this, a wrong path would silently launch the from-ep8 SLiM step
    # instead of reusing the existing tree sequence.
    TS="$STAGE1_SRC/farm_selection_from_ep8.Q_${Q}.L_${L}.cb_${CB}.seed_${SEED}.full.ts"
    if [[ ! -f "$TS" ]]; then
        echo "ERROR: missing x35 stage-1 tree sequence for ${ID}: $TS" >&2
        echo "Aborting -- fix the reuse path before submitting (no jobs launched)." >&2
        exit 1
    fi

    mkdir -p "$WD"
    JID=$(sbatch --parsable --job-name="${ID}_x20" \
        --output="${WD}/controller_%j.out" --error="${WD}/controller_%j.err" \
        --export=ALL,CONFIGFILE="${CFG}",SEED="${SEED}",WD="${WD}",PD="${PD}",STAGE1_SRC="${STAGE1_SRC}",EXTRA_CONFIG="" \
        controller_2Mb_x10.sbatch)
    echo "  ${ID} (seed ${SEED}): controller $JID  [reuses $TS]"
done

echo "Submitted category E x20: ${NREP} replicates = ${NREP} controllers (none in priority/long)."
