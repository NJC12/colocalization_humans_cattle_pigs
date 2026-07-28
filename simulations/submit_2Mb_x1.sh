#!/bin/bash
# Launch the x1 phenotype-scaling variant of category A: 10 replicates.
#
#   A  human, directional-negative trait selection  -> human_2Mb_x1.yaml
#      (identical to human_2Mb.yaml but gwas_scaling=gtex_scaling=1, not 35)
#
# Seeds match the x35 A run exactly (A{N} -> seed 10+N -> 11..20) so each x1
# replicate REUSES the corresponding x35 stage-1 tree sequence
# (A{N}/stage1/hts_{seed}.ts) via stage1_search_dirs -- no genetic simulation
# reruns. Only stages 2-5 (trait creation, DAP-G, fastEnloc, plink GLM) run.
#
# Outputs land under a separate ..._2Mb_x1 root so x35 and x1 sit side by side
# and never collide (stages 3-5 are named only by the stage-1 filename, which has
# no scaling in it, so a shared workdir would overwrite the x35 results).
#
# Nothing lands in `priority` or `long`: the controller is on `medium`, its child
# stage2 is forced to `medium`, and stages 3-5 are `short` (see
# controller_2Mb_x10.sbatch). No stage1_human job runs (stage 1 is reused).
#
# Run once from a login node:  bash submit_2Mb_x1.sh

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake

# Source of the reused x35 stage-1 tree sequences (do NOT write here).
X35_SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb
# Destination roots for the x1 outputs.
XOUT_SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb_x1
XOUT_PUBLISH=/n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake/simulations_round_2_2Mb_x1
cd "$REPO"

NREP=10
BASE=10                 # A{N} -> seed BASE+N -> 11..20
CFG=config/human_2Mb_x1.yaml

echo "== A (x1): human, directional-negative, gwas_scaling=gtex_scaling=1 =="
for N in $(seq 1 "$NREP"); do
    SEED=$((BASE + N))
    ID="A${N}"
    STAGE1_SRC="$X35_SCRATCH/${ID}/stage1"
    WD="$XOUT_SCRATCH/${ID}"
    PD="$XOUT_PUBLISH/${ID}"

    # Pre-submit guard: the reused x35 stage-1 tree sequence must already exist.
    # Without this, a wrong path would silently launch a multi-day stage1_human
    # genetic simulation instead of reusing the existing one.
    TS="$STAGE1_SRC/hts_${SEED}.ts"
    if [[ ! -f "$TS" ]]; then
        echo "ERROR: missing x35 stage-1 tree sequence for ${ID}: $TS" >&2
        echo "Aborting -- fix the reuse path before submitting (no jobs launched)." >&2
        exit 1
    fi

    mkdir -p "$WD"
    JID=$(sbatch --parsable --job-name="${ID}_x1" \
        --output="${WD}/controller_%j.out" --error="${WD}/controller_%j.err" \
        --export=ALL,CONFIGFILE="${CFG}",SEED="${SEED}",WD="${WD}",PD="${PD}",STAGE1_SRC="${STAGE1_SRC}",EXTRA_CONFIG="" \
        controller_2Mb_x10.sbatch)
    echo "  ${ID} (seed ${SEED}): controller $JID  [reuses $TS]"
done

echo "Submitted category A x1: ${NREP} replicates = ${NREP} controllers (none in priority/long)."
