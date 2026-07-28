#!/bin/bash
# Re-launch ONLY the cattle categories (E/F/G) after the SLiM 5 .genomes ->
# .haplosomes fix in farm_selection_from_ep8.slim. The human A/B controllers
# from the first launch are still running fine and are NOT touched here.
#
#   E -> 51..60   F -> 61..70   G -> 71..80
#
# Their workdirs (E*/F*/G*) should be wiped before this runs (the first launch
# left failed stage-1 logs + per-workdir .snakemake). Run once from a login node:
#   bash submit_2Mb_cattle_rerun.sh
# Nothing lands in priority.

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb
PUBLISH=/n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake/simulations_round_2_2Mb
cd "$REPO"

NREP=10

# submit_cat PREFIX CONFIGFILE SEED_BASE
submit_cat() {
    local PREFIX=$1 CFG=$2 BASE=$3
    for N in $(seq 1 "$NREP"); do
        local SEED=$((BASE + N)) ID="${PREFIX}${N}"
        local WD="$SCRATCH/${ID}" PD="$PUBLISH/${ID}"
        mkdir -p "$WD"
        local JID
        JID=$(sbatch --parsable --job-name="${ID}" \
            --output="${WD}/controller_%j.out" --error="${WD}/controller_%j.err" \
            --export=ALL,CONFIGFILE="${CFG}",SEED="${SEED}",WD="${WD}",PD="${PD}",EXTRA_CONFIG="" \
            controller_2Mb.sbatch)
        echo "  ${ID} (seed ${SEED}): controller $JID"
    done
}

echo "== E: cattle baseline =="
submit_cat E config/cattle_baseline_from_midpoint_2Mb.yaml 50
echo "== F: cattle + positive selection, bottlenecked =="
submit_cat F config/cattle_sel_bottlenecked_2Mb.yaml 60
echo "== G: cattle + positive selection, not bottlenecked =="
submit_cat G config/cattle_sel_not_bottlenecked_2Mb.yaml 70

echo "Submitted 3 cattle categories x ${NREP} = $((3 * NREP)) controllers (none in priority)."
