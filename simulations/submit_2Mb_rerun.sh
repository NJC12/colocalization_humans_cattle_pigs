#!/bin/bash
# Launch the full 2 Mb study: 5 categories x 10 replicates = 50 controllers.
#
# Categories (see README_snakemake.md):
#   A  human, directional-negative trait selection      -> human_2Mb.yaml
#   B  human, neutral trait variants                     -> human_2Mb.yaml + neutral_trait_vars=True
#   E  cattle baseline (no artificial selection)         -> cattle_baseline_from_midpoint_2Mb.yaml
#   F  cattle + positive selection, WITH bottlenecking   -> cattle_sel_bottlenecked_2Mb.yaml
#   G  cattle + positive selection, NO  bottlenecking    -> cattle_sel_not_bottlenecked_2Mb.yaml
#
# Seeds are banded by category (alphabet position x10, ones digit = replicate):
#   A -> 11..20   B -> 21..30   E -> 51..60   F -> 61..70   G -> 71..80
#
# Every replicate runs under the new trait scheme (50 central GWAS loci = 50
# shared GTEx loci, + 50 GTEx-only flank loci -> 100 GTEx), with GTEx fanned out
# to 1000 + 500 donors (gtex_size=-1). The cattle categories (E/F/G) all continue
# from the SAME shared 2 Mb epoch_8 checkpoint (no extra baseline run).
#
# Run once from a login node:  bash submit_2Mb_rerun.sh
# Nothing lands in priority (controller on `long`; child stage2 -> medium,
# human stage1 -> long via --set-resources in controller_2Mb.sbatch).

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
echo "== G: cattle + positive selection, not bottlenecked =="
submit_cat G config/cattle_sel_not_bottlenecked_2Mb.yaml 70

echo "Submitted 5 categories x ${NREP} replicates = $((5 * NREP)) controllers (none in priority)."
