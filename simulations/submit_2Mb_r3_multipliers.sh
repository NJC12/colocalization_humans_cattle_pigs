#!/bin/bash
# Round-3 phenotype-scaling variants of the A1 + E1 pilot: x10 and x20.
#
# The x35 pilot (submit_2Mb_r3_pilot.sh) established that cattle colocalizes more
# than human under round-3 code. These runs sweep the effect-size multiplier down
# to see how that gap behaves as the traits get weaker. Effect size is
#   beta = sqrt(|selection_coeff|) * scaling
# so x10 and x20 rescale every phenotype by 10/35 and 20/35 respectively while
# leaving the genetics completely untouched.
#
# Stage 1 is NOT rerun. Each variant reuses the x35 replicate's own stage-1 tree
# sequence via stage1_search_dirs, so the stage-1 rule collapses to an `ln -sf`
# and only stages 2-5 execute. That reuse is also what makes the comparison
# clean: x10, x20 and x35 share identical genetics and differ only in phenotype
# scaling.
#
# Outputs land under separate ..._2Mb_x10 / ..._2Mb_x20 roots. This is
# load-bearing, not tidiness: stages 3-5 are named only by the stage-1 filename,
# which carries no scaling, so a shared workdir would silently overwrite the x35
# results that are already on disk.
#
# Nothing lands in priority or long -- controller_2Mb_x10.sbatch runs on medium
# and forces its children off priority.
#
# Run once from a login node:  bash submit_2Mb_r3_multipliers.sh
# Override the set with e.g.:  MULTIPLIERS="10" bash submit_2Mb_r3_multipliers.sh

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
X35_SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_round_3_2Mb
cd "$REPO"

# Cells to run, as space-separated LABEL tokens. Each LABEL must have matching
# config/{human,cattle_baseline_from_midpoint}_2Mb_${LABEL}_r3.yaml files and
# gets its own ..._2Mb_${LABEL} output root.
#
#   x5 x10 x15 x20  symmetric: gwas_scaling == gtex_scaling
#   g20t5           asymmetric: gwas_scaling=20, gtex_scaling=5
#
# Asymmetric cells work without pipeline changes: helpers/paths.py selects
# gwas_scaling for *gwas categories and gtex_scaling for *gtex ones, and stage-2
# output is already namespaced gwas_{gw}_gtex_{gt}_maf_{maf}. Nothing assumes
# the two are equal.
CELLS="${CELLS:-${MULTIPLIERS:-10 20}}"

# Constants embedded in the reused cattle stage-1 filename (must match the x35
# config -- see helpers/paths.py::stage1_cattle_baseline_from_midpoint_full).
Q=0.01
L=2000000
CB=20250303             # cattle_baseline_seed

# id | seed | config-stem (LBL is replaced by the cell label) | reused stage-1 basename
REPLICATES=(
  "A1|11|human_2Mb_LBL_r3.yaml|hts_11.ts"
  "E1|51|cattle_baseline_from_midpoint_2Mb_LBL_r3.yaml|farm_selection_from_ep8.Q_${Q}.L_${L}.cb_${CB}.seed_51.full.ts"
)

# Verify every reused tree sequence and every config exists BEFORE submitting
# anything, so a bad path cannot leave a half-submitted set behind. Without the
# .ts guard a wrong path would silently launch a real genetic simulation instead
# of reusing the existing one.
for CELL in $CELLS; do
    for row in "${REPLICATES[@]}"; do
        IFS='|' read -r ID SEED CFGSTEM TSNAME <<< "$row"
        CFG="config/${CFGSTEM/LBL/$CELL}"
        TS="$X35_SCRATCH/${ID}/stage1/${TSNAME}"
        [[ -f "$REPO/$CFG" ]] || { echo "ERROR: missing config: $REPO/$CFG" >&2; exit 1; }
        [[ -f "$TS" ]] || {
            echo "ERROR: missing round-3 x35 stage-1 tree sequence for ${ID}: $TS" >&2
            echo "Run submit_2Mb_r3_pilot.sh first. Aborting -- no jobs launched." >&2
            exit 1; }
    done
done
echo "Pre-flight OK: all configs and reused stage-1 tree sequences present."
echo

for CELL in $CELLS; do
    OUT_SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_round_3_2Mb_${CELL}
    OUT_PUBLISH=/n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake/simulations_round_3_2Mb_${CELL}
    # Report the scalings from the config itself rather than parsing the label,
    # so the banner cannot drift from what actually runs.
    SC=$(grep -hE "^(gwas|gtex)_scaling:" "$REPO/config/human_2Mb_${CELL}_r3.yaml" \
         | awk '{printf "%s%s ", $1, $2}')
    echo "== ${CELL} (${SC}) =="
    for row in "${REPLICATES[@]}"; do
        IFS='|' read -r ID SEED CFGSTEM TSNAME <<< "$row"
        CFG="config/${CFGSTEM/LBL/$CELL}"
        STAGE1_SRC="$X35_SCRATCH/${ID}/stage1"
        WD="$OUT_SCRATCH/${ID}"
        PD="$OUT_PUBLISH/${ID}"
        mkdir -p "$WD"
        JID=$(sbatch --parsable --job-name="${ID}_${CELL}" \
            --output="${WD}/controller_%j.out" --error="${WD}/controller_%j.err" \
            --export=ALL,CONFIGFILE="${CFG}",SEED="${SEED}",WD="${WD}",PD="${PD}",STAGE1_SRC="${STAGE1_SRC}",EXTRA_CONFIG="" \
            controller_2Mb_x10.sbatch)
        printf "  %-3s %-6s (seed %s): controller %s  [reuses %s]\n" "$ID" "$CELL" "$SEED" "$JID" "$TSNAME"
    done
done

echo
echo "Watch with: squeue -u \$USER -o '%.10i %.20j %.9P %.2t %.10M %R'"
