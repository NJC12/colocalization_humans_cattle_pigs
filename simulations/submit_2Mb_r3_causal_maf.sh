#!/bin/bash
# Causal-MAF-floor arms of the round-3 A1 + E1 grid.
#
# Every run so far drew its causative variants from MAF >= 0.01 -- the same floor
# used for the GWAS and (once fm_min_maf arrived) for fine-mapping. This script
# lowers ONLY the causative floor, leaving the tested set alone:
#
#   causal_min_maf=0.001   causal variants drawn from MAF >= 0.001
#   fm_min_maf=0.01        but only MAF >= 0.01 variants are fine-mapped
#   min_maf=0.01           and only MAF >= 0.01 variants enter the GWAS
#
# The floor has NO ceiling, so this does not make the run unfindable -- most of the
# MAF >= 0.001 pool is still above 0.01. It splits the loci in two:
#
#   causal MAF >= 0.01        in both tested sets; behaves exactly as at baseline
#   causal MAF in [.001,.01)  in NEITHER; findable only through a variant that tags it
#
# The second group is the question these arms exist to answer, and it is why
# fm_min_maf must stay at 0.01 -- drop it to 0 and you are instead asking "what
# happens when we fine-map rare variants", which the maf01 arms of
# submit_2Mb_r3_finemap_variants.sh already covered from the other direction.
#
# SIZING: n_central_traits is drawn from the LOWERED pool, so only
# n_central_traits * P(MAF >= 0.01 | MAF >= 0.001) loci land in the findable group.
# That fraction is much smaller for human than for bottlenecked cattle (88% vs 9%
# of variants below MAF 0.01 per +/-250 kb window), so at the usual 50 the human
# arm's findable subset is the smallest number in the comparison that matters.
# N_CENTRAL below raises it; the cost is more DAP-G loci per replicate.
#
# THREE DIFFERENCES FROM ITS SIBLING SCRIPT, all deliberate:
#
#  1. stage2_search_dirs is NOT set, and must not be. The finemap-variant script
#     reuses stage 2 wholesale because ld_ctrl and fm_min_maf do not change the
#     phenotypes. causal_min_maf does -- it decides which variants become traits at
#     all -- so stage 2 has to rerun. The Snakefile enforces this: the stage-2
#     directory is named for the causal floor, so a baseline directory does not
#     match, and its provenance guard refuses a stage-2 dir whose recorded floor
#     disagrees.
#  2. stage1_search_dirs still points at the replicate's own stage-1 tree
#     sequence, so the genetics are shared with the baseline cell. Without that
#     the arms are not comparable: any difference could be demography rather than
#     the MAF floor.
#  3. Output roots are suffixed _cmaf_<floor>, matching the path segment
#     helpers/paths.py generates, so the two arms coexist even if someone points
#     them at one root.
#
# Default is deliberately narrow -- one cell, one floor, two arms = 2 controllers.
# The finemap script's 6-cell default exists because that was a planned complete
# sweep; this axis is new, so run Phase 3's single small run first, look at the
# pool sizes in logs/stage2_split_pheno.log, then widen.
#
# Run from a login node:
#   bash submit_2Mb_r3_causal_maf.sh                          # x35, floor 0.001
#   FLOORS="0.001 0.005" bash submit_2Mb_r3_causal_maf.sh     # two floors
#   CELLS="x35 g5t20" FLOORS="0.001" bash submit_2Mb_r3_causal_maf.sh

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SCRATCH=/n/scratch/users/n/njc12/snakemake
PUBLISH=/n/data2/hms/dbmi/sunyaev/lab/nconnally/snakemake
cd "$REPO"

CELLS="${CELLS:-x35}"
FLOORS="${FLOORS:-0.001}"
# The fine-mapping floor to hold the tested set at. Keep it at or above the
# highest FLOORS value or the experiment collapses into an ordinary rare-variant
# run -- see the header.
FM_FLOOR="${FM_FLOOR:-0.01}"
# Central trait loci. Empty = leave the config's value (50) alone. Set it higher to
# compensate for the share of loci that land below FM_FLOOR and are therefore
# unfindable -- see SIZING in the header.
N_CENTRAL="${N_CENTRAL:-}"

# Every cell's stage-1 tree sequence lives in the x35 root and nowhere else (the
# multiplier runs were launched with stage1_search_dirs pointing here), so stage 1
# is always sourced from there regardless of which cell is being rerun.
STAGE1_ROOT="$SCRATCH/simulations_round_3_2Mb"

# Optional SLURM --begin for the controllers, e.g. BEGIN=now+90minutes. Each
# controller runs Snakemake with `jobs: 200`, so keep concurrent waves modest
# against the 10000 MaxSubmit association limit.
BEGIN="${BEGIN:-}"

# id | seed | config stem; full path is config/<stem><cell suffix>_r3.yaml
REPLICATES=(
  "A1|11|human_2Mb"
  "E1|51|cattle_baseline_from_midpoint_2Mb"
)

# x35 is the original pilot: its configs and output root carry no cell label.
cell_suffix() {
    case "$1" in
        x35) echo "" ;;
        *)   echo "_$1" ;;
    esac
}

# ---------------------------------------------------------------- pre-flight
# Validate every config and every reused stage-1 tree sequence BEFORE submitting.
# Without the stage-1 guard a wrong path silently launches a real genetic
# simulation -- days of compute -- instead of reusing one.
for FLOOR in $FLOORS; do
    awk -v f="$FLOOR" 'BEGIN{ exit !(f+0 > 0 && f+0 < 0.5) }' </dev/null || {
        echo "ERROR: FLOORS entry '$FLOOR' is not in (0, 0.5)" >&2; exit 1; }
    awk -v f="$FLOOR" -v m="$FM_FLOOR" 'BEGIN{ exit !(f+0 <= m+0) }' </dev/null || {
        echo "WARNING: causal floor $FLOOR is ABOVE fm_min_maf $FM_FLOOR, so every" >&2
        echo "         causal variant is still fine-mapped. That is a valid run, but" >&2
        echo "         it is not the untested-causal-variant experiment." >&2; }
done

for CELL in $CELLS; do
    SFX=$(cell_suffix "$CELL")
    for row in "${REPLICATES[@]}"; do
        IFS='|' read -r ID SEED CFGSTEM <<< "$row"
        CFG="config/${CFGSTEM}${SFX}_r3.yaml"
        [[ -f "$REPO/$CFG" ]] || { echo "ERROR: missing config: $REPO/$CFG" >&2; exit 1; }

        S1_DIR="$STAGE1_ROOT/${ID}/stage1"
        compgen -G "$S1_DIR/*.ts" > /dev/null || {
            echo "ERROR: no stage-1 tree sequence for ${ID} in $S1_DIR" >&2
            echo "Run submit_2Mb_r3_pilot.sh first. Aborting -- no jobs launched." >&2
            exit 1; }
    done
done
N_EXPECT=$(( $(echo $CELLS | wc -w) * $(echo $FLOORS | wc -w) * ${#REPLICATES[@]} ))
echo "Pre-flight OK: $(echo $CELLS | wc -w) cells x $(echo $FLOORS | wc -w) floors x ${#REPLICATES[@]} arms = ${N_EXPECT} controllers."
echo "Stage 2 WILL rerun (that is the point); stage 1 is reused from $STAGE1_ROOT."
echo

# ------------------------------------------------------------------- submit
N=0
for FLOOR in $FLOORS; do
    echo "=== causal_min_maf=${FLOOR}  (fm_min_maf=${FM_FLOOR}) ==="
    for CELL in $CELLS; do
        SFX=$(cell_suffix "$CELL")
        OUT_SCRATCH="$SCRATCH/simulations_round_3_2Mb${SFX}_cmaf_${FLOOR}"
        OUT_PUBLISH="$PUBLISH/simulations_round_3_2Mb${SFX}_cmaf_${FLOOR}"
        for row in "${REPLICATES[@]}"; do
            IFS='|' read -r ID SEED CFGSTEM <<< "$row"
            CFG="config/${CFGSTEM}${SFX}_r3.yaml"
            WD="$OUT_SCRATCH/${ID}"
            PD="$OUT_PUBLISH/${ID}"
            mkdir -p "$WD"
            # EXTRA_CONFIG is expanded unquoted by the controller so each
            # space-separated token becomes its own --config arg. No token may
            # contain a space. Note the absence of stage2_search_dirs.
            EXTRA="causal_min_maf=${FLOOR} fm_min_maf=${FM_FLOOR}"
            [[ -n "$N_CENTRAL" ]] && EXTRA="${EXTRA} n_central_traits=${N_CENTRAL}"
            JID=$(sbatch --parsable --job-name="${ID}_${CELL}_cmaf${FLOOR}" \
                ${BEGIN:+--begin="$BEGIN"} \
                --output="${WD}/controller_%j.out" --error="${WD}/controller_%j.err" \
                --export=ALL,CONFIGFILE="${CFG}",SEED="${SEED}",WD="${WD}",PD="${PD}",STAGE1_SRC="${STAGE1_ROOT}/${ID}/stage1",EXTRA_CONFIG="${EXTRA}" \
                controller_2Mb_x10.sbatch)
            printf "  %-3s %-6s cmaf=%-7s controller %s\n" "$ID" "$CELL" "$FLOOR" "$JID"
            N=$((N + 1))
        done
    done
    echo
done
echo "Submitted ${N} controllers."
echo "Watch with: squeue -u \$USER -o '%.10i %.24j %.9P %.2t %.10M %R'"
echo
echo "First thing to check once stage 2 lands, per replicate:"
echo "  grep -E 'causative|neutral|flank' \$WD/logs/stage2_split_pheno.log"
echo "The causative-pool count must exceed the baseline cell's. If it matches,"
echo "the floor did not take effect and the arms are not testing anything."
