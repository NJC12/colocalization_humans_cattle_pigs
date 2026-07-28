#!/bin/bash
#SBATCH --job-name=handoff_val
#SBATCH --partition=medium
#SBATCH --time=2-00:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=2
#
# Checklist item 9: does the split-Q architecture give the same answer as running
# everything at one Q?
#
# Two arms, run from independent msprime starting populations at the same seed:
#   ARM=pure    orig_pop(Q=$QS) -> farm_selection.slim epochs 2-12 entirely at Q=$QS
#   ARM=hybrid  orig_pop(Q=1)   -> farm_selection.slim epochs 2-7 at Q=1
#                               -> rescale_checkpoint.py --from-Q 1 --to-Q $QS
#                               -> farm_selection_from_ep8.slim epochs 8-12 at Q=$QS
# Both end at ep_sizes[12]/$QS individuals and are compared by
# helpers/compare_handoff.py on SFS, pi, Tajima's D and LD decay.
#
# WHY QS=0.1 AND NOT 0.01. A pure Q=0.01 arm is 335,400 ticks at 1.7M individuals
# -- roughly 39 days even at 500 kb, so the test would outlive the thing it
# validates. Q=0.1 is 33,540 ticks and runs in hours. The substitution is sound
# because the handoff artifact is set by the HYBRID's parental pool, which is
# ep_sizes[7] = 1500 individuals at Q_deep=1 in BOTH the real config and this
# test -- one extra WF generation at 1500 parents is 1/(2*1500) = 3.33e-4 of
# drift either way. What changes is only the pure arm it is measured against:
# 1/(2*150000) = 3.33e-6 at Q=0.01 versus 1/(2*15000) = 3.33e-5 at Q=0.1. Both
# are negligible against 3.33e-4, so this test recovers ~91% of the real
# configuration's handoff signal.
#
# Required via --export: ARM (pure|hybrid), SEED. Optional: QS, LEN, OUT.

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
PY=/home/njc12/miniconda3/envs/coloc_sims/bin/python
SLIM="bash /home/njc12/bin/slim_5/slim_wrapper.sh"

: "${ARM:?set ARM=pure or ARM=hybrid}" "${SEED:?set SEED}"
QS=${QS:-0.1}
LEN=${LEN:-500000}
RECOMB=${RECOMB:-1.3e-8}
OUT=${OUT:-/n/scratch/users/n/njc12/snakemake/simulations_round_3_2Mb/handoff_validation}

WD="$OUT/$ARM/seed_$SEED"
mkdir -p "$WD"
cd "$WD"

echo "=== handoff validation: ARM=$ARM SEED=$SEED QS=$QS LEN=$LEN ==="
date

if [ "$ARM" = "pure" ]; then
    ORIG="$WD/orig.Q_${QS}.ts"
    $PY "$REPO/farm_create_orig_pop_e2.py" \
        --seed "$SEED" --Q "$QS" --length "$LEN" --recomb_rate "$RECOMB" \
        --out_dir "$WD" --out_name "$(basename "$ORIG")"

    # tick=1 because we read the orig_pop, not a burn-in checkpoint.
    # farm_selection.slim schedules its epochs relative to sim.cycle after
    # readFromPopulationFile, so `tick` is informational only.
    $SLIM -m -t -l -s "$SEED" \
        -d Q_scaling="$QS" -d length="$LEN" -d recomb_rate="$RECOMB" -d tick=1 \
        -d file_in="\"$ORIG\"" -d dir_out="\"$WD/\"" -d file_out="\"sel\"" \
        "$REPO/farm_selection.slim"

elif [ "$ARM" = "hybrid" ]; then
    ORIG="$WD/orig.Q_1.ts"
    $PY "$REPO/farm_create_orig_pop_e2.py" \
        --seed "$SEED" --Q 1 --length "$LEN" --recomb_rate "$RECOMB" \
        --out_dir "$WD" --out_name "$(basename "$ORIG")"

    $SLIM -m -t -l -s "$SEED" \
        -d Q_scaling=1 -d length="$LEN" -d recomb_rate="$RECOMB" -d tick=1 \
        -d file_in="\"$ORIG\"" -d dir_out="\"$WD/\"" -d file_out="\"deep\"" \
        "$REPO/farm_selection.slim"

    $PY "$REPO/helpers/rescale_checkpoint.py" \
        --input "$WD/deep.ep7.ts" --output "$WD/deep.ep7.Q_${QS}.ts" \
        --from-Q 1 --to-Q "$QS"

    # Same seed as the pure arm on purpose: epochs 8-12 then draw the same
    # random stream in both arms, so a difference in the final samples is
    # attributable to the handoff rather than to Monte Carlo noise in the tail.
    $SLIM -m -t -l -s "$SEED" \
        -d Q_scaling="$QS" -d genome_length="$LEN" -d recomb_rate="$RECOMB" \
        -d num_muts_selected=0 \
        -d file_in="\"$WD/deep.ep7.Q_${QS}.ts\"" -d dir_out="\"$WD/\"" -d file_out="\"sel\"" \
        "$REPO/farm_selection_from_ep8.slim"
else
    echo "ERROR: ARM must be 'pure' or 'hybrid', got '$ARM'" >&2
    exit 1
fi

echo "=== done ==="
date
ls -l "$WD"/*.ts
