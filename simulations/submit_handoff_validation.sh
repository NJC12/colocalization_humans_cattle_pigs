#!/bin/bash
# Checklist item 9: validate the split-Q cattle deep history against a pure
# single-Q run. Launches both arms across several seeds.
#
# The comparison is PAIRED across seeds, which is the point: a single pure-vs-
# hybrid pair cannot distinguish a handoff artifact from ordinary Monte Carlo
# noise. With N seeds we can ask whether the between-arm difference exceeds the
# between-seed spread within an arm, which is the question that actually matters.
#
# The hybrid arm is minutes (epochs 2-7 at Q=1 is ~3,330 ticks at 17,000-1,500
# individuals); the pure arm is hours (33,540 ticks at up to 170,000). They are
# submitted together and the pure arm sets the wall clock.
#
# See validate_handoff.sbatch for why QS=0.1 rather than 0.01.
#
# Run once from a login node:  bash submit_handoff_validation.sh
# Then, when everything has finished:
#   python helpers/compare_handoff.py --root <OUT> \
#       --out-tsv  figures_and_tables/burnin_diagnostics/handoff_validation.tsv \
#       --out-plot figures_and_tables/burnin_diagnostics/handoff_validation.png

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
OUT=${OUT:-/n/scratch/users/n/njc12/snakemake/simulations_round_3_2Mb/handoff_validation}
SEEDS=${SEEDS:-"901 902 903 904 905"}
QS=${QS:-0.1}
LEN=${LEN:-500000}
cd "$REPO"

mkdir -p "$OUT/logs"
echo "handoff validation -> $OUT   (QS=$QS, L=$LEN, seeds: $SEEDS)"
echo

for SEED in $SEEDS; do
    for ARM in hybrid pure; do
        JID=$(sbatch --parsable --job-name="hv_${ARM}_${SEED}" \
            --output="$OUT/logs/${ARM}_${SEED}_%j.out" \
            --error="$OUT/logs/${ARM}_${SEED}_%j.err" \
            --export=ALL,ARM="$ARM",SEED="$SEED",QS="$QS",LEN="$LEN",OUT="$OUT" \
            validate_handoff.sbatch)
        printf "  %-6s seed %s: job %s\n" "$ARM" "$SEED" "$JID"
    done
done

echo
echo "Watch with: squeue -u \$USER -o '%.10i %.18j %.9P %.2t %.10M %R'"
