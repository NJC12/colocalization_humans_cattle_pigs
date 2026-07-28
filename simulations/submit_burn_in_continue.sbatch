#!/bin/bash
#SBATCH --job-name=cattle_burnin_hedge
#SBATCH --partition=long
#SBATCH --time=30-00:00:00
#SBATCH --mem=150G
#SBATCH --cpus-per-task=4
#SBATCH --output=/n/scratch/users/n/njc12/snakemake/simulations_round_3_2Mb/cattle_burnin_q001/logs/burn_in_continue.slurm_%j.out
#SBATCH --error=/n/scratch/users/n/njc12/snakemake/simulations_round_3_2Mb/cattle_burnin_q001/logs/burn_in_continue.slurm_%j.err
#
# Round-3 Q=0.01 burn-in HEDGE (Track 0).
#
# The primary round-3 cattle deep history runs at Q_scaling=1 (see the plan);
# this job exists so the Q=0.01 architecture is available to swap in later
# without waiting a month for it. Nothing downstream depends on it. To adopt it,
# point `cattle_baseline_search_dirs` in the from-midpoint config at whatever
# this produces instead of the Q=1 handoff checkpoint.
#
# 2 Mb rather than 10 Mb on purpose: ~5x cheaper per tick, so ~5x more burn-in
# depth per wall-clock, and it matches the human arm, which is simulated
# standalone at 2 Mb rather than carved out of a larger region.
#
# CHAINING: this is one 30-day link. farm_burn_in_continue.slim writes
# `.cycle_<N>.ts` every 1000 ticks (10/Q), so the next link just takes the
# highest checkpoint as its FILE_IN and a larger END_TICK. Override either from
# the command line:
#     sbatch --export=ALL,FILE_IN=<...cycle_130000.ts>,END_TICK=260000 \
#            submit_burn_in_continue.sbatch
#
# PRECONDITION: the 2 Mb starting population must exist at FILE_IN. For the
# first link that is either
#   (a) helpers/extract_subregion.py --end 2000000 applied to the existing
#       10 Mb farm_orig_pop.Q_0.01.L_10000000.seed_20250303.ts -- still valid
#       at the round-3 rate, because farm_create_orig_pop_e2.py already used
#       1.3e-8; or
#   (b) a fresh farm_create_orig_pop_e2.py --length 2000000 --recomb_rate 1.3e-8
#       run if that file has been purged from scratch.
# Starting from the ROUND-2 `.cycle_25000.ts` burn-in checkpoint is NOT valid:
# it was produced at 9.26e-9.

set -euo pipefail

ROOT=${ROOT:-/n/scratch/users/n/njc12/snakemake/simulations_round_3_2Mb/cattle_burnin_q001}
STAGE1="$ROOT/stage1"
LOGDIR="$ROOT/logs"
mkdir -p "$STAGE1" "$LOGDIR"

SEED=${SEED:-20250303}
Q=${Q:-0.01}
L=${L:-2000000}
RECOMB=${RECOMB:-1.3e-8}
OUT_BASE=${OUT_BASE:-farm_burn_in.Q_${Q}.L_${L}.seed_${SEED}}

# Absolute tick to stop at. ~20 s/tick at 2 Mb / 1.7M individuals puts roughly
# 130k ticks (~1300 real generations at Q=0.01) inside a 30-day wall.
END_TICK=${END_TICK:-130000}

# Checkpoint every 5000 ticks (50 real generations) rather than the script's
# 10/Q default of 1000. Each checkpoint carries all 1.7M individuals and 3.4M
# nodes, so they are hundreds of MB apiece; 1000-tick spacing would write 130 of
# them. 5000 gives ~26 time points -- the same number the round-2 convergence
# analysis used -- for a fifth of the scratch.
CKPT_EVERY=${CKPT_EVERY:-5000}
FILE_IN=${FILE_IN:-$STAGE1/farm_orig_pop.Q_${Q}.L_${L}.seed_${SEED}.ts}

# FILE_IN=AUTO resolves to the highest-numbered .cycle_<N>.ts in $STAGE1, at JOB
# START rather than at submit time. That is what makes this script chainable
# under `sbatch --dependency=afterany:<prev>`: when the link is queued, the
# checkpoint it will resume from does not exist yet and its tick number is not
# knowable, so a literal FILE_IN cannot be supplied in advance.
if [ "$FILE_IN" = "AUTO" ]; then
    FILE_IN=$(ls -1 "$STAGE1"/${OUT_BASE}.cycle_*.ts 2>/dev/null \
              | sed 's/.*\.cycle_\([0-9]*\)\.ts/\1 &/' | sort -n -k1,1 | tail -1 | cut -d' ' -f2-)
    if [ -z "$FILE_IN" ]; then
        echo "ERROR: FILE_IN=AUTO found no ${OUT_BASE}.cycle_*.ts in $STAGE1" >&2
        echo "The previous link produced no checkpoint; nothing to resume from." >&2
        exit 1
    fi
    echo "FILE_IN=AUTO resolved to: $FILE_IN"
    # A previous link that already reached END_TICK leaves a checkpoint at or past
    # it. farm_burn_in_continue.slim would stop() with a non-obvious message, so
    # say plainly that the chain is complete and exit clean rather than failing.
    GOT=$(basename "$FILE_IN" | sed 's/.*\.cycle_\([0-9]*\)\.ts/\1/')
    if [ "$GOT" -ge "$END_TICK" ]; then
        echo "Nothing to do: latest checkpoint is tick $GOT, END_TICK is $END_TICK."
        echo "Raise END_TICK to extend the chain further."
        exit 0
    fi
fi

if [ ! -f "$FILE_IN" ]; then
    echo "ERROR: FILE_IN does not exist: $FILE_IN" >&2
    echo "See the PRECONDITION note at the top of this script." >&2
    exit 1
fi

echo "burn-in hedge: $FILE_IN -> tick $END_TICK (Q=$Q, L=$L, r=$RECOMB)"

bash /home/njc12/bin/slim_5/slim_wrapper.sh -m -t -l \
  -s "$SEED" \
  -d Q_scaling="$Q" \
  -d length="$L" \
  -d recomb_rate="$RECOMB" \
  -d end_tick="$END_TICK" \
  -d checkpoint_every="$CKPT_EVERY" \
  -d file_in='"'"$FILE_IN"'"' \
  -d dir_out='"'"$STAGE1/"'"' \
  -d file_out='"'"$OUT_BASE"'"' \
  /n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake/farm_burn_in_continue.slim \
  > "$LOGDIR/burn_in_continue.log" 2>&1
