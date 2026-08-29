#!/bin/bash
#
# Rebuild the cattle deep history from seeds: orig_pop -> burn-in -> epochs 1-11,
# then rescale the epoch-7 checkpoint from Q=1 to Q=0.01 for the replicates to
# resume from.
#
#   bash rebuild_cattle_deep_history.sh            # submit
#   bash rebuild_cattle_deep_history.sh --rescale  # after it finishes
#   bash rebuild_cattle_deep_history.sh --check    # compare against the archive
#
# WHY REBUILD AT ALL. Every cattle replicate (E, F, G, L) resumes from one shared
# farm_selection_*.ep7.ts. The chain that produced the archived copy was seeded
# throughout EXCEPT for one line -- farm_create_orig_pop_e2.py's
# `np.random.permutation`, which decides which mutation receives which selection
# coefficient. So the archived checkpoint could not be regenerated from its seed,
# and neither could anything downstream of it.
#
# That line is now `rng.permutation` against the module's seeded generator, which
# makes the whole chain reproducible -- and rebuilding is cheap. Measured:
#
#   farm_create_orig_pop_e2.py    13-16 min   (sacct cattle_origpop_q001)
#   burn-in, 30,000 ticks @ Q=1   38m 48s     (sacct 48879426 cattle_burnin_q1)
#   epochs 1-11 @ Q=1             ~40 min     (epoch 1 is 29,800 ticks at the
#                                              same 17,000 individuals as the
#                                              burn-in's 30,000; epochs 2-11 are
#                                              3,354 ticks at 10,000->90 and
#                                              their checkpoints land 23s apart)
#   rescale Q=1 -> Q=0.01         seconds
#
# i.e. about two hours, not the ~60 days the 30-day walltime declarations in
# rules/stage1_cattle_baseline.smk suggest. Those exist for the Q=0.01 case --
# the `cattle_burnin_hedge` run, which has consumed ~32 days to reach 6.6% of a
# 2,500,000-tick target and is a different simulation entirely.
#
# It writes to a NEW workdir. The archived round-3 deep history is left exactly
# where it is, so the old results stay reproducible from it and the two can be
# compared (--check).

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${REPO:-$HERE}"
source "$REPO/lib/publication_common.sh"

SNAKEMAKE="${SNAKEMAKE:-/home/njc12/miniconda3/envs/coloc_sims/bin/snakemake}"
PYTHON="${PYTHON:-/home/njc12/miniconda3/envs/coloc_sims/bin/python}"
CFG="${CFG:-config/cattle_baseline_2Mb_r3.yaml}"
# The archived deep history, for --check only. Never written to.
ARCHIVE_WD="${ARCHIVE_WD:-/n/scratch/users/n/njc12/snakemake/simulations_round_3_2Mb/cattle_baseline}"

# ACCEPTED, not just derived. The seed used to be read out of the config and used
# only to spell filenames for --rescale/--check -- it was never passed to
# Snakemake at all (the --config list below carried workdir and publishdir and
# nothing else), so this script could build exactly one deep history: the one its
# configfile named. Four independent histories need four seeds.
CB_SEED="${CB_SEED:-$(awk '/^stage1_seed:/{print $2; exit}' "$REPO/$CFG")}"

# Per-seed workdir. The tree-sequence FILENAMES already carry the seed, but the
# workdir does not, and everything else under it is shared: .snakemake/locks,
# .snakemake/iocache, slurm_logs, params/, and the three rules' log: paths
# (rules/stage1_cattle_baseline.smk:86,152,210) which are workdir-relative and
# NOT seed-namespaced. Concurrent builds in one workdir would collide on all of
# them. Same shape as validate_handoff.sbatch:43, which already builds five
# independent histories this way.
WD="${WD:-$SCRATCH_ROOT/deep_history/cattle_baseline_seed_$CB_SEED}"
Q_DEEP=$(awk '/^Q_scaling:/{print $2; exit}' "$REPO/$CFG")
L_VAL=$(awk '/^L:/{print $2; exit}' "$REPO/$CFG")
EP7="farm_selection_Q_${Q_DEEP}.L_${L_VAL}.seed_${CB_SEED}.ep7.ts"
EP7_RESCALED="farm_selection_Q_0.01.L_${L_VAL}.seed_${CB_SEED}.ep7.ts"

if [[ "${1:-}" == "--rescale" ]]; then
    IN="$WD/stage1/$EP7"
    OUT="$WD/stage1/$EP7_RESCALED"
    [[ -f "$IN" ]] || { echo "ERROR: $IN not there yet -- has the DAG finished?" >&2; exit 1; }
    # Q=1 deep history -> Q=0.01 epochs 8-12. rescale_checkpoint multiplies the
    # selection coefficients by to_Q/from_Q and rescales the time axis, which is
    # what makes the piecewise time scale (handoff_ticks / deep_Q_scaling in the
    # from-midpoint configs) line up.
    "$PYTHON" "$REPO/helpers/rescale_checkpoint.py" \
        --input "$IN" --output "$OUT" --from-Q "$Q_DEEP" --to-Q 0.01 || exit 1
    echo "wrote $OUT"
    echo
    echo "Next: copy it into the publication stage-1 inputs, which is what every"
    echo "cattle replicate's cattle_baseline_search_dirs will point at:"
    echo "  mkdir -p $PUBLISH_ROOT/stage1_inputs"
    echo "  cp -p $OUT $PUBLISH_ROOT/stage1_inputs/"
    exit 0
fi

if [[ "${1:-}" == "--check" ]]; then
    echo "Fresh vs archived deep history (both are legitimate; they differ because"
    echo "the archived one was built with the UNSEEDED permutation):"
    for f in "farm_orig_pop.Q_${Q_DEEP}.L_${L_VAL}.seed_${CB_SEED}.ts" "$EP7"; do
        a="$WD/stage1/$f"; b="$ARCHIVE_WD/stage1/$f"
        if [[ -f "$a" && -f "$b" ]]; then
            "$PYTHON" "$REPO/helpers/compare_tree_sequences.py" "$a" "$b" 2>&1 | tail -1 \
                | awk -F'\t' -v n="$f" '{printf "  %-56s %s  %s\n", n, $1, $4}'
        else
            printf "  %-56s %s\n" "$f" "$([[ -f "$a" ]] && echo "archive copy missing" || echo "not built yet")"
        fi
    done
    exit 0
fi

mkdir -p "$WD"
echo "repo:    $REPO ($(git -C "$REPO" describe --exact-match --tags 2>/dev/null || git -C "$REPO" rev-parse --short HEAD))"
echo "config:  $CFG   (Q=$Q_DEEP, L=$L_VAL, seed=$CB_SEED)"
echo "workdir: $WD"
echo "archive: $ARCHIVE_WD  (read-only; not touched)"
echo

# One controller for the whole chain: the three rules are a linear dependency, so
# Snakemake sequences them and -j 2 is plenty. `long` because the burn-in rule
# declares a 30-day walltime for the Q=0.01 case even though this Q=1 run needs
# under an hour -- the declaration is what SLURM sees.
# stage1_seed IS the deep history: farm_create_orig_pop_e2.py seeds the coalescent
# genealogy with it, the mutation overlay with seed+1, and the DFE permutation from
# the same generator, and both SLiM stages take it as -s. Without it on BOTH
# --config lists every invocation rebuilds the configfile's 20250303.
JID=$(sbatch --parsable \
    --job-name="deep_history_$CB_SEED" \
    --partition=medium --time=1-00:00:00 --mem=8G --cpus-per-task=1 \
    --output="$WD/rebuild_%j.out" --error="$WD/rebuild_%j.err" \
    --wrap="cd '$WD' && '$SNAKEMAKE' --snakefile '$REPO/Snakefile' \
        --configfile '$REPO/$CFG' --unlock \
        --config stage1_seed=$CB_SEED stage2_seed=$CB_SEED workdir='$WD' publishdir='$WD' || true; \
        '$SNAKEMAKE' --snakefile '$REPO/Snakefile' --configfile '$REPO/$CFG' \
        --profile '$REPO/profiles/o2' --rerun-triggers mtime -j 2 --until stage1 \
        --config stage1_seed=$CB_SEED stage2_seed=$CB_SEED workdir='$WD' publishdir='$WD'")

echo "controller $JID"
echo
echo "Expect ~2 h. Then:"
echo "  bash rebuild_cattle_deep_history.sh --rescale"
echo "  bash rebuild_cattle_deep_history.sh --check     # against the archived one"
