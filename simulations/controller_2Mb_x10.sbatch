#!/bin/bash
#SBATCH --partition=medium
#SBATCH --time=5-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#
# x10 phenotype-scaling controller for one 2 Mb replicate -- category A (human)
# or category E (cattle baseline from midpoint). Same controller for both; the
# config + STAGE1_SRC passed in decide the species.
#
# The controller runs on `medium` (NOT long, NOT priority): stage 1 (the genetic
# simulation) is REUSED from the x35 run via stage1_search_dirs, so NO real
# genetic simulation runs here. The stage-1 rule still executes, but only as a
# trivial `ln -sf` symlink job -- forced to `short` for human (its default is
# `priority`), already `short` for cattle. The other child jobs are stage 2
# (trait creation, 24 h -> medium) and stages 3-5 (DAP-G / fastEnloc / plink GLM,
# all short). The whole DAG finishes in ~1-2 days, so `medium`'s 5-day ceiling is
# ample and we never occupy `long` or `priority`. The controller itself is
# lightweight (4 GB / 1 CPU, mostly idle polling) and exits when the pipeline
# finishes.
#
# Run AS A SLURM JOB (no login-node tmux). Parameterized via --export environment
# variables set by submit_2Mb_x10.sh (A) / submit_2Mb_x10_E.sh (E):
#   CONFIGFILE  - config/human_2Mb_x10.yaml (A) or
#                 config/cattle_baseline_from_midpoint_2Mb_x10.yaml (E)
#   SEED        - this replicate's stage1_seed / stage2_seed
#   WD          - this replicate's x10 workdir
#   PD          - this replicate's x10 publishdir
#   STAGE1_SRC  - the x35 stage-1 dir to reuse (this replicate's A{N}/ or E{N} stage1 dir)
#   EXTRA_CONFIG- optional extra --config tokens (empty for A and E)

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SNAKEMAKE=/home/njc12/miniconda3/envs/coloc_sims/bin/snakemake

: "${CONFIGFILE:?}" "${SEED:?}" "${WD:?}" "${PD:?}"
# Optional extra --config tokens (e.g. "neutral_trait_vars=True" for category B).
# Intentionally unquoted on the snakemake lines so multiple space-separated
# key=value pairs expand as separate args.
EXTRA_CONFIG="${EXTRA_CONFIG:-}"

# Reuse this replicate's x35 stage-1 tree sequence instead of rerunning it. The
# search-dir lookup is seed-strict, so STAGE1_SRC must be THIS replicate's own
# A{N}/stage1 dir (submit_2Mb_x10.sh guarantees that). Passed as a single quoted
# --config token; Snakemake YAML-parses the value into a one-element list.
#
# EMPTY means "build stage 1 here". It used to be mandatory (`${STAGE1_SRC:?}`),
# which was right while every caller was a multiplier cell reusing an x35 run.
# Category H has no x35 predecessor to reuse -- it is a different stage-1
# pipeline, and a coalescent cheap enough that reuse buys nothing -- so the
# variable is now optional and the search-dir token is omitted when it is unset.
STAGE1_SRC="${STAGE1_SRC:-}"
if [[ -n "$STAGE1_SRC" ]]; then
    STAGE1_CFG="stage1_search_dirs=['$STAGE1_SRC']"
else
    STAGE1_CFG="stage1_search_dirs=[]"
fi

mkdir -p "$WD"

# Run FROM the replicate's own workdir so Snakemake creates its .snakemake/
# (lock + iocache + slurm_logs) under $WD, NOT a shared $REPO/.snakemake, which
# would race across concurrent controllers. SIM_REPO_DIR is derived from the
# absolute --snakefile path, so every script/include still resolves against $REPO.
cd "$WD"

# --unlock is a no-op on a fresh workdir; harmless if a prior controller died.
"$SNAKEMAKE" --snakefile "$REPO/Snakefile" --configfile "$REPO/$CONFIGFILE" --unlock \
    --config stage1_seed="$SEED" stage2_seed="$SEED" workdir="$WD" publishdir="$PD" \
    "$STAGE1_CFG" $EXTRA_CONFIG || true

# No --set-resources overrides. They existed to force stage2_split_pheno (24h
# request) and stage1_human (7-day request) off `priority`, since neither fit
# `short`'s 12h cap. Both are now sized from measured round-3 runs and declare
# `short` themselves, so nothing lands in priority. In these x-cell runs
# stage 1 is reused and does not execute at all.
"$SNAKEMAKE" --snakefile "$REPO/Snakefile" --configfile "$REPO/$CONFIGFILE" \
    --profile "$REPO/profiles/o2" --rerun-triggers mtime \
    --config stage1_seed="$SEED" stage2_seed="$SEED" workdir="$WD" publishdir="$PD" \
    "$STAGE1_CFG" $EXTRA_CONFIG
