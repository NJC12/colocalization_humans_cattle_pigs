#!/bin/bash
#SBATCH --partition=medium
#SBATCH --time=1-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#
# Controller for ONE replicate of ONE arm of the publication simulation set.
#
# A copy of controller_2Mb_x10.sbatch, which it deliberately does not replace:
# twenty category-O/P controllers were in flight against that file when the
# publication run was prepared, and every SLURM child job reads the repo at job
# start, so editing it under them is not safe. Three differences:
#
#   JOBS   REQUIRED, and passed as -j. The other controller passes no -j at all,
#          so `jobs: 200` from profiles/o2/config.yaml applies. That is correct
#          for a handful of controllers and wrong for this run: 120 controllers x
#          200 is 24,000 queued jobs against an account MaxSubmit of 10,000, on a
#          SHARED account. A refused sbatch is a failed job to Snakemake's SLURM
#          executor; with restart-times: 2 and keep-going: True it retries twice,
#          gives up, and the replicate finishes with no .enloc.sig.out -- which
#          reads downstream as "this arm colocalized nothing" rather than as a
#          gap. The publish gate is what catches it, hours later. Keeping -j low
#          is how it does not happen.
#   UNTIL  optional --until target, so a wave can stop at stage 2 in the FINAL
#          workdir and nothing is thrown away when it continues.
#   time   1 day, not 5. Stage 1 is always adopted here (a symlink job), and the
#          per-replicate DAG is a few hours at these -j values. Sized to still
#          cover the worst case: _dapg_runtime scales with `attempt`, so one
#          walltime-killed hgwas locus can add 60+120+180 min on its own.
#
# Parameterised via --export, set by submit_publication.sh:
#   CONFIGFILE  - this category's config, relative to $REPO
#   SEED        - stage1_seed and stage2_seed for this replicate
#   WD          - this replicate's workdir for THIS arm
#   PD          - where the publish script will copy it (recorded only)
#   STAGE1_SRC  - directory holding this replicate's stage-1 tree sequence
#   EXTRA_CONFIG- the arm: floors, trait counts, sampling scheme, pairing rule
#   JOBS        - max concurrent SLURM jobs for this controller
#   UNTIL       - optional --until target (e.g. stage2)

set -euo pipefail
REPO="${REPO:-/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/publication_repo/simulations}"
SNAKEMAKE="${SNAKEMAKE:-/home/njc12/miniconda3/envs/coloc_sims/bin/snakemake}"

: "${CONFIGFILE:?}" "${SEED:?}" "${WD:?}" "${PD:?}" "${JOBS:?}"

# Intentionally unquoted on the snakemake lines so multiple space-separated
# key=value pairs expand as separate --config args. NO TOKEN MAY CONTAIN A SPACE.
EXTRA_CONFIG="${EXTRA_CONFIG:-}"
UNTIL="${UNTIL:-}"

# Every category in this run adopts its stage-1 tree sequence -- A/E/F/G from the
# archived copy, K and L from A's and E's. So STAGE1_SRC is REQUIRED here, unlike
# in the controller this was copied from. An empty value does not fail loudly: it
# makes Snakemake decide it has to BUILD the tree sequence, and a multi-hour
# genetic simulation starts instead of a symlink.
: "${STAGE1_SRC:?STAGE1_SRC is required -- an empty value silently rebuilds stage 1}"
STAGE1_CFG="stage1_search_dirs=['$STAGE1_SRC']"

mkdir -p "$WD"

# Run FROM the replicate's own workdir so Snakemake's .snakemake/ (lock, iocache,
# slurm_logs) lands under $WD rather than a shared $REPO/.snakemake, which would
# race across concurrent controllers. SIM_REPO_DIR is derived from the absolute
# --snakefile path, so every script and include still resolves against $REPO.
cd "$WD"

# --unlock is a no-op on a fresh workdir; harmless if a prior controller died.
"$SNAKEMAKE" --snakefile "$REPO/Snakefile" --configfile "$REPO/$CONFIGFILE" --unlock \
    --config stage1_seed="$SEED" stage2_seed="$SEED" workdir="$WD" publishdir="$PD" \
    "$STAGE1_CFG" $EXTRA_CONFIG || true

"$SNAKEMAKE" --snakefile "$REPO/Snakefile" --configfile "$REPO/$CONFIGFILE" \
    --profile "$REPO/profiles/o2" --rerun-triggers mtime \
    -j "$JOBS" ${UNTIL:+--until "$UNTIL"} \
    --config stage1_seed="$SEED" stage2_seed="$SEED" workdir="$WD" publishdir="$PD" \
    "$STAGE1_CFG" $EXTRA_CONFIG
