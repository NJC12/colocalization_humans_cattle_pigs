#!/bin/bash
#SBATCH --partition=long
#SBATCH --time=30-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#
# The controller runs on `long` (NOT priority): priority is capped at 2 jobs
# per user and is reserved for the slow cattle stage-1 selection only. The
# controller is lightweight (4 GB/1 CPU, mostly idle polling) and exits as soon
# as its replicate's pipeline finishes, so the 30-day ceiling is just a safety
# margin. The two child rules that otherwise request priority (stage2_split_pheno,
# stage1_human) are forced off priority via --set-resources below.
#
# Lightweight Snakemake controller for one 2 Mb replicate, run AS A SLURM JOB
# (no login-node tmux). It submits the per-rule jobs for this replicate and
# exits when the pipeline finishes. Parameterized via --export environment
# variables set by submit_2Mb.sh:
#   CONFIGFILE  - config/cattle_baseline_from_midpoint_2Mb.yaml or config/human_2Mb.yaml
#   SEED        - this replicate's stage1_seed / stage2_seed
#   WD          - this replicate's workdir
#   PD          - this replicate's publishdir
# The expensive shared cattle selection is NOT here -- it is a separate
# standalone job -- so this controller only manages cheap, short work and never
# needs to outlive a multi-week job.

set -euo pipefail
REPO=/n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake
SNAKEMAKE=/home/njc12/miniconda3/envs/coloc_sims/bin/snakemake

: "${CONFIGFILE:?}" "${SEED:?}" "${WD:?}" "${PD:?}"
# Optional extra --config tokens (e.g. "neutral_trait_vars=True" for category B).
# Intentionally unquoted on the snakemake lines so multiple space-separated
# key=value pairs expand as separate args.
EXTRA_CONFIG="${EXTRA_CONFIG:-}"
mkdir -p "$WD"

# Run FROM the replicate's own workdir so Snakemake creates its .snakemake/
# (lock + iocache + slurm_logs) under $WD, NOT a shared $REPO/.snakemake. When
# many controllers launched at once and all ran from $REPO they raced on that
# shared dir (this killed human_E6). SIM_REPO_DIR is derived from the absolute
# --snakefile path (Snakefile:17), and rule scripts are os.path.join'd onto it,
# so every script/include still resolves against $REPO regardless of CWD.
cd "$WD"

# --unlock is a no-op on a fresh workdir; harmless if a prior controller died.
"$SNAKEMAKE" --snakefile "$REPO/Snakefile" --configfile "$REPO/$CONFIGFILE" --unlock \
    --config stage1_seed="$SEED" stage2_seed="$SEED" workdir="$WD" publishdir="$PD" $EXTRA_CONFIG || true

# No --set-resources overrides. They existed because stage2_split_pheno asked
# for 24h and stage1_human for 7 days, neither of which fits `short`'s 12h cap,
# so both declared `priority` and had to be forced elsewhere. Both are now sized
# from measured round-3 runs (stage2 max 2:17 -> 1h; human stage-1 at L<=4 Mb
# max 1:05 -> 2h) and declare `short` themselves. Nothing lands in priority.
#
# If a 10 Mb config is ever run again, stage1_human's L-aware tier returns it to
# priority/7 days on its own -- reinstate the override then.
"$SNAKEMAKE" --snakefile "$REPO/Snakefile" --configfile "$REPO/$CONFIGFILE" \
    --profile "$REPO/profiles/o2" --rerun-triggers mtime \
    --config stage1_seed="$SEED" stage2_seed="$SEED" workdir="$WD" publishdir="$PD" $EXTRA_CONFIG
