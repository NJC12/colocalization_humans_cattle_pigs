#!/bin/bash
# Fetch the LARGE 2 Mb result files that exceed the remote-bridge 1 MB download
# cap: *_pheno.sbams and *_vars_*_maf_*.tsv, for categories A, B, E, F (10 reps
# each). Pulls them FLAT into the matching local replicate directory, e.g.
#   simulation_data/r2_25_2Mb/B/B4/hgwas_scaling_35_pheno.sbams
#
# The small files (*.enloc.sig.out, *_glm_lead_snps.tsv) were already fetched
# via the bridge; this script only handles the big ones.
#
# Transfers go through the O2 transfer node (transfer.rc.hms.harvard.edu), which
# requires ONE interactive Duo push -- an SSH ControlMaster is opened up front so
# you authenticate once, then all 40 replicates reuse that connection.
#
# Run locally:  bash simulations/fetch_big_results_2Mb.sh

set -uo pipefail

USER=njc12
HOST=transfer.rc.hms.harvard.edu
SCRATCH=/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb
SUB=stage2/*/gwas_35_gtex_35_maf_0.01
LOCAL=/Users/noah/colocalization_humans_cattle_pigs/simulation_data/r2_25_2Mb

SOCK=$(mktemp -u "${TMPDIR:-/tmp}/o2_xfer_XXXXXX.sock")
cleanup() { ssh -S "$SOCK" -O exit "$USER@$HOST" 2>/dev/null; }
trap cleanup EXIT

echo "Opening SSH master connection (approve the Duo push)..."
ssh -fNM -S "$SOCK" -o ControlPersist=30m -o ServerAliveInterval=30 "$USER@$HOST" || {
  echo "ERROR: could not establish SSH connection to $HOST" >&2; exit 1; }

RSYNC_E="ssh -S $SOCK"
ok=0; warn=0
for CAT in A B E F; do
  for N in $(seq 1 10); do
    ID="${CAT}${N}"
    DEST="$LOCAL/$CAT/$ID"
    mkdir -p "$DEST"
    SRC="$SCRATCH/$ID/$SUB/*_pheno.sbams $SCRATCH/$ID/$SUB/*_vars_*_maf_*.tsv"
    if rsync -a --info=name0 -e "$RSYNC_E" "$USER@$HOST:$SRC" "$DEST/"; then
      cnt=$(ls "$DEST"/*_pheno.sbams "$DEST"/*_vars_*_maf_*.tsv 2>/dev/null | wc -l | tr -d ' ')
      echo "$ID: $cnt big files"
      ok=$((ok+1))
    else
      echo "$ID: WARN rsync returned non-zero" >&2
      warn=$((warn+1))
    fi
  done
done
echo
echo "Done. reps_ok=$ok warn=$warn"
echo "Expected per rep: 4 *_pheno.sbams + 4 *_vars_*_maf_*.tsv = 8 big files."
