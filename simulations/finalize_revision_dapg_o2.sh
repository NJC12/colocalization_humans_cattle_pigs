#!/bin/bash
#SBATCH -J dapg_finalize
#SBATCH -p short
#SBATCH -t 0-02:00
#SBATCH -c 1
#SBATCH --mem=4G
#
# Publish dap-g outputs from /n/scratch back to /n/data2 once the array job
# is done. Runs with afterany so partial results from a failed array still
# get published. Idempotent (rsync -a --update).
#
# Required env (passed via sbatch --export):
#   STAGE_DIR        /n/scratch/.../<sim>/<cat_abbv>
#   PUBLISH_DIR      /n/data2/.../<sim>/<cat_abbv>

set -euo pipefail

: "${STAGE_DIR:?STAGE_DIR not set}"
: "${PUBLISH_DIR:?PUBLISH_DIR not set}"

mkdir -p "$PUBLISH_DIR/outputs" "$PUBLISH_DIR/logs"

echo "[$(date -Is)] rsync $STAGE_DIR/outputs/ -> $PUBLISH_DIR/outputs/"
rsync -a --update "$STAGE_DIR/outputs/" "$PUBLISH_DIR/outputs/"

echo "[$(date -Is)] rsync $STAGE_DIR/logs/ -> $PUBLISH_DIR/logs/"
rsync -a --update "$STAGE_DIR/logs/" "$PUBLISH_DIR/logs/"

n_out=$(find "$PUBLISH_DIR/outputs" -type f -name '*.dapg.out' | wc -l)
echo "[$(date -Is)] published $n_out output files to $PUBLISH_DIR/outputs/"
