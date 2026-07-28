#!/bin/bash
#SBATCH --job-name=cattle_2Mb_sel
#SBATCH --partition=priority
#SBATCH --time=30-00:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --output=/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb/cattle_baseline/logs/stage1_cattle_selection.slurm_%j.out
#SBATCH --error=/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb/cattle_baseline/logs/stage1_cattle_selection.slurm_%j.err
#
# Shared cattle baseline 2 Mb selection (standalone; not under a controller).
# Reads the EXTRACTED 2 Mb burn-in (helpers/extract_subregion.py output) and
# runs farm_selection.slim epochs 2-12, writing epoch_8.ts (consumed by the 8
# from-midpoint replicates) and full.ts. ~3-4 weeks. Same command shape as the
# 10 Mb run, with length=2000000 and the 2 Mb workdir.

set -euo pipefail
STAGE1=/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb/cattle_baseline/stage1
LOGDIR=/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb/cattle_baseline/logs
mkdir -p "$STAGE1" "$LOGDIR"

bash /home/njc12/bin/slim_5/slim_wrapper.sh -m -t -l \
  -s 20250303 \
  -d Q_scaling=0.01 \
  -d length=2000000 \
  -d tick=25000 \
  -d file_in='"/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb/cattle_baseline/stage1/farm_burn_in.Q_0.01.L_2000000.seed_20250303.cycle_25000.ts"' \
  -d dir_out='"/n/scratch/users/n/njc12/snakemake/simulations_round_2_2Mb/cattle_baseline/stage1/"' \
  -d file_out='"farm_selection_Q_0.01.L_2000000.seed_20250303"' \
  /n/data2/hms/dbmi/sunyaev/lab/nconnally/slim_simulations/snakemake/farm_selection.slim \
  > "$LOGDIR/stage1_cattle_selection.log" 2>&1
