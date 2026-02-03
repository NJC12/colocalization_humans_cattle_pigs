#!/bin/sh
#SBATCH -c 4
#SBATCH -t 30-00:00:00
#SBATCH -p priority
#SBATCH --mem=150G
#SBATCH --output=/home/njc12/slim_simulations/eo/farm_selection.%j.out
#SBATCH --error=/home/njc12/slim_simulations/eo/farm_selection.%j.err

# Args
# 1: seed
# 2: Q
# 3: length
# No 4: file_in (no quotes)
# No 5: file_out (no quotes)
# Tick for starting

/home/njc12/bin/slim/build/slim -m -t -l -s 20250303 -d Q_scaling=0.01 -d length=10000000 -d tick=25000 -d file_in=\"/n/scratch/users/n/njc12/sims/slim/farm_burn_in.Q_0.01.L_10000000.seed_20250303.cycle_25000.ts\" -d dir_out=\"/n/scratch/users/n/njc12/sims/slim/\" -d file_out=\"farm_selection.Q_0.01.L_10000000.seed_20250303\" /home/njc12/slim_simulations/scripts/farm_selection.slim
