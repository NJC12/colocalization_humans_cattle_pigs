#!/bin/sh
#SBATCH -c 4
#SBATCH -t 30-00:00:00
#SBATCH -p priority
#SBATCH --mem=150G
#SBATCH --output=/home/njc12/slim_simulations/eo/farm_burn_in_e2.%j.out
#SBATCH --error=/home/njc12/slim_simulations/eo/farm_burn_in_e2.%j.err

# Args
# 1: seed
# 2: Q
# 3: length
# 4: file_in (no quotes)
# 5: dir_out (no quotes)
# 6: file_out (no quotes, and leave off .txt)

# /home/njc12/bin/slim/build/slim -m -t -l -s $1 -d Q_scaling=$2 -d length=$3 -d file_in=\"$4\" -d dir_out=\"$5\" -d file_out=\"$6\" /home/njc12/slim_simulations/scripts/farm_burn_in_no_ts_e2.slim
/home/njc12/bin/slim/build/slim -m -t -l -s 20250303 -d Q_scaling=0.01 -d length=10000000 -d file_in=\"/home/njc12/slim_simulations/farm_slim_outputs/farm_orig_pop.Q_0.01.L_10000000.seed_20250303.ts\" -d 
dir_out=\"/n/scratch/users/n/njc12/sims/slim/\" -d file_out=\"farm_burn_in.Q_0.01.L_10000000.seed_20250303\" /home/njc12/slim_simulations/scripts/farm_burn_in_e2.slim
