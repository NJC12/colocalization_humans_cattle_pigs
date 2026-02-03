#!/bin/sh
#SBATCH -c 4
#SBATCH -t 3-00:00:00
#SBATCH -p priority
#SBATCH --mem=50G
#SBATCH --output=/home/njc12/slim_simulations/eo/farm_create_orig_pop.%j.out
#SBATCH --error=/home/njc12/slim_simulations/eo/farm_create_orig_pop.%j.err

# module load gcc/9.2.0
# module load python/3.9.14

source /home/njc12/miniconda3/etc/profile.d/conda.sh
conda activate coloc_sims

# echo $(which python)

echo "Seed: $1"

/home/njc12/miniconda3/envs/coloc_sims/bin/python /home/njc12/slim_simulations/scripts/farm_create_orig_pop.py --seed $1 --Q 0.01 --length $2
