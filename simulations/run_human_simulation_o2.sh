#!/bin/sh
#SBATCH -c 4
#SBATCH -t 30-00:00:00
#SBATCH -p priority
#SBATCH --mem=200G
#SBATCH --output=/home/njc12/slim_simulations/eo/human_simulation.%j.out
#SBATCH --error=/home/njc12/slim_simulations/eo/human_simulation.%j.err

module load gcc/9.2.0
module load python/3.9.14
module load java/jdk-11.0.11
module load bcftools/1.14

source /home/njc12/miniconda3/etc/profile.d/conda.sh
conda activate coloc_sims

HAIL_TMPDIR='/n/scratch/users/n/njc12/sims/tmp'
export HAIL_TMPDIR

echo "Seed: $1"

/home/njc12/miniconda3/envs/coloc_sims/bin/python /home/njc12/slim_simulations/scripts/human_simulation_o2.py --seed $1 --gwas_h2 0.01 --gtex_h2 0.4 --length 1000000

