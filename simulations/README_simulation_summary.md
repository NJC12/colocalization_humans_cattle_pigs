# List of simulations



# Baseline cattle population

First, create a large population, add neutral mutations (using coalesence), then add selection coefficients to the formerly neutral mutations.

`/home/njc12/miniconda3/envs/coloc_sims/bin/python /home/njc12/slim_simulations/scripts/farm_create_orig_pop.py --seed 20250303 --Q 0.01 --length $2`

Then burn in the population until it reaches mutation-selection-drift equilibrium.

`/home/njc12/bin/slim/build/slim -m -t -l -s 20250303 -d Q_scaling=0.01 -d length=10000000 -d file_in=\"/home/njc12/slim_simulations/farm_slim_outputs/farm_orig_pop.Q_0.01.L_10000000.seed_20250303.ts\" -d dir_out=\"/n/scratch/users/n/njc12/sims/slim/\" -d file_out=\"farm_burn_in.Q_0.01.L_10000000.seed_20250303\" /home/njc12/slim_simulations/scripts/farm_burn_in_e2.slim`

Now run the rest of the demographic history under selection.

`/home/njc12/bin/slim/build/slim -m -t -l -s 20250303 -d Q_scaling=0.01 -d length=10000000 -d tick=25000 -d file_in=\"/n/scratch/users/n/njc12/sims/slim/farm_burn_in.Q_0.01.L_10000000.seed_20250303.cycle_25000.ts\" -d dir_out=\"/n/scratch/users/n/njc12/sims/slim/\" -d file_out=\"farm_selection.Q_0.01.L_10000000.seed_20250303\" /home/njc12/slim_simulations/scripts/farm_selection.slim`

# Original simulation

These were the simulations included in the preprint. For revision, we also repeated these simulations with nothing changed but the number of simulated GTEx samples.

## Humans

Create the original population.

`/home/njc12/miniconda3/envs/coloc_sims/bin/python /home/njc12/slim_simulations/scripts/human_simulation_o2.py --seed $1 --gwas_h2 0.01 --gtex_h2 0.4 --length 1000000`

Create the GWAS and GTEx data. (N.b., some changes have been made to this code. If there are any issues running it, try with `--already_includes_neutral False`.)

`python /Users/noah/colocalization_humans_cattle_pigs/simulations/create_gwas_files_and_phenotypes.py --gtex_size 1000 --gwas_scaling 35 --gtex_scaling 35 --r2_value 0.2 --min_maf 0.01 --human_ts_file '/Users/noah/comparative_colocalization/data/simulations/demographic_simulations/human/human_selection_Q_1.L10000000.seed_20250521.full.ts' --out_dir '/Users/noah/tmp/selsims/outputs/' --already_includes_neutral True`

## Cattle

Use the data from "Baseline cattle population."

`python /Users/noah/colocalization_humans_cattle_pigs/simulations/create_gwas_files_and_phenotypes.py --gtex_size 1000 --gwas_scaling 35 --gtex_scaling 35 --r2_value 0.2 --min_maf 0.01 --cattle_ts_file '/Users/noah/comparative_colocalization/data/simulations/demographic_simulations/cattle/farm_selection_Q_0.01.L_10000000.seed_20250303.full.ts' --out_dir '/Users/noah/tmp/selsims/outputs/' --already_includes_neutral True`

## Humans (500 GTEx subjects)

`python /Users/noah/colocalization_humans_cattle_pigs/simulations/create_gwas_files_and_phenotypes.py --gtex_size 500 --gwas_scaling 35 --gtex_scaling 35 --r2_value 0.2 --min_maf 0.01 --human_ts_file '/Users/noah/comparative_colocalization/data/simulations/demographic_simulations/human/human_selection_Q_1.L10000000.seed_20250521.full.ts' --out_dir '/Users/noah/tmp/selsims/outputs/' --already_includes_neutral True`

## Cattle (500 GTEx subjects)

`python /Users/noah/colocalization_humans_cattle_pigs/simulations/create_gwas_files_and_phenotypes.py --gtex_size 500 --gwas_scaling 35 --gtex_scaling 35 --r2_value 0.2 --min_maf 0.01 --cattle_ts_file '/Users/noah/comparative_colocalization/data/simulations/demographic_simulations/cattle/farm_selection_Q_0.01.L_10000000.seed_20250303.full.ts' --out_dir '/Users/noah/tmp/selsims/outputs/' --already_includes_neutral True`

# Larger human dataset

Create the original population.

`/home/njc12/miniconda3/envs/coloc_sims/bin/python /home/njc12/slim_simulations/scripts/human_simulation_o2_larger.py --seed 19930224 --gwas_h2 0.01 --gtex_h2 0.4 --length 1000000`

Create the GWAS and GTEx data.
`python /Users/noah/colocalization_humans_cattle_pigs/simulations/create_gwas_files_and_phenotypes.py --gtex_size 50000 --gwas_scaling 35 --gtex_scaling 35 --r2_value 0.2 --min_maf 0.01 --human_ts_file '/Users/noah/tmp/selsims/largesim/hts_19930224_large.ts' --out_dir '/Users/noah/tmp/selsims/largesim/outputs/' --already_includes_neutral True`

# Cattle simulations including positive selection

## With continued bottlenecking

This simulation uses the data from "Baseline cattle population" through epoch 8, at which point it changes the selection.

`/home/njc12/bin/slim/build/slim -m -t -l -s 24 -d selection_multiplier=100 -d selected_generations=23 -d continue_bottlenecking=1 -d num_muts_selected=26 -d Q_scaling=0.01 -d genome_length=10000000 -d file_in=\"/n/data2/hms/dbmi/sunyaev/lab/nconnally/comparative_colocalization/simulations/selection/farm_selection_Q_0.01.L_10000000.seed_20250303.epoch_8.ts\" -d dir_out=\"/n/data2/hms/dbmi/sunyaev/lab/nconnally/tmp/\" -d file_out=\"revision_farm_selection_mult_100_gen_23_muts_26_bottlenecked_sd24\" ../slim_simulations/scripts/revision_pos_sel_farm_selection.slim`

Simulate the traits using the specific mutations whose selection coefficients were changes.

`python /Users/noah/colocalization_humans_cattle_pigs/simulations/create_gwas_files_and_phenotypes.py --gwas_scaling 35 --gtex_scaling 35 --r2_value 0.2 --min_maf 0.01 --cattle_ts_file '/Users/noah/tmp/selsims/revision_farm_selection_mult_100_gen_23_muts_26_bottlenecked_sd24.full.ts' --cattle_m4_file '/Users/noah/tmp/selsims/revision_farm_selection_mult_100_gen_23_muts_26_bottlenecked_sd24.m4_marks.tsv' --out_dir '/Users/noah/tmp/selsims/outdir/'`

## Without continued bottlenecking

This simulation uses the data from "Baseline cattle population" through epoch 8, at which point it changes the selection.

`/home/njc12/bin/slim/build/slim -m -t -l -s 24 -d selection_multiplier=100 -d selected_generations=23 -d continue_bottlenecking=0 -d num_muts_selected=26 -d Q_scaling=0.01 -d genome_length=10000000 -d file_in=\"/n/data2/hms/dbmi/sunyaev/lab/nconnally/comparative_colocalization/simulations/selection/farm_selection_Q_0.01.L_10000000.seed_20250303.epoch_8.ts\" -d dir_out=\"/n/data2/hms/dbmi/sunyaev/lab/nconnally/tmp/\" -d file_out=\"revision_farm_selection_mult_100_gen_23_muts_26_not_bottlenecked_sd24\" ../slim_simulations/scripts/revision_pos_sel_farm_selection.slim`

Simulate the traits using the specific mutations whose selection coefficients were changes.

`python /Users/noah/colocalization_humans_cattle_pigs/simulations/create_gwas_files_and_phenotypes.py --gwas_scaling 35 --gtex_scaling 35 --r2_value 0.2 --min_maf 0.01 --cattle_ts_file '/Users/noah/tmp/selsims/selection_not_bottlenecked/revision_farm_selection_mult_100_gen_23_muts_26_not_bottlenecked_sd24.full.ts' --cattle_m4_file '/Users/noah/tmp/selsims/selection_not_bottlenecked/revision_farm_selection_mult_100_gen_23_muts_26_not_bottlenecked_sd24.m4_marks.tsv' --out_dir '/Users/noah/tmp/selsims/selection_not_bottlenecked/outputs'`