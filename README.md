Code for the preprint [Farm animal evolution demonstrates hidden molecular basis of human traits](https://www.biorxiv.org/content/10.64898/2026.02.02.703413v1).

All python and packages for this analysis are in the `comparative_colocalization.yaml` file.

All other software is linked below.

# Simulations

Slim can be downloaded from the [Messer Lab](https://messerlab.org/slim/).

## Simulating genome for human demography

This simulation was run by `run_human_simulation_o2.sh`.
- Most of the purpose of the shell script is to handle modules needed on the compute cluster I used.

The actual simulation occurs in `human_simulation_o2.py`.

Simulation summary
- Simulate selected mutations using the python package stdpopsim.
- Scaling was 10, so the simulated population was 10x smaller and there were 10x fewer generations.
- Add neutral mutations using the coalescence python package msprime.
- Clean up the simulated genomes (e.g., remove triallelic sites).
- Separate into 8,000 subjects for GWAS and 1,000 subjects for GTEx eQTL-mapping, and save each as a vcf.
<!-- - Slightly reformat the files to work more easily with hail.
- Use python package hail to read in files as matrix tables.
- Extract information including the ancestral allele and the selection coefficient.
- For negative selection coefficients, use the python package tstrait to simulate phenotypes.
- Use the genomes and phenotypes to run GWAS and eQTL mapping. Save the results. -->

## Simulating genome for cattle-like demography

Simulation of cattle-like genomes and traits took more steps than human genomes
- We based our population history on McLeod et al. 2013, but used selection coefficients, mutation rates, and recombination rates from the human simulation.
- The final effective population size of Holstein cattle is very small (Ne = 90).
- To simulate GWAS and GTEx, we had to scale the population up 100x, which meant that the number of generations also had to increase 10x.
- The first population epoch of cattle has an Ne of 62,000.
- We could not simulate 6,200,000 cattle.
- So we started in the second epoch, which has an Ne of 17,000.
- To start the simulation, we used msprime coalesence to simulate a neutral population of 1,700,000 cattle (`run_farm_create_orig_pop.sh` and `farm_create_orig_pop_e2.py`).
- We then added selection coefficients to the neutral variants, and ran a forward slim simulation to until these variants reached an equilibrium (`run_farm_burn_in_e2.sh` and `farm_burn_in_e2.slim`).
- We then ran a forward slim simulation to cover all subsequent epoch in the cattle model (`run_farm_selection.sh` and `farm_selection.slim`).

## Simulating traits, GWAS, and eQTL-mapping

Run using `create_gwas_files_and_phenotypes.py`
- Reads in vcfs.
- Extracts relevant information from variants (derived allele frequency, selection coefficient, etc.)
- Uses python package tstrait to simulate phenotypes
- Saves in sbam format for dap-g
- Saves a bunch of shell and plink2 commands

Shell and plink commands
- Slightly reformat data and convert for plink
- Run PCA
- Run GWAS and eQTL-mapping using the genomes and traits from above
 
## Fine-mapping

DAP-G can be downloaded from [its github](https://github.com/xqwen/dap/tree/master?tab=readme-ov-file)
- See installation instructions there.

The script `run_regional_dapg_portion.sh`
- Takes a set of traits to run (a range of lines in sbam files)
- Extracts the relevant region (1 Mb on either side of the selected SNP)
- Runs dap-g
- I was having issues with dap-g grapping as many cores as possible, even when limiting the number of tasks. So the script assigns it specific cores. This might not be necessary on other clusters.
- Commands looked like this `qsub ../run_regional_dapg_portion.sh $PWD cgwas_scaling_35.0 0 100 0.75 10`

## Running fastEnloc

FastEnloc can be downloaded from [its github](https://github.com/xqwen/fastenloc/tree/master?tab=readme-ov-file).
- See installation instructions there.

Steps
- In some cases I created filenames with periods (beyond the period before the file extension)
- This creates problems with fastEnloc, which identifies traits based on the name of the file up to the first period.
- Reformat names in with `for fl in $(ls); do outfl=$(echo $fl | sed -r 's/35.0/35/g'); mv $fl $outfl; done`
- Convert dap-g format to fastEnloc format using the built-in fastEnloc script `dap2enloc`
- `~/bin/fastenloc/utility/dap2enloc -dir cgtex/outputs/ -vcf vcf/cgtex.dap.vcf.gz | gzip - > cgtex_scaling_35.dap.annotation.vcf.gz`
    - If the output is empty, this may be due to a mismatch in SNP IDs. A simple solution is to just made each SNP ID its position.
    - `~/bin/fastenloc/utility/dap2enloc -dir cgtex/outputs/ -vcf <(zcat vcf/cgtex.dap.vcf.gz | awk 'BEGIN{OFS="\t"} {$3=$2; print $0}' | gzip) | gzip - > cgtex_scaling_35.dap.annotation.vcf.gz`
- This produces files ending in `.dap.annotation.vcf.gz`
- We then run fastEnloc using the commands
- Cattle demography: `fastenloc -eqtl cgtex_scaling_35.dap.annotation.vcf.gz -gwas cgwas_scaling_35.dap.annotation.vcf.gz -total_variants [see below]`
- Human demography: `fastenloc -eqtl hgtex_scaling_35.dap.annotation.vcf.gz -gwas hgwas_scaling_35.dap.annotation.vcf.gz -total_variants [see below]`
- The total number of variants is the number of variants shared by the GWAS and eQTL files, and can be found with a few command line tools. (Some systems require `zcat` instead of `gzcat`)
- `comm -12 <(gzcat cgtex_scaling_35.dap.annotation.vcf.gz | awk 'NR > 1 {print $3}' | sort | uniq) <(gzcat cgwas_scaling_35.dap.annotation.vcf.gz | awk 'NR > 1 {print $3}' | sort | uniq) | wc -l`
- `comm -12 <(gzcat hgtex_scaling_35.dap.annotation.vcf.gz | awk 'NR > 1 {print $3}' | sort | uniq) <(gzcat hgwas_scaling_35.dap.annotation.vcf.gz | awk 'NR > 1 {print $3}' | sort | uniq) | wc -l`

# Figures

Each main-text figure has a jupyter notebook (using R).
- Figures 3 and 4 are combined in one notebook, because they mostly use the same data.
- Supplementary figures 1, 2, and 3 are contained in the notebook for figure 2 (so is the allelic age supplementary figure, which at this moment may or may not be included in the final draft).
- Supplementary figure 4 is contained in the notebook for figures 3 and 4.

Each notebook contains:
- What data were used for the figure, and where they can be downloaded
- What processing steps are required post-download (there are *many*)
- The code for the analyses and figures
- The figures themselves
- Notes about the analyses and figures

# License

Copyright 2026 Noah Connally

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

