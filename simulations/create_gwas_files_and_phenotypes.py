import tskit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pyslim
import msprime
import tstrait
import argparse
import os
import subprocess
import math
import gzip

def remove_fixed(ts):
    # Also removes triallelic sites
    bad_sites = []
    for tree in ts.trees():
        for site in tree.sites():
            if len(site.mutations) > 1:
                bad_sites.append(site.id)
            if len(site.mutations) > 0:
                mut = site.mutations[0]
                daf = tree.num_samples(mut.node) / ts.num_samples
                if daf == 0 or daf == 1:
                    bad_sites.append(site.id)
    return ts.delete_sites(bad_sites)

def relabel_ag_variants(ts, marks_file):
    marks = {}
    with open(marks_file) as f:
        next(f)  # header
        for line in f:
            pos_str, s_str = line.strip().split("\t")
            marks[int(float(pos_str))] = float(s_str)

    tables = ts.dump_tables()
    schema = tables.mutations.metadata_schema

    new_metadata = []
    relabeled = 0
    for row in tables.mutations:
        site_pos = int(ts.site(row.site).position)
        md = row.metadata
        if site_pos in marks:
            target_s = marks[site_pos]
            for ml in md["mutation_list"]:
                if math.isclose(abs(ml["selection_coeff"]), target_s, rel_tol=1e-5):
                    ml["mutation_type"] = 4
                    relabeled += 1
        new_metadata.append(md)

    tables.mutations.packset_metadata(
        [schema.validate_and_encode_row(md) for md in new_metadata]
    )
    ts_relabeled = tables.tree_sequence()
    print('Mutations relabeled as m4: ' + str(relabeled))
    return ts_relabeled

def add_neutral(sts, Q=1, seed=0):
    print(f'The tree sequence began with {sts.num_mutations} mutations')
    mut_rate = 2.36e-8 / Q
    # From https://tskit.dev/pyslim/docs/latest/tutorial.html#sec-tutorial-adding-neutral-mutations
    next_id = pyslim.next_slim_mutation_id(sts)
    ts = msprime.sim_mutations(
        sts,
        rate=mut_rate,
        model=msprime.SLiMMutationModel(type=0, next_id=next_id),
        keep=True,
        random_seed=seed
    )

    print(f"The tree sequence now has {ts.num_mutations} mutations,\n"
        f"and mean pairwise nucleotide diversity is {ts.diversity():0.3e}.")
    
    return ts

def get_vars_df(ts, Q=1):
    tree = pd.DataFrame({
        'id': [v.site.id for v in ts.variants()],
        'type': [v.site.mutations[0].metadata['mutation_list'][0]['mutation_type'] for v in ts.variants()],
        'selco': [v.site.mutations[0].metadata['mutation_list'][0]['selection_coeff'] for v in ts.variants()],
        'daf': [sum(v.genotypes)/18000 for v in ts.variants()],
        'time': [v.site.mutations[0].time for v in ts.variants()],
        'position': [v.site.position for v in ts.variants()]
    })

    if Q != 1:
        tree['selco'] = tree['selco']*Q
        tree['time'] = tree['time']/Q

    tree['maf'] = tree['daf'].apply(lambda x: min([x, 1-x]))
    tree['Vs'] = np.abs(tree['selco']) * tree['daf'] * (1-tree['daf'])
    
    return tree

def generate_phenos_from_pos(position, ts, vars, scaling=1):
    site_id = vars[vars['position'] == position]['id']
    site_id = int(site_id.iloc[0])
    selco = vars[vars['position'] == position]['selco']
    selco = float(selco.iloc[0])
    # beta = selco_to_beta(selco)
    beta = np.sqrt(np.abs(selco)) * scaling
    position = int(position)

    pheno = tstrait.sim_phenotype(ts, model=tstrait.TraitModelFixed(beta, random_sign=True), causal_sites=[site_id], h2=1, random_seed=position)
    pheno.phenotype['environmental_noise'] = np.random.normal(0, 1, pheno.phenotype.shape[0])
    pheno.phenotype['trait_id'] = 'tr' + str(int(position))
    pheno.trait['trait_id'] = 'tr' + str(int(position))
    pheno.phenotype['phenotype'] = pheno.phenotype['genetic_value'] + pheno.phenotype['environmental_noise']

    return pheno

def combine_phenos_to_df(positions, ts, vars, scaling):
    id_list = ['tsk_' + str(x.id) for x in ts.individuals()]
    phenos = pd.DataFrame({'FID': 0, 'IID': id_list})
    snp_key = pd.DataFrame({'position': [], 'site_id': [], 'effect_size': [], 'causal_allele': [], 'allele_freq': [], 'trait_id': []})

    for pos in positions['position']:
        tr = generate_phenos_from_pos(pos, ts, vars, scaling)
        tr_id = tr.trait.trait_id.iloc[0]
        tr_vals = tr.phenotype.phenotype
        # phenos[tr_id] = tr_vals
        # phenos = pd.concat([phenos, pd.DataFrame({'fid': 0, 'iid': id_list, tr_id: tr_vals})], axis=1)
        phenos = pd.concat([phenos, pd.DataFrame({tr_id: tr_vals})], axis=1)
        snp_key = pd.concat([snp_key, tr.trait])
    
    return snp_key, phenos

def write_ts_as_sbams(ts, output_path):
    # Stream genotypes directly to an SBAMS file without materializing the full matrix.
    with open(output_path, 'w') as f:
        for var in ts.variants():
            g = var.genotypes                     # 1-D array, haplotypes for this site only
            diploid = g[0::2] + g[1::2]           # combine pairs into per-individual dosages
            snp_id = 'snp' + str(int(var.site.position))
            # Format: geno <tab> snp_id <tab> group_id=1 <tab> dosages...
            f.write('geno\t')
            f.write(snp_id)
            f.write('\t1\t')
            f.write('\t'.join(map(str, diploid.tolist())))
            f.write('\n')


def write_traits_as_sbams(df, output_path):
    # Stream the transposed phenotype table without building the transpose in memory.
    trait_cols = [c for c in df.columns if c not in ('FID', 'IID')]
    with open(output_path, 'w') as f:
        for trait in trait_cols:
            values = df[trait].to_numpy()         # one column at a time
            f.write('pheno\t')
            f.write(str(trait))
            f.write('\t1\t')
            f.write('\t'.join(map(str, values.tolist())))
            f.write('\n')

# def traits_to_sbams(df):
#     sbams = df.rename(columns={'IID': 'pheno_id'}).drop('FID', axis=1).T.reset_index()
#     colnames = sbams.iloc[0].values
#     sbams = sbams.iloc[1:, ]
#     sbams.columns = colnames
#     sbams.insert(0, 'pheno', 'pheno')
#     sbams.insert(2, 'group_id', 1)
#     return sbams

# def ts_to_sbams(ts):
#     gm = ts.genotype_matrix()
#     # This has each chromosome, which we want to combine into values for individuals
#     chrom1 = gm[:, 0::2]
#     chrom2 = gm[:, 1::2]
#     inds = pd.DataFrame(chrom1 + chrom2)
#     inds.columns = ['tsk_' + str(i) for i in range(inds.shape[1])]
#     inds.insert(0, 'group_id', 1)
#     inds.insert(0, 'snp_id', ['snp' + str(int(v.site.position)) for v in ts.variants()])
#     inds.insert(0, 'geno', 'geno')
#     return(inds)

def create_pca(ts):
    pca = ts.pca(20, individuals=[i.id for i in ts.individuals()])
    pca = pd.DataFrame(pca.factors.T)
    pca.insert(0, 'group', 1)
    pca.insert(0, 'variable_name', ['pca'+str(i) for i in range(20)])
    pca.insert(0, 'controlled', 'controlled')
    return pca

def get_commands(gwas_scaling, gtex_scaling, min_maf, sim_dir, run_hums, run_cows):
    commands = []
    commands.append('# To run locally')
    commands.append(f'cd {sim_dir}')
    commands.append('mkdir -p plink_analysis/glm')
    commands.append('cp *.vcf plink_analysis')
    commands.append('mv *.pheno plink_analysis')
    if run_cows:
        commands.append('awk \'BEGIN{OFS="\\t"} $1 == 1 {$3="snp"$2; print $0}\' cgwas.vcf | gzip > cgwas.dap.vcf.gz')
        commands.append('awk \'BEGIN{OFS="\\t"} $1 == 1 {$3="snp"$2; print $0}\' cgtex.vcf | gzip > cgtex.dap.vcf.gz')
    if run_hums:
        commands.append('awk \'BEGIN{OFS="\\t"} $1 == 1 {$3="snp"$2; print $0}\' hgwas.vcf | gzip > hgwas.dap.vcf.gz')
        commands.append('awk \'BEGIN{OFS="\\t"} $1 == 1 {$3="snp"$2; print $0}\' hgtex.vcf | gzip > hgtex.dap.vcf.gz')
    commands.append('scp -A -P 2222 -J njc12@transfer02.o2.rc.hms.harvard.edu * nconnally@127.0.0.1:/net/home/nconnally/comparative_colocalization/fastenloc_simulations/tmp_upload')
    commands.append('scp -A -P 2222 -J njc12@transfer01.o2.rc.hms.harvard.edu * nconnally@127.0.0.1:/net/home/nconnally/comparative_colocalization/fastenloc_simulations/tmp_upload')
    commands.append('')
    commands.append('# To run on muhee')
    commands.append(f'cd /net/home/nconnally/comparative_colocalization/fastenloc_simulations')
    commands.append(f'mkdir gwas_scaling_{gwas_scaling}_gtex_scaling_{gtex_scaling}_min_maf_{min_maf}')
    commands.append(f'cd gwas_scaling_{gwas_scaling}_gtex_scaling_{gtex_scaling}_min_maf_{min_maf}')
    commands.append(f'mkdir cgtex cgwas hgtex hgwas eo pca vcf')
    commands.append('mv ../tmp_upload/*vcf.gz vcf')
    commands.append('mv ../tmp_upload/*pca pca')
    if run_cows:
        commands.append('mv ../tmp_upload/cgtex*.sbams cgtex')
        commands.append('mv ../tmp_upload/cgwas*.sbams cgwas')
    if run_hums:
        commands.append('mv ../tmp_upload/hgtex*.sbams hgtex')
        commands.append('mv ../tmp_upload/hgwas*.sbams hgwas')
    base_submission = f'qsub ../run_dapg_portion.sh $PWD'
    if run_cows:
        commands.append(f'{base_submission} cgtex_scaling_{gtex_scaling} 0 200 0.75 16')
        commands.append(f'{base_submission} cgtex_scaling_{gtex_scaling} 200 400 0.75 24')
        commands.append(f'{base_submission} cgwas_scaling_{gwas_scaling} 0 200 0.75 32')
        commands.append(f'{base_submission} cgwas_scaling_{gwas_scaling} 200 400 0.75 40')
    if run_hums:
        commands.append(f'{base_submission} hgtex_scaling_{gtex_scaling} 0 100 0.75 48')
        commands.append(f'{base_submission} hgwas_scaling_{gwas_scaling} 0 100 0.75 56')
    commands.append('')
    commands.append('# plink commands')
    commands.append('# "plink2" must call the plink2 binary')
    commands.append(f'cd {sim_dir}/plink_analysis')
    if run_cows:
        commands.append('# cgwas')
        commands.append('plink2 --vcf cgwas.vcf --make-pgen --out cgwas')
        commands.append('plink2 --pfile cgwas --maf 0.01 --indep-pairwise 500 0.5 --out cgwas_ldpruned')
        commands.append('plink2 --pfile cgwas --extract cgwas_ldpruned.prune.in --pca approx --out cgwas_pca')
        commands.append(f'plink2 --pfile cgwas --glm hide-covar --covar cgwas_pca.eigenvec --pheno cgwas_traits.scaling_{gwas_scaling}.pheno --out glm/cgwas_glm_scaling_{gwas_scaling}')
        commands.append('# cgtex')
        commands.append('plink2 --vcf cgtex.vcf --make-pgen --out cgtex')
        commands.append('plink2 --pfile cgtex --maf 0.01 --indep-pairwise 500 0.5 --out cgtex_ldpruned')
        commands.append('plink2 --pfile cgtex --extract cgtex_ldpruned.prune.in --pca approx --out cgtex_pca')
        commands.append(f'plink2 --pfile cgtex --glm hide-covar --covar cgtex_pca.eigenvec --pheno cgtex_traits.scaling_{gtex_scaling}.pheno --out glm/cgtex_glm_scaling_{gtex_scaling}')
    if run_hums:
        commands.append('# hgwas')
        commands.append('plink2 --vcf hgwas.vcf --make-pgen --out hgwas')
        commands.append('plink2 --pfile hgwas --maf 0.01 --indep-pairwise 500 0.5 --out hgwas_ldpruned')
        commands.append('plink2 --pfile hgwas --extract hgwas_ldpruned.prune.in --pca approx --out hgwas_pca')
        commands.append(f'plink2 --pfile hgwas --glm hide-covar --covar hgwas_pca.eigenvec --pheno hgwas_traits.scaling_{gwas_scaling}.pheno --out glm/hgwas_glm_scaling_{gwas_scaling}')
        commands.append('# hgtex')
        commands.append('plink2 --vcf hgtex.vcf --make-pgen --out hgtex')
        commands.append('plink2 --pfile hgtex --maf 0.01 --indep-pairwise 500 0.5 --out hgtex_ldpruned')
        commands.append('plink2 --pfile hgtex --extract hgtex_ldpruned.prune.in --pca approx --out hgtex_pca')
        commands.append(f'plink2 --pfile hgtex --glm hide-covar --covar hgtex_pca.eigenvec --pheno hgtex_traits.scaling_{gtex_scaling}.pheno --out glm/hgtex_glm_scaling_{gtex_scaling}')
    command_string = '\n'.join(commands)
    return command_string

def get_arguments():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--gwas_scaling', type=int)
    parser.add_argument('--gtex_scaling', type=int)
    parser.add_argument('--r2_value', type=float) # Note to self: I think I no longer use this. Confirm, then delete.
    parser.add_argument('--min_maf', type=float, help='The minor allele frequency variants must reach to be tested')
    parser.add_argument('--human_ts_file', type=str)
    parser.add_argument('--cattle_ts_file', type=str)
    parser.add_argument('--cattle_m4_file', type=str, help='A file keeping track of the cattle variants under intense agricultural selection')
    parser.add_argument('--out_dir', type=str)
    parser.add_argument('--n_samples', type=float, required=False, default=None, help='The number of individuals in the GTEx sample and GWAS sample combined (default is the full number of simulated individuals)')
    parser.add_argument('--gtex_size', type=float, required=False, default=500, help='The number of individuals in the GTEx sample (default 500)')
    parser.add_argument('--already_includes_neutral', type=bool, required=False, default=False, help='Only needed if you already added neutral mutations in the input data')
    return parser.parse_args()

if __name__ == '__main__':
    print('Getting arguments')
    args = get_arguments()
    gwas_scaling = args.gwas_scaling
    gtex_scaling = args.gtex_scaling
    r2_value = args.r2_value
    min_maf = args.min_maf

    run_hums = not args.human_ts_file is None
    run_cows = not args.cattle_ts_file is None

    # Clean up the data
    print('Loading data')
    if run_cows:
        # cows = tskit.load('/Users/noah/comparative_colocalization/data/simulations/demographic_simulations/cattle/farm_selection_Q_0.01.L_10000000.seed_20250303.full.ts')
        cows = tskit.load(args.cattle_ts_file)
        if args.cattle_m4_file is not None:
            cows = relabel_ag_variants(cows, args.cattle_m4_file)
        print('Cattle demography')
        cows = remove_fixed(add_neutral(cows, Q=100, seed=20250303))
        cows = pyslim.convert_alleles(pyslim.generate_nucleotides(cows))
    if run_hums:
        # hums = tskit.load('/Users/noah/comparative_colocalization/data/simulations/demographic_simulations/human/human_selection_Q_1.L10000000.seed_20250521.full.ts')
        hums = tskit.load(args.human_ts_file)
        print('Human demography')
        if not args.already_includes_neutral:
            hums = remove_fixed(add_neutral(hums, seed=20250521))
            hums = pyslim.convert_alleles(pyslim.generate_nucleotides(hums))

    # Split into GWAS and GTEx
    print('Splitting into GWAS and GTEx')
    # n_samples = 9000
    n_samples = args.n_samples
    if run_hums:
        if n_samples is None:
            n_samples = hums.num_individuals
        print('Human samples: ' + str(n_samples))
        samps = [i for i in range(0, 2*n_samples, int(n_samples/(args.gtex_size/2)))] + [i+1 for i in range(0, 2*n_samples, int(n_samples/(args.gtex_size/2)))]
        samps.sort()
        hgtex = hums.simplify(hums.samples()[samps])
        hgwas = hums.simplify(hums.samples()[[i for i in range(0, 2*n_samples,) if i not in samps]])
    if run_cows:
        if n_samples is None:
            n_samples = cows.num_individuals
        print('Cattle samples: ' + str(n_samples))
        samps = [i for i in range(0, 2*n_samples, int(n_samples/(args.gtex_size/2)))] + [i+1 for i in range(0, 2*n_samples, int(n_samples/(args.gtex_size/2)))]
        samps.sort()
        cgtex = cows.simplify(cows.samples()[samps])
        cgwas = cows.simplify(cows.samples()[[i for i in range(0, 2*n_samples,) if i not in samps]])

    # Get the variants
    print('Getting variants')
    if run_cows:
        cvars = get_vars_df(cows, Q=100)
        cgwas_vars = get_vars_df(cgwas, Q=100)
        cgtex_vars = get_vars_df(cgtex, Q=100)
        print('Number of cattle vars: ' + str(len(cvars)))
    if run_hums:
        hvars = get_vars_df(hums)
        hgwas_vars = get_vars_df(hgwas)
        hgtex_vars = get_vars_df(hgtex)
        print('Number of human vars: ' + str(len(hvars)))

    # Select the relevant alleles
    print('Selecting relevant alleles')
    if run_hums:
        hgtex_maf_dict = {var.position: var.maf for _, var in hgtex_vars.iterrows()}
        hcausative_maf01 = hgwas_vars.loc[(hgwas_vars['maf'] >= min_maf) & (hgwas_vars['selco'] != 0) & (hgwas_vars['position'] > 5e5) & (hgwas_vars['position'] < 9.5e6) & (hgwas_vars.position.isin(hgtex_vars.position))]
        print(f'Number of causative human variants: {hcausative_maf01.shape[0]}')
        hcausative_maf01 = hcausative_maf01[hcausative_maf01.apply(lambda row: hgtex_maf_dict[row.position] >= min_maf, axis=1)]
    if run_cows:
        cgtex_maf_dict = {var.position: var.maf for _, var in cgtex_vars.iterrows()}
        ccausative_maf01 = cgwas_vars.loc[(cgwas_vars['maf'] >= min_maf) & (cgwas_vars['selco'] != 0) & (cgwas_vars['position'] > 5e5) & (cgwas_vars['position'] < 9.5e6) & (cgwas_vars.position.isin(cgtex_vars.position))]
        print(f'Number of causative cattle variants: {ccausative_maf01.shape[0]}')
        ccausative_maf01 = ccausative_maf01[ccausative_maf01.apply(lambda row: cgtex_maf_dict[row.position] >= min_maf, axis=1)]

    # Creating phenotypes
    print('Creating phenotypes')
    if run_cows:
        cgwas_key_maf01, cgwas_traits_maf01 = combine_phenos_to_df(ccausative_maf01, cgwas, cgwas_vars, scaling=gwas_scaling)
        cgtex_key_maf01, cgtex_traits_maf01 = combine_phenos_to_df(ccausative_maf01, cgtex, cgtex_vars, scaling=gtex_scaling)
    if run_hums:
        hgwas_key_maf01, hgwas_traits_maf01 = combine_phenos_to_df(hcausative_maf01, hgwas, hgwas_vars, scaling=gwas_scaling)
        hgtex_key_maf01, hgtex_traits_maf01 = combine_phenos_to_df(hcausative_maf01, hgtex, hgtex_vars, scaling=gtex_scaling)

    # Write everything out
    print('Writing everything out')
    # sim_dir = f'/Users/noah/comparative_colocalization/data/simulations/gwas_and_eqtl_mapping/gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}'
    sim_dir = f'{args.out_dir}/gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}'
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    if run_cows:
        cgwas_vars.to_csv(f'{sim_dir}/cgwas_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)
        cgtex_vars.to_csv(f'{sim_dir}/cgtex_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)
        cgwas_traits_maf01.to_csv(f'{sim_dir}/cgwas_traits_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)
        cgtex_traits_maf01.to_csv(f'{sim_dir}/cgtex_traits_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)

        with open (f'{sim_dir}/cgwas.vcf', 'w') as f:
            cgwas.write_vcf(f)
        with open (f'{sim_dir}/cgtex.vcf', 'w') as f:
            cgtex.write_vcf(f)
        
        # ts_to_sbams(cgtex).to_csv(f'{sim_dir}/cgtex_scaling_{gtex_scaling}_geno.sbams', sep='\t', header=False, index=False)
        # traits_to_sbams(cgtex_traits_maf01).to_csv(f'{sim_dir}/cgtex_scaling_{gtex_scaling}_pheno.sbams', sep='\t', header=False, index=False)
        # ts_to_sbams(cgwas).to_csv(f'{sim_dir}/cgwas_scaling_{gwas_scaling}_geno.sbams', sep='\t', header=False, index=False)
        # traits_to_sbams(cgwas_traits_maf01).to_csv(f'{sim_dir}/cgwas_scaling_{gwas_scaling}_pheno.sbams', sep='\t', header=False, index=False)
        write_ts_as_sbams(cgwas, f'{sim_dir}/cgwas_scaling_{gwas_scaling}_geno.sbams')
        write_traits_as_sbams(cgwas_traits_maf01, f'{sim_dir}/cgwas_scaling_{gwas_scaling}_pheno.sbams')
        write_ts_as_sbams(cgtex, f'{sim_dir}/cgtex_scaling_{gtex_scaling}_geno.sbams')
        write_traits_as_sbams(cgtex_traits_maf01, f'{sim_dir}/cgtex_scaling_{gtex_scaling}_pheno.sbams')

        create_pca(cgtex).to_csv(f'{sim_dir}/cgtex.pca', sep='\t', header=False, index=False)
        create_pca(cgwas).to_csv(f'{sim_dir}/cgwas.pca', sep='\t', header=False, index=False)

        # And we're also going to create the output for plink -> sumstats -> coloc
        cgwas_traits_maf01.to_csv(f'{sim_dir}/cgwas_traits.scaling_{gwas_scaling}.pheno', sep='\t', index=False)
        cgtex_traits_maf01.to_csv(f'{sim_dir}/cgtex_traits.scaling_{gtex_scaling}.pheno', sep='\t', index=False)

    if run_hums:
        hgwas_vars.to_csv(f'{sim_dir}/hgwas_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)
        hgtex_vars.to_csv(f'{sim_dir}/hgtex_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)
        hgwas_traits_maf01.to_csv(f'{sim_dir}/hgwas_traits_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)
        hgtex_traits_maf01.to_csv(f'{sim_dir}/hgtex_traits_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)

        with gzip.open(f'{sim_dir}/hgwas.vcf.gz', 'wt') as f:
            hgwas.write_vcf(f)
        with gzip.open(f'{sim_dir}/hgtex.vcf.gz', 'wt') as f:
            hgtex.write_vcf(f)

        # ts_to_sbams(hgtex).to_csv(f'{sim_dir}/hgtex_scaling_{gtex_scaling}_geno.sbams', sep='\t', header=False, index=False)
        # traits_to_sbams(hgtex_traits_maf01).to_csv(f'{sim_dir}/hgtex_scaling_{gtex_scaling}_pheno.sbams', sep='\t', header=False, index=False)
        # ts_to_sbams(hgwas).to_csv(f'{sim_dir}/hgwas_scaling_{gwas_scaling}_geno.sbams', sep='\t', header=False, index=False)
        # traits_to_sbams(hgwas_traits_maf01).to_csv(f'{sim_dir}/hgwas_scaling_{gwas_scaling}_pheno.sbams', sep='\t', header=False, index=False)
        write_ts_as_sbams(hgwas, f'{sim_dir}/hgwas_scaling_{gwas_scaling}_geno.sbams')
        write_traits_as_sbams(hgwas_traits_maf01, f'{sim_dir}/hgwas_scaling_{gwas_scaling}_pheno.sbams')
        write_ts_as_sbams(hgtex, f'{sim_dir}/hgtex_scaling_{gtex_scaling}_geno.sbams')
        write_traits_as_sbams(hgtex_traits_maf01, f'{sim_dir}/hgtex_scaling_{gtex_scaling}_pheno.sbams')

        create_pca(hgtex).to_csv(f'{sim_dir}/hgtex.pca', sep='\t', header=False, index=False)
        create_pca(hgwas).to_csv(f'{sim_dir}/hgwas.pca', sep='\t', header=False, index=False)

        hgwas_traits_maf01.to_csv(f'{sim_dir}/hgwas_traits.scaling_{gwas_scaling}.pheno', sep='\t', index=False)
        hgtex_traits_maf01.to_csv(f'{sim_dir}/hgtex_traits.scaling_{gtex_scaling}.pheno', sep='\t', index=False)

    # Write out a long string to the file sim_dir/gwas_scaling_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.txt
    with open(f'{sim_dir}/gwas_scaling_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.txt', 'w') as f:
        f.write(get_commands(gwas_scaling, gtex_scaling, min_maf, sim_dir, run_hums, run_cows))



# Included in final analysis

# Not included in final analysis

# Too soon to tell
# python ../create_gwas_files_and_phenotypes.py --gwas_scaling 35 --gtex_scaling 35 --r2_value 0.2 --min_maf 0.01 --cattle_ts_file '/Users/noah/tmp/selsims/revision_farm_selection_mult_100_gen_23_muts_26_bottlenecked_sd24.full.ts' --cattle_m4_file '/Users/noah/tmp/selsims/revision_farm_selection_mult_100_gen_23_muts_26_bottlenecked_sd24.m4_marks.tsv' --out_dir '/Users/noah/tmp/selsims/outdir/'
# python ../create_gwas_files_and_phenotypes.py --gtex_size 50000 --gwas_scaling 35 --gtex_scaling 35 --r2_value 0.2 --min_maf 0.01 --human_ts_file '/Users/noah/tmp/selsims/largesim/hts_19930224_large.ts' --out_dir '/Users/noah/tmp/selsims/largesim/outputs/' --already_includes_neutral True
