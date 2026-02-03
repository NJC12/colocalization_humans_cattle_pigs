# module load gcc/9.2.0
# module load bcftools/1.14
# module load java/jdk-11.0.11

import stdpopsim
import numpy as np
import msprime
import tstrait
import tskit
import subprocess
import hail as hl
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, required=True, help='Random seed')
parser.add_argument('--gwas_h2', type=float, required=True, help='Per-SNP heritability of the GWAS trait')
parser.add_argument('--gtex_h2', type=float, required=True, help='Per-SNP heritability of the GTEx trait')
parser.add_argument('--length', type=float, required=True, help='Length (in basepairs) of the region to simulate')
args = parser.parse_args()
seed = args.seed
gwas_h2 = args.gwas_h2
gtex_h2 = args.gtex_h2
L = args.length
print(f'Using seed: {seed}')
print(f'GWAS h2: {gwas_h2}')
print(f'GTEx h2: {gtex_h2}')

temp_dir='/n/scratch/users/n/njc12/sims/tmp'
hl.init(global_seed=seed, tmp_dir=temp_dir, local_tmpdir=temp_dir, spark_conf={"spark.local.dir": temp_dir})

Q = 10
# L = 10000000
n_samples = 9000

##### Running slim ####
if not os.path.exists(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_bad_ids_{seed}.vcf') or not os.path.exists(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_bad_ids_{seed}.vcf'):
    print('Running slim')
    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model("OutOfAfrica_2T12")
    engine = stdpopsim.get_engine("slim")

    # contig = species.get_contig("chr22", left=1e7, right=1.1e7)

    contig = species.get_contig(length=L)
    samples = {"AFR": 0, "EUR": n_samples}
    dfe = species.get_dfe("Gamma_K17")
    contig.add_dfe(intervals=np.array([[0, L]]), DFE=dfe)

    ts = engine.simulate(model, contig, samples, slim_scaling_factor=Q)


    def remove_fixed(ts):
        # Also removes triallelic
        print(f'Initial number of sites: {ts.num_sites}')
        bad_sites = []
        for tree in ts.trees():
            for site in tree.sites():
                if len(site.mutations) > 1 or len(site.mutations) == 0:
                    bad_sites.append(site.id)
                elif len(site.mutations) == 1:
                    mut = site.mutations[0]
                    # daf = tree.num_samples(mut.node) / ts.num_samples
                    # if daf == 0 or daf == 1:
                    #     bad_sites.append(site.id)
                    if tree.num_samples(mut.node) == 0 or tree.num_samples(mut.node) == ts.num_samples:
                        bad_sites.append(site.id)
                        # print('Bad: ' + str(tree.num_samples(mut.node)))
                    else:
                        # print('Good: ' + str(tree.num_samples(mut.node)))
                        # Add the allele frequency to the mutation metadata
                        site.metadata = {'allele_frequency': tree.num_samples(mut.node) / ts.num_samples}
                # else:
                #     site
        # print(bad_sites[:10])
        ts = ts.delete_sites(bad_sites)
        print(f"Removing {len(bad_sites)} fixed sites")
        print(f"Number of remaining sites: {ts.num_sites}")
        return ts


    mut_model = msprime.SLiMMutationModel(type=2)

    hts =  msprime.sim_mutations(ts, rate = 8.4e-9 * Q, model=mut_model, keep=True, discrete_genome=False, random_seed=seed)
    hts = remove_fixed(hts)


    samps = [i for i in range(0, 2*n_samples, int(n_samples/500))] + [i+1 for i in range(0, 2*n_samples, int(n_samples/500))]
    samps.sort()
    hgtex = hts.simplify(hts.samples()[samps])
    hgwas = hts.simplify(hts.samples()[[i for i in range(0, 2*n_samples,) if i not in samps]])

    # Save the ts files
    hts.dump(f'/n/scratch/users/n/njc12/sims/tmp/hts_{seed}.ts')
    hgtex.dump(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_{seed}.ts')
    hgwas.dump(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_{seed}.ts')

    print(f'hgtex has {hgtex.num_sites} sites and {hgtex.num_mutations} mutations')
    print(f'hgwas has {hgwas.num_sites} sites and {hgwas.num_mutations} mutations')

    with open(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_bad_ids_{seed}.vcf', 'w') as output:
        hgtex.write_vcf(output, allow_position_zero=True)

    with open(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_bad_ids_{seed}.vcf', 'w') as output:
        hgwas.write_vcf(output, allow_position_zero=True)


if not os.path.exists(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_{seed}.vcf'):
    subprocess.run(f'/n/app/bcftools/1.14/bin/bcftools annotate --set-id "%CHROM\\_%POS" /n/scratch/users/n/njc12/sims/tmp/hgtex_bad_ids_{seed}.vcf -o /n/scratch/users/n/njc12/sims/tmp/hgtex_{seed}.vcf', env=os.environ, shell=True)

if not os.path.exists(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_{seed}.vcf'):
    subprocess.run(f'/n/app/bcftools/1.14/bin/bcftools annotate --set-id "%CHROM\\_%POS" /n/scratch/users/n/njc12/sims/tmp/hgwas_bad_ids_{seed}.vcf -o /n/scratch/users/n/njc12/sims/tmp/hgwas_{seed}.vcf', env=os.environ, shell=True)

# if not os.path.exists(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_{seed}.mt'):
#     hl.import_vcf(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_bad_ids_{seed}.vcf').write(f'/n/scratch/users/n/njc12/sims/tmp/tmp_hgtex_{seed}.mt', overwrite=True)
#     hgtex_mt = hl.read_matrix_table(f'/n/scratch/users/n/njc12/sims/tmp/tmp_hgtex_{seed}.mt')
#     pruned_hgtex_rows = hl.ld_prune(hgtex_mt.GT, r2=0.2, bp_window_size=500000)
#     pruned_hgtex_mt = hgtex_mt.filter_rows(hl.is_defined(pruned_hgtex_rows[hgtex_mt.row_key]))
#     eigenvalues, pcs, _ = hl.hwe_normalized_pca(pruned_hgtex_mt.GT) # The unused third output is the loadings
#     hgtex_mt = hgtex_mt.annotate_cols(scores=pcs[hgtex_mt.s].scores)
#     hgtex_mt.write(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_{seed}.mt', overwrite=True)

# if not os.path.exists(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_{seed}.mt'):
#     hl.import_vcf(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_bad_ids_{seed}.vcf').write(f'/n/scratch/users/n/njc12/sims/tmp/tmp_hgwas_{seed}.mt', overwrite=True)
#     hgwas_mt = hl.read_matrix_table(f'/n/scratch/users/n/njc12/sims/tmp/tmp_hgwas_{seed}.mt')
#     pruned_hgwas_rows = hl.ld_prune(hgwas_mt.GT, r2=0.2, bp_window_size=500000)
#     pruned_hgwas_mt = hgwas_mt.filter_rows(hl.is_defined(pruned_hgwas_rows[hgwas_mt.row_key]))
#     eigenvalues, pcs, _ = hl.hwe_normalized_pca(pruned_hgwas_mt.GT) # The unused third output is the loadings
#     hgwas_mt = hgwas_mt.annotate_cols(scores=pcs[hgwas_mt.s].scores)
#     hgwas_mt.write(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_{seed}.mt', overwrite=True)


# # Read in the files
# # hts = tstrait.load_tree_sequence(f'/n/scratch/users/n/njc12/sims/tmp/hts_{seed}.ts')
# # hgtex = tstrait.load_tree_sequence(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_{seed}.ts')
# # hgwas = tstrait.load_tree_sequence(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_{seed}.ts')
# hgtex = tskit.load(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_{seed}.ts')
# hgwas = tskit.load(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_{seed}.ts')

# hgtex_mt = hl.read_matrix_table(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_{seed}.mt')
# hgwas_mt = hl.read_matrix_table(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_{seed}.mt')

# # def get_selco_from_site(site):
# #     return site.mutations[0].metadata['mutation_list'][0]['selection_coeff']

# def get_selco_from_var(var):
#     return var.site.mutations[0].metadata['mutation_list'][0]['selection_coeff']

# def get_daf1_freq_from_var(var):
#     ancestral_state = var.states('ancestral_state')[0]
#     return var.frequencies()[ancestral_state]

# # sites = [x for x in hts.sites()]

# def get_vars_in_range(v, a, b, maf):
#     return [x for x in v if \
#             (a <= get_selco_from_var(x) <= b or a >= get_selco_from_var(x) >= b) and \
#             get_daf1_freq_from_var(x) >= maf and get_daf1_freq_from_var(x) <= 1-maf and \
#             x.site.position > 500000 and x.site.position < L - 500000]

# # This is important because it sets the range of selection coefficients we examine
# candidate_gtex_vars = get_vars_in_range(hgtex.variants(), -0.00001, -5, 0.001)
# candidate_gwas_vars = get_vars_in_range(hgwas.variants(), -0.00001, -5, 0.001)

# causative_vars = [x for x in candidate_gtex_vars if x in candidate_gwas_vars]
# print(f'Number of causative vars: {len(causative_vars)}')

# for var in causative_vars:
#     position = var.site.position
#     lower_bound = position - 500000
#     upper_bound = position + 500000
#     print(f'Testing position {position}')
#     id = var.site.id
#     if not os.path.exists(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_trait_{seed}_pos_{position}.tsv'):
#         # Simulate the phenotype
#         hgtex_traits = tstrait.sim_phenotype(hgtex, tstrait.TraitModelNormal(0,1), causal_sites=[id], h2=gtex_h2)
#         hgtex_traits.phenotype['individual_id'] = ['tsk_' + str(i) for i in hgtex_traits.phenotype['individual_id']]
#         hgtex_traits.phenotype.to_csv(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_traits_{seed}_pos_{position}.tsv', sep='\t', index=False)
#         # Add traits to (temporary copy of) the matrix table
#         hgtex_tab = hl.import_table(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_traits_{seed}_pos_{position}.tsv', impute=True).key_by('individual_id')
#         tmp_hgtex_mt = hgtex_mt.annotate_cols(pheno=hgtex_tab[hgtex_mt.s])
#         # Extract a more manageable region
#         tmp_hgtex_mt = tmp_hgtex_mt.filter_rows((tmp_hgtex_mt.locus.position >= lower_bound) & (tmp_hgtex_mt.locus.position <= upper_bound))
#         # Run and save the GWAS
#         tmp_hgtex_mapping = hl.linear_regression_rows(y=tmp_hgtex_mt.pheno.phenotype, x=tmp_hgtex_mt.GT.n_alt_alleles(), covariates=[1.0, tmp_hgtex_mt.scores[0], tmp_hgtex_mt.scores[1], tmp_hgtex_mt.scores[2], tmp_hgtex_mt.scores[3], tmp_hgtex_mt.scores[4], tmp_hgtex_mt.scores[5], tmp_hgtex_mt.scores[6], tmp_hgtex_mt.scores[7], tmp_hgtex_mt.scores[8], tmp_hgtex_mt.scores[9]])
#         tmp_hgtex_mapping.row.export(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_mapping_{seed}_pos_{position}.tsv')

#     if not os.path.exists(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_trait_{seed}_pos_{position}.tsv'):
#         # Simulate the phenotype
#         hgwas_traits = tstrait.sim_phenotype(hgwas, tstrait.TraitModelNormal(0,1), causal_sites=[id], h2=gwas_h2)
#         hgwas_traits.phenotype['individual_id'] = ['tsk_' + str(i) for i in hgwas_traits.phenotype['individual_id']]
#         hgwas_traits.phenotype.to_csv(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_traits_{seed}_pos_{position}.tsv', sep='\t', index=False)
#         # Add traits to (temporary copy of) the matrix table
#         hgwas_tab = hl.import_table(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_traits_{seed}_pos_{position}.tsv', impute=True).key_by('individual_id')
#         tmp_hgwas_mt = hgwas_mt.annotate_cols(pheno=hgwas_tab[hgwas_mt.s])
#         # Extract a more manageable region
#         tmp_hgwas_mt = tmp_hgwas_mt.filter_rows((tmp_hgwas_mt.locus.position >= lower_bound) & (tmp_hgwas_mt.locus.position <= upper_bound))
#         # Run and save the GWAS
#         tmp_hgwas_mapping = hl.linear_regression_rows(y=tmp_hgwas_mt.pheno.phenotype, x=tmp_hgwas_mt.GT.n_alt_alleles(), covariates=[1.0, tmp_hgwas_mt.scores[0], tmp_hgwas_mt.scores[1], tmp_hgwas_mt.scores[2], tmp_hgwas_mt.scores[3], tmp_hgwas_mt.scores[4], tmp_hgwas_mt.scores[5], tmp_hgwas_mt.scores[6], tmp_hgwas_mt.scores[7], tmp_hgwas_mt.scores[8], tmp_hgwas_mt.scores[9]])
#         tmp_hgwas_mapping.row.export(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_mapping_{seed}_pos_{position}.tsv')

