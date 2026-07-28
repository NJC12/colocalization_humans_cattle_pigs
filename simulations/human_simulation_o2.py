# module load gcc/9.2.0
# module load bcftools/1.14
# module load java/jdk-11.0.11

import stdpopsim
import numpy as np
import msprime
import tstrait
import tskit
# import hail as hl
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, required=True, help='Random seed')
parser.add_argument('--gwas_h2', type=float, required=True, help='Per-SNP heritability of the GWAS trait')
parser.add_argument('--gtex_h2', type=float, required=True, help='Per-SNP heritability of the GTEx trait')
parser.add_argument('--length', type=float, required=True, help='Length (in basepairs) of the region to simulate')
parser.add_argument('--Q', type=int, required=True, help='SLiM scaling factor (e.g. 10 for ~9k indiv runs, 1 for ~150k indiv runs)')
parser.add_argument('--n_samples', type=int, required=True, help='Number of EUR individuals to simulate; should equal GWAS_size + largest GTEx_size for the downstream simulation')
parser.add_argument('--tmp_dir', type=str, required=False, default='/n/scratch/users/n/njc12/sims/tmp',
                    help='Directory for SLiM/intermediate .ts/.vcf files')
parser.add_argument('--recomb_rate', type=float, required=True,
                    help='Per-bp per-generation recombination rate (config key '
                         '`recombination_rate`). Required rather than defaulted: leaving it '
                         'to stdpopsim gave a generic-contig flat map of 1.282e-8 while the '
                         'cattle arm ran at 9.26e-9, so the two arms differed in recombination '
                         'as well as demography. Must match the --recomb_rate passed to '
                         'farm_create_orig_pop_e2.py.')
args = parser.parse_args()
seed = args.seed
gwas_h2 = args.gwas_h2
gtex_h2 = args.gtex_h2
L = args.length
Q = args.Q
n_samples = args.n_samples
recomb_rate = args.recomb_rate
print(f'Using seed: {seed}')
print(f'GWAS h2: {gwas_h2}')
print(f'GTEx h2: {gtex_h2}')
print(f'Q: {Q}')
print(f'n_samples: {n_samples}')
print(f'recombination rate: {recomb_rate}')

temp_dir = args.tmp_dir.rstrip('/')
os.makedirs(temp_dir, exist_ok=True)
# hl.init(global_seed=seed, tmp_dir=temp_dir, local_tmpdir=temp_dir, spark_conf={"spark.local.dir": temp_dir})

##### Running slim ####
# Gate on the only output anything downstream consumes. rules/stage1_human.smk
# copies hts_{seed}.ts and nothing else; stage 2 re-does the GWAS/GTEx split
# from it. This used to gate on the per-sample VCFs written at the end of this
# block, which are no longer produced.
if not os.path.exists(f'{temp_dir}/hts_{seed}.ts'):
    print('Running slim')
    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model("OutOfAfrica_2T12")
    engine = stdpopsim.get_engine("slim")

    # contig = species.get_contig("chr22", left=1e7, right=1.1e7)

    # recombination_rate is pinned rather than inherited: without it stdpopsim
    # hands back a flat map at the HomSap generic-contig rate (1.282e-8), which
    # is not what the cattle arm uses. Demography is meant to be the only
    # difference between the arms.
    contig = species.get_contig(length=L, recombination_rate=recomb_rate)
    contig.mutation_rate = 5.6e-9  # unscaled; stdpopsim's slim_scaling_factor multiplies by Q
    samples = {"AFR": 0, "EUR": n_samples}

    mt1 = stdpopsim.MutationType(
        dominance_coeff=0.5,
        distribution_type="e",
        distribution_args=[-4e-4],
    )
    mt2 = stdpopsim.MutationType(
        dominance_coeff=0.5,
        distribution_type="g",
        distribution_args=[-0.03, 0.206],
    )
    mt3 = stdpopsim.MutationType(
        dominance_coeff=0.5,
        distribution_type="e",
        distribution_args=[0.005],
    )
    dfe = stdpopsim.DFE(
        id="custom_3comp",
        description="m1 exp del, m2 gamma del, m3 exp pos",
        long_description="Shared DFE across cattle and human simulations",
        mutation_types=[mt1, mt2, mt3],
        proportions=[0.975, 0.024, 0.001],
    )
    contig.add_dfe(intervals=np.array([[0, L]]), DFE=dfe)

    ts = engine.simulate(model, contig, samples, slim_scaling_factor=Q)
    print(ts)

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

    # discrete_genome=True: mutations placed on integer positions. Collisions
    # become multi-allelic sites, which get filtered by remove_fixed() below.
    # Avoids the ~6% duplicate-integer-position problem that discrete_genome=False
    # caused (float positions rounded to ints by tskit.write_vcf).
    # Rate is deliberately NOT multiplied by Q. engine.simulate() has already run
    # stdpopsim's _recap_and_rescale, which does `table.time *= slim_scaling_factor`
    # on the node, migration and mutation tables, so branch lengths here are in
    # real generations and want the real per-generation rate. The old `* Q` gave
    # 10x too many neutral sites (208,464 variants in a 2 Mb sample where ~32,000
    # is correct). Selected mutations come from SLiM at the correct rate and were
    # unaffected, so the error surfaced as a distorted selected:neutral ratio
    # rather than an obviously broken run.
    hts = msprime.sim_mutations(ts, rate=8.4e-9, model=mut_model, keep=True, discrete_genome=True, random_seed=seed)
    hts = remove_fixed(hts)

    hts.dump(f'{temp_dir}/hts_{seed}.ts')
    print(f'hts has {hts.num_sites} sites and {hts.num_mutations} mutations')
else:
    print(f'{temp_dir}/hts_{seed}.ts already exists; skipping SLiM')
