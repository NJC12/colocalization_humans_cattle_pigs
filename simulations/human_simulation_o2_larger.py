# module load gcc/9.2.0
# module load bcftools/1.14
# module load java/jdk-11.0.11

import stdpopsim
import numpy as np
import msprime
import tstrait
import tskit
import subprocess
# import hail as hl
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
# hl.init(global_seed=seed, tmp_dir=temp_dir, local_tmpdir=temp_dir, spark_conf={"spark.local.dir": temp_dir})

Q = 1
# L = 10000000
n_samples = 370_000

##### Running slim ####
if not os.path.exists(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_bad_ids_{seed}_large.vcf') or not os.path.exists(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_bad_ids_{seed}_large.vcf'):
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

    hts =  msprime.sim_mutations(ts, rate = 8.4e-9 * Q, model=mut_model, keep=True, discrete_genome=False, random_seed=seed)
    hts = remove_fixed(hts)


    # samps = [i for i in range(0, 2*n_samples, int(n_samples/5000))] + [i+1 for i in range(0, 2*n_samples, int(n_samples/5000))]
    # samps = [i for i in range(0, 2*n_samples, 18)] + [i+1 for i in range(0, 2*n_samples, 18)]
    samps = [i for i in range(0, 2*n_samples, 10)] + [i+1 for i in range(0, 2*n_samples, 10)]
    samps.sort()
    hgtex = hts.simplify(hts.samples()[samps])
    hgwas = hts.simplify(hts.samples()[[i for i in range(0, 2*n_samples,) if i not in samps]])

    # Save the ts files
    hts.dump(f'/n/scratch/users/n/njc12/sims/tmp/hts_{seed}_large.ts')
    hgtex.dump(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_{seed}_large.ts')
    hgwas.dump(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_{seed}_large.ts')

    print(f'hgtex has {hgtex.num_sites} sites and {hgtex.num_mutations} mutations')
    print(f'hgwas has {hgwas.num_sites} sites and {hgwas.num_mutations} mutations')

    with open(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_bad_ids_{seed}_large.vcf', 'w') as output:
        hgtex.write_vcf(output, allow_position_zero=True)

    with open(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_bad_ids_{seed}_large.vcf', 'w') as output:
        hgwas.write_vcf(output, allow_position_zero=True)


if not os.path.exists(f'/n/scratch/users/n/njc12/sims/tmp/hgtex_{seed}.vcf'):
    subprocess.run(f'/n/app/bcftools/1.21-gcc-14.2.0/bin/bcftools annotate --set-id "%CHROM\\_%POS" /n/scratch/users/n/njc12/sims/tmp/hgtex_bad_ids_{seed}_large.vcf -o /n/scratch/users/n/njc12/sims/tmp/hgtex_{seed}_large.vcf', env=os.environ, shell=True)

if not os.path.exists(f'/n/scratch/users/n/njc12/sims/tmp/hgwas_{seed}.vcf'):
    subprocess.run(f'/n/app/bcftools/1.21-gcc-14.2.0/bin/bcftools annotate --set-id "%CHROM\\_%POS" /n/scratch/users/n/njc12/sims/tmp/hgwas_bad_ids_{seed}_large.vcf -o /n/scratch/users/n/njc12/sims/tmp/hgwas_{seed}_large.vcf', env=os.environ, shell=True)
