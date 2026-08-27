# Stage 1 for the `human_neutral` pipeline (category H).
#
# A coalescent simulation with no selection anywhere in it. Category A runs the
# same demography forward in SLiM under a three-component DFE and then overlays
# neutral mutations; this runs the demography backwards in msprime and overlays
# ONE class of mutations, none of them selected. The causal variants and their
# effect-size parameters are drawn in stage 2 instead -- see
# create_gwas_files_and_phenotypes.py --synthetic_dfe_effects and
# helpers/synthetic_dfe.py.
#
# WHY THE MUTATION RATE IS THE TOTAL, NOT THE NEUTRAL COMPONENT
#
# Category A mutates at 1.4e-8 per bp per generation, split between a
# DFE-applied 5.6e-9 in SLiM and an 8.4e-9 neutral overlay in msprime. Both
# classes end up in the tree sequence stage 2 reads, so A's variant density is
# set by the total. Simulating only the 8.4e-9 neutral component here would give
# ~40% fewer variants than A, and variant density drives LD, the size of the
# DAP-G candidate set, and therefore fine-mapping difficulty -- so the arms would
# differ in something other than the genetic model. --mut_rate is required rather
# than defaulted for the same reason --recomb_rate is: it must not become an
# invisible literal that can drift away from the arm it is meant to match.
#
# WHY NOT stdpopsim's msprime ENGINE
#
# engine.simulate() would put mutations on with an infinite-sites/JC69 model,
# whose tskit metadata schema is not the SLiM one. get_vars_df reads
# `mut.metadata['mutation_list'][0]['selection_coeff']` and stage 2 calls
# pyslim.generate_nucleotides, both of which need SLiM-shaped mutation metadata.
# So the demography comes from stdpopsim and the ancestry/mutations from msprime
# directly -- the same shape as farm_create_orig_pop_e2.py, which has been doing
# this for the cattle epoch-1 coalescent all along.

import argparse
import os
import resource
import time

import msprime
import numpy as np
import pyslim
import stdpopsim

# The same demographic model the category-A arm passes to the SLiM engine
# (human_simulation_o2.py). Held here as a constant rather than a CLI flag: the
# whole point of this arm is that its demography is A's, so making it settable
# would only create a way for the two to diverge unnoticed.
#
# The same three constants are now also in helpers/human_demography.py, which is
# where human_simulation_o2.py reads them and where the AFR arm (category J)
# gets its population from. They are duplicated rather than imported for the
# reason above: this arm's identity is "A's demography, no selection", and
# H should not follow A if A's population is ever changed.
DEMOGRAPHIC_MODEL = "OutOfAfrica_2T12"
SPECIES = "HomSap"
SAMPLED_POPULATION = "EUR"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, required=True, help='Random seed')
    parser.add_argument('--length', type=float, required=True,
                        help='Length (bp) of the region to simulate')
    parser.add_argument('--n_samples', type=int, required=True,
                        help='Number of EUR individuals to simulate; should equal '
                             'gwas_size + the largest GTEx size the run needs')
    parser.add_argument('--recomb_rate', type=float, required=True,
                        help='Per-bp per-generation recombination rate (config key '
                             '`recombination_rate`). Required, not defaulted: it must '
                             'match the rate the arm this is compared against ran at.')
    parser.add_argument('--mut_rate', type=float, required=True,
                        help='Per-bp per-generation mutation rate (config key '
                             '`mutation_rate`). Use the TOTAL rate of the selected arm '
                             '(1.4e-8), not its neutral component (8.4e-9) -- see the '
                             'module docstring.')
    parser.add_argument('--tmp_dir', type=str, required=False,
                        default='/n/scratch/users/n/njc12/sims/tmp',
                        help='Directory for the output .ts')
    return parser.parse_args()


def remove_fixed(ts):
    """Drop sites that are multiallelic or not segregating in the sample.

    Same rule as human_simulation_o2.py:remove_fixed. `discrete_genome=True`
    places mutations on integer positions, so two mutations can collide onto one
    site; that site is multiallelic and stage 2's biallelic genotype-sum
    assumption would silently mis-count its allele frequency.
    """
    print(f'Initial number of sites: {ts.num_sites}')
    bad_sites = []
    for tree in ts.trees():
        for site in tree.sites():
            if len(site.mutations) != 1:
                bad_sites.append(site.id)
                continue
            mut = site.mutations[0]
            n_derived = tree.num_samples(mut.node)
            if n_derived == 0 or n_derived == ts.num_samples:
                bad_sites.append(site.id)
    ts = ts.delete_sites(bad_sites)
    print(f'Removing {len(bad_sites)} fixed, absent or multiallelic sites')
    print(f'Number of remaining sites: {ts.num_sites}')
    return ts


def simulate(seed, length, n_samples, recomb_rate, mut_rate):
    species = stdpopsim.get_species(SPECIES)
    model = species.get_demographic_model(DEMOGRAPHIC_MODEL)
    print(f'Demography: {DEMOGRAPHIC_MODEL}, populations '
          f'{[p.name for p in model.model.populations]}, sampling '
          f'{n_samples} {SAMPLED_POPULATION} individuals')

    start = time.time()
    ts = msprime.sim_ancestry(
        samples={SAMPLED_POPULATION: n_samples},
        demography=model.model,
        recombination_rate=recomb_rate,
        sequence_length=length,
        random_seed=seed,
    )
    print(f'\tAncestry: {ts.num_trees} trees, {ts.num_samples} sample nodes, '
          f'{time.strftime("%H:%M:%S", time.gmtime(time.time() - start))}')

    # Gives the tree sequence the top-level SLiM metadata that
    # pyslim.generate_nucleotides looks for in stage 2. tick=1 / stage="late"
    # mirror farm_create_orig_pop_e2.py; nothing downstream reads either, since
    # this arm has no forward phase whose tick count would mean anything.
    ts = pyslim.annotate(ts, model_type="WF", tick=1, stage="late")

    start = time.time()
    # type=0 is SLiM's neutral mutation type, so every mutation lands with
    # selection_coeff 0.0 and mutation_type 0 -- exactly what add_neutral()
    # produces for the neutral overlay in the selected arms, which is what makes
    # get_vars_df report selco == 0 across the board here.
    #
    # No Q anywhere in this file: there is no forward phase, so branch lengths are
    # already in generations and the rate is the real per-generation rate. The
    # `* Q` that belongs in the SLiM arms would be a 10x error here.
    ts = msprime.sim_mutations(
        ts, rate=mut_rate, model=msprime.SLiMMutationModel(type=0),
        keep=True, discrete_genome=True, random_seed=seed,
    )
    print(f'\tMutations: {ts.num_mutations} at {ts.num_sites} sites, '
          f'{time.strftime("%H:%M:%S", time.gmtime(time.time() - start))}')

    ts = remove_fixed(ts)

    n = ts.num_samples
    theta_w = ts.segregating_sites(span_normalise=False) / sum(1.0 / i for i in range(1, n))
    print(f'Watterson theta_W = {theta_w / length:.4e} per bp, '
          f'pi = {ts.diversity():.4e}, {ts.num_sites / (length / 1e6):.0f} sites/Mb')
    return ts


def main():
    args = parse_args()
    print(f'Using seed: {args.seed}')
    print(f'L: {args.length:g}, n_samples: {args.n_samples}')
    print(f'recombination rate: {args.recomb_rate:g}, mutation rate: {args.mut_rate:g}')

    tmp_dir = args.tmp_dir.rstrip('/')
    os.makedirs(tmp_dir, exist_ok=True)
    out = f'{tmp_dir}/nts_{args.seed}.ts'

    # Gate on the output stage 2 consumes, matching the human arm's behaviour --
    # a re-run of a partially-completed workdir should not redo the simulation.
    if os.path.exists(out):
        print(f'{out} already exists; skipping the simulation')
        return

    ts = simulate(args.seed, args.length, args.n_samples,
                  args.recomb_rate, args.mut_rate)
    ts.dump(out)
    print(f'Wrote {out}: {ts.num_sites} sites, {ts.num_mutations} mutations')


if __name__ == '__main__':
    started = time.time()
    main()
    usage_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 ** 2)
    print(f'Total time: {time.strftime("%H:%M:%S", time.gmtime(time.time() - started))}')
    print(f'Peak memory usage: {usage_gb:.4f} GB')
