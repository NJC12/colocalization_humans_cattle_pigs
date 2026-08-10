# Stage 1 for the `cattle_neutral` pipeline (category I).
#
# A coalescent simulation of the cattle demography with no selection anywhere in
# it. Categories E/F/G run that demography forward in SLiM under a
# three-component DFE; this runs it backwards in msprime and puts on ONE class of
# mutations, none of them selected. The causal variants and their effect-size
# parameters are drawn in stage 2 instead -- see
# create_gwas_files_and_phenotypes.py --synthetic_dfe_effects and
# helpers/synthetic_dfe.py, both of which category H already exercises. The DFE
# is shared between the two species (human_simulation_o2.py:67-92 and
# farm_selection.slim:40-45 define the same three components with the same
# proportions), so H and I draw effect sizes from the SAME distribution and
# differ in exactly one thing: the genealogy.
#
# THE WHOLE HISTORY, NOT JUST EPOCHS 8-12
#
# E/F/G resume from a shared `.ep7.ts` checkpoint because 29,800 ticks of burn-in
# at N=17,000 plus epochs 2-7 is unaffordable to re-run per replicate in a
# forward simulator. A coalescent has no such constraint: epoch 1 is simply the
# ancestral state, reached exactly rather than approached over a finite burn-in.
# So this simulates all twelve epochs, including the full bottleneck from 17,000
# down to 90, and needs no checkpoint and no `cattle_baseline_search_dirs`.
#
# WHY THE MUTATION RATE IS 5.6e-9 AND NOT 1.4e-8
#
# This is the opposite convention from human_neutral_simulation.py, and the
# difference is real. The cattle arms split their 1.4e-8 into a DFE-applied
# 5.6e-9 that SLiM puts on in stage 1 and a neutral 8.4e-9 that
# create_gwas_files_and_phenotypes.py overlays in stage 2 (add_neutral, called
# unconditionally on the --cattle_ts_file branch). The human arms do both halves
# in stage 1, which is why category H asks for the total. Since the cattle
# stage-2 overlay runs whether or not this arm wants it, the only way to land on
# 1.4e-8 without touching shared code is to reproduce the same split. Both
# classes are neutral here, so the split carries no meaning beyond that -- but it
# does make the total provably identical to E/F/G's.
#
# WHY Q APPEARS HERE AND NOT IN human_neutral_simulation.py
#
# Under the Hudson coalescent, running at Q=0.01 in ticks and at Q=1 in
# generations are the same distribution: scaling every size and duration by 1/Q
# and every rate by Q leaves the integral of dt/2N, and mu*t and r*t, unchanged.
# The output would be identical either way. What differs is which process the run
# CLAIMS to approximate. E/F/G are Wright-Fisher populations of 9,000 whose
# entire membership is sampled; Hudson is the standard approximation to that. At
# Q=1 the same 9,000-individual sample would be drawn from a population of 90 --
# 100x oversampled, and an approximation to a process nobody is simulating. So
# the rescaled units are the honest ones, and Q enters exactly where it does in
# farm_create_orig_pop_e2.py, which has been running the cattle epoch-1
# coalescent this way all along.
#
# WHY NOT stdpopsim
#
# There is no BosTau model here to use: the demography is hand-written (see
# helpers/cattle_demography.py). And even for human, stdpopsim's msprime engine
# puts mutations on with a JC69 model whose tskit metadata is not the SLiM
# schema, which get_vars_df and pyslim.generate_nucleotides both require. Hence
# msprime directly, with a SLiMMutationModel.

import argparse
import os
import resource
import sys
import time

import msprime
import pyslim

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from helpers import cattle_demography  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, required=True, help='Random seed')
    parser.add_argument('--length', type=float, required=True,
                        help='Length (bp) of the region to simulate')
    parser.add_argument('--n_samples', type=int, required=True,
                        help='Number of individuals to sample; should equal '
                             'gwas_size + the largest GTEx size the run needs. At '
                             'Q_scaling=0.01 the terminal population is 9,000, which '
                             'is exactly the 8,000 + 1,000 the round-3 cattle arms use '
                             '-- i.e. the sample is the whole population, as in E/F/G.')
    parser.add_argument('--recomb_rate', type=float, required=True,
                        help='Per-bp per-REAL-GENERATION recombination rate (config key '
                             '`recombination_rate`). Q-scaled here, exactly as the '
                             'farm_*.slim scripts do it.')
    parser.add_argument('--mut_rate', type=float, required=True,
                        help='Per-bp per-REAL-GENERATION mutation rate for THIS stage '
                             '(config key `mutation_rate`). Use 5.6e-9, the cattle '
                             'arms\' stage-1 component; stage 2 overlays the remaining '
                             '8.4e-9 -- see the module docstring.')
    parser.add_argument('--Q', type=float, required=True,
                        help='Stage-1 rescaling factor (config key `Q_scaling`, 0.01). '
                             'Sizes and durations are divided by it, rates multiplied '
                             'by it. Must match what stage 2 is told, or get_vars_df '
                             'un-rescales the times by the wrong factor.')
    parser.add_argument('--tmp_dir', type=str, required=False,
                        default='/n/scratch/users/n/njc12/sims/tmp',
                        help='Directory for the output .ts')
    return parser.parse_args()


def remove_fixed(ts):
    """Drop sites that are multiallelic or not segregating in the sample.

    Same rule as human_neutral_simulation.py and human_simulation_o2.py.
    `discrete_genome=True` places mutations on integer positions, so two can
    collide onto one site; that site is multiallelic and stage 2's biallelic
    genotype-sum assumption would silently mis-count its allele frequency.
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


def build_demography(Q):
    """The twelve-epoch schedule as an msprime.Demography, in Q-rescaled ticks."""
    demography = msprime.Demography()
    demography.add_population(name=cattle_demography.POPULATION,
                              initial_size=cattle_demography.terminal_size(Q))
    for time_ago, size, _epoch in cattle_demography.size_changes(Q):
        demography.add_population_parameters_change(
            time=time_ago, initial_size=size,
            population=cattle_demography.POPULATION)
    return demography


def simulate(seed, length, n_samples, recomb_rate, mut_rate, Q):
    print(cattle_demography.describe(Q))

    # Exactly the rescaling farm_create_orig_pop_e2.py:25 and the four
    # farm_*.slim scripts use. Not the linear r*Q approximation: the two agree to
    # ~15 digits at these rates, but the exact form is what SLiM applies and
    # there is no reason for the two engines to differ even in the last bit.
    recomb_per_tick = (1 - (1 - 2 * recomb_rate) ** Q) / 2
    mut_per_tick = mut_rate * Q
    print(f'Per-tick rates at Q={Q:g}: recombination {recomb_per_tick:.6e}, '
          f'mutation {mut_per_tick:.6e} '
          f'(per real generation: {recomb_rate:.6e}, {mut_rate:.6e})')
    print(f'Sampling {n_samples} individuals from a terminal population of '
          f'{cattle_demography.terminal_size(Q):,.0f}')

    start = time.time()
    ts = msprime.sim_ancestry(
        samples={cattle_demography.POPULATION: n_samples},
        demography=build_demography(Q),
        recombination_rate=recomb_per_tick,
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
    # selection_coeff 0.0 and mutation_type 0 -- the same shape add_neutral()
    # produces, which is what makes get_vars_df report selco == 0 across the
    # board here. (farm_create_orig_pop_e2.py uses type=2 because it then
    # overwrites each coefficient from the DFE; there is nothing to overwrite.)
    ts = msprime.sim_mutations(
        ts, rate=mut_per_tick, model=msprime.SLiMMutationModel(type=0),
        keep=True, discrete_genome=True, random_seed=seed,
    )
    print(f'\tMutations: {ts.num_mutations} at {ts.num_sites} sites, '
          f'{time.strftime("%H:%M:%S", time.gmtime(time.time() - start))}')

    ts = remove_fixed(ts)

    # theta_W and pi are 4*Ne*mu per bp -- diversity, not a rate, so neither can
    # be compared against mut_per_tick directly. What IS comparable is the Ne
    # they imply, and it is the sharpest available check on the Q arithmetic:
    # theta_W = 4*(Ne_real/Q)*(mu_real*Q) = 4*Ne_real*mu_real, so dividing by the
    # REAL rate gives a Q-INVARIANT number. Run the same seed at Q=1 and it must
    # come out the same; a `* Q` where a `/ Q` belongs moves it ~100x.
    n = ts.num_samples
    theta_w = ts.segregating_sites(span_normalise=False) / sum(1.0 / i for i in range(1, n))
    theta_w_per_bp = theta_w / length
    print(f'Watterson theta_W = {theta_w_per_bp:.4e} per bp, pi = {ts.diversity():.4e}, '
          f'{ts.num_sites / (length / 1e6):.0f} sites/Mb')
    print(f'Implied harmonic-mean Ne = {theta_w_per_bp / (4 * mut_rate):,.0f} real '
          f'individuals (Q-invariant; this arm\'s stage-1 rate is only 5.6e-9 of '
          f'the eventual 1.4e-8, so it is not comparable to a finished tree sequence)')
    return ts


def main():
    args = parse_args()
    print(f'Using seed: {args.seed}')
    print(f'L: {args.length:g}, n_samples: {args.n_samples}, Q: {args.Q:g}')
    print(f'recombination rate: {args.recomb_rate:g}, mutation rate: {args.mut_rate:g}')

    tmp_dir = args.tmp_dir.rstrip('/')
    os.makedirs(tmp_dir, exist_ok=True)
    out = f'{tmp_dir}/cnts_{args.seed}.ts'

    # Gate on the output stage 2 consumes, matching the other arms' behaviour --
    # a re-run of a partially-completed workdir should not redo the simulation.
    if os.path.exists(out):
        print(f'{out} already exists; skipping the simulation')
        return

    ts = simulate(args.seed, args.length, args.n_samples,
                  args.recomb_rate, args.mut_rate, args.Q)
    ts.dump(out)
    print(f'Wrote {out}: {ts.num_sites} sites, {ts.num_mutations} mutations')


if __name__ == '__main__':
    started = time.time()
    main()
    usage_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 ** 2)
    print(f'Total time: {time.strftime("%H:%M:%S", time.gmtime(time.time() - started))}')
    print(f'Peak memory usage: {usage_gb:.4f} GB')
