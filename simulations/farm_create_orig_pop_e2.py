import numpy as np
import pyslim
import msprime
import argparse
import time
import resource

# Code adapted from pyslim and stdpopsim websites
# def epoch1_demographic_model(Q, seq_len, farm_recomb_rate):
def epoch1_demographic_model(seq_len, farm_recomb_rate, pop_size=17000):
    demog_model = msprime.Demography()
    demog_model.add_population(initial_size=pop_size // Q)
    ots = msprime.sim_ancestry(
            samples=pop_size / Q,
            demography=demog_model,
            recombination_rate= (1 - (1 - 2 * farm_recomb_rate)**Q)/2,
            sequence_length=seq_len)
    # "Adds default metadata to everything that needs it"
    # nts: the tick is 1 because the corresponds to slim only, right? msprime has no ticks
    ots = pyslim.annotate(ots, model_type="WF", tick=1, stage="late")
    return ots

# def add_muts(ots, Q, farm_mutation_rate):
def add_muts(ots, farm_mutation_rate):
    mut_model = msprime.SLiMMutationModel(type=2)
    ots = msprime.sim_mutations(
            ots,
            rate=farm_mutation_rate * Q,
            model=mut_model,
            keep=True)
    print(f'\tThe tree sequence now has {ots.num_mutations} mutations, at {ots.num_sites} distinct sites.', flush=True)
    if ots.num_mutations != ots.num_sites:
        ots = remove_bad_sites(ots)
    return ots

def m1_muts(n):
    muts = rng.exponential(scale=0.5e-4, size=n)
    return -1*Q*muts

def m2_muts(n):
    muts = rng.gamma(shape=0.206, scale=0.03/0.206, size=n)
    return -1*Q*muts

def m3_muts(n):
    muts = rng.exponential(scale=0.005, size=n)
    return Q*muts

def selection_coeff_bulk(num_muts):
    mut_types = rng.multinomial(num_muts, [0.83, 0.16, 0.01])

    m1s = m1_muts(mut_types[0])
    m2s = m2_muts(mut_types[1])
    m3s = m3_muts(mut_types[2])

    type_list = np.array([int(1)]*mut_types[0] + [int(2)]*mut_types[1] + [int(3)]*mut_types[2])
    mutation_selection_coefficients = np.concatenate((m1s, m2s, m3s))
    p = np.random.permutation(len(type_list))
    
    return mutation_selection_coefficients[p], type_list[p]

def add_selection_coeffs_bulk(ots):
    new_tables = ots.dump_tables()
    new_tables.mutations.clear()
    
    selection_coefficients, type_list = selection_coeff_bulk(ots.num_mutations)

    for mut, sel, mtype in zip(ots.mutations(), selection_coefficients, type_list):
        md_list = mut.metadata['mutation_list'].copy()
        md_list[0]['selection_coeff'] = sel
        md_list[0]['mutation_type'] = int(mtype)

        _ = new_tables.mutations.append(mut.replace(metadata={'mutation_list': md_list}))

    return new_tables.tree_sequence()

def bad_site(var):
    # Is it triallelic?
    if var.num_alleles > 2:
        return True
    # Is there a mutation that has become fixed?
    if any(f == 0 for f in list(var.frequencies().values())):
        return True
    
    return False

def remove_bad_sites(tree):
    vars = tree.variants()
    bad_sites = [v.site.id for v in vars if bad_site(v)]
    ts = tree.delete_sites(bad_sites)
    print(f'\tRemoving {len(bad_sites)} fixed or triallelic sites', flush=True)
    print(f'\tNumber of remaining sites: {ts.num_sites}', flush=True)
    return ts

def main():
    # Read in the command line arguments seed and Q
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, required=True, help='Random seed')
    parser.add_argument('--Q', type=float, required=True, help='Slim scaling factor')
    parser.add_argument('--length', type=int, required=True, help='Sequence length')
    args = parser.parse_args()
    seed = args.seed
    length = args.length

    global Q
    Q = args.Q

    print(f'Using seed {seed} and Q {Q}', flush=True)

    global rng
    rng = np.random.default_rng(seed)
    # farm_mutation_rate = 8.4e-9 / 40
    farm_mutation_rate = 2.36e-8 / 40
    farm_recomb_rate = 1.3e-8
    # Default drift barrier value is 0

    print('Creating demographic model', flush=True)
    start_time = time.time()
    # Measure the memory usage but still assign the output to ots
    ots = epoch1_demographic_model(length, farm_recomb_rate)
    # Print the time taken in hours, minutes, and seconds
    print(f'\tTime taken: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}', flush=True)
       
    print('Adding mutations', flush=True)
    start_time = time.time()
    ots = add_muts(ots, farm_mutation_rate)
    print(f'\tTime taken: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}', flush=True)
        
    print('Adding selection coefficients')
    start_time = time.time()
    ots = add_selection_coeffs_bulk(ots)
    print(f'\tTime taken: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}', flush=True)

    out_name = f'/home/njc12/slim_simulations/farm_slim_outputs/farm_orig_pop.Q_{Q}.L_{length}.seed_{seed}.ts'
    ots.dump(out_name)


if __name__ == "__main__":
    full_start_time = time.time()
    main()
    usage_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    usage_gb = usage_kb / (1024 ** 2)
    print('Finished creating epoch 1')
    print(f'\tTotal time: {time.strftime("%H:%M:%S", time.gmtime(time.time() - full_start_time))}')
    print(f'\tPeak memory usage: {usage_gb:.4f} GB')

