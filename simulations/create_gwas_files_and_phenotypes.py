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


def str2bool(s):
    # argparse's `type=bool` is broken: bool("False") is True. Use this instead.
    if isinstance(s, bool):
        return s
    sl = str(s).strip().lower()
    if sl in ('true', '1', 'yes', 'y', 't'):
        return True
    if sl in ('false', '0', 'no', 'n', 'f'):
        return False
    raise argparse.ArgumentTypeError(f'Boolean value expected, got {s!r}')

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

def remove_position_zero(ts):
    """Drop any site at integer position 0.

    tskit's `write_vcf` refuses to emit a variant at VCF position 0 (the
    VCF spec is 1-indexed); with `discrete_genome=True`, msprime can place a
    mutation at integer position 0 (seed-dependent). Causative trait variants
    are restricted to [0.5 Mb, 9.5 Mb] downstream, so a position-0 variant
    can never become a trait variant -- filtering it out upstream is harmless
    and keeps every downstream `write_vcf` valid. See notebook entry
    `b3-stage2-position-zero` (B3 / seed 23).
    """
    bad_sites = [s.id for s in ts.sites() if int(s.position) == 0]
    if bad_sites:
        print(f'Dropping {len(bad_sites)} site(s) at position 0')
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

def add_neutral(sts, Q=1, seed=0, handoff_ticks=None, deep_Q_scaling=None):
    """Overlay neutral mutations on a raw SLiM tree sequence.

    `Q` is the RECIPROCAL of the config's Q_scaling (callers pass 1/q_scaling),
    so mut_rate = 8.4e-9/Q == 8.4e-9 * Q_scaling. Raw SLiM branch lengths are in
    ticks, and at Q_scaling=0.01 one tick is 0.01 real generations, so the rate
    per tick has to shrink by the same factor for the product to come out at
    8.4e-9 per bp per real generation.

    `handoff_ticks` switches this to a two-pass overlay, for tree sequences that
    span the Q=1 -> Q=0.01 deep-history handoff. Such a sequence has a piecewise
    time scale -- one tick is one real generation before the handoff and 0.01
    after it -- so a single rate is wrong for one of the two segments. Since the
    deep genealogy carries most of the total branch length, applying the recent
    (0.01-scaled) rate throughout would lose ~99% of the neutral variation: the
    same class of error as the human arm's `8.4e-9 * Q` overlay.

    Pass `handoff_ticks` = the number of ticks between the handoff and the end of
    the simulation (epochs 8-12 = 600+600+600+300+300 = 2400 at Q_scaling=0.01),
    and `deep_Q_scaling` = the Q_scaling the deep phase ran under (1). Mutations
    younger than the handoff get 8.4e-9 * Q_scaling per tick; older ones get
    8.4e-9 * deep_Q_scaling.

    Verify with Watterson's theta on the final tree sequence: it should land at
    8.4e-9 per bp per real generation either way.
    """
    print(f'The tree sequence began with {sts.num_mutations} mutations')
    recent_rate = 8.4e-9 / Q
    # From https://tskit.dev/pyslim/docs/latest/tutorial.html#sec-tutorial-adding-neutral-mutations
    next_id = pyslim.next_slim_mutation_id(sts)

    if handoff_ticks is None:
        ts = msprime.sim_mutations(
            sts,
            rate=recent_rate,
            model=msprime.SLiMMutationModel(type=0, next_id=next_id),
            keep=True,
            random_seed=seed
        )
    else:
        if deep_Q_scaling is None:
            raise ValueError('handoff_ticks requires deep_Q_scaling')
        deep_rate = 8.4e-9 * deep_Q_scaling
        print(f'Piecewise neutral overlay at the Q handoff {handoff_ticks} ticks before '
              f'the end: {recent_rate:.3e}/tick more recently, {deep_rate:.3e}/tick deeper')
        # end_time / start_time bound the branch material each pass can mutate,
        # measured as time-ago, so the two passes partition the tree exactly.
        ts = msprime.sim_mutations(
            sts,
            rate=recent_rate,
            end_time=handoff_ticks,
            model=msprime.SLiMMutationModel(type=0, next_id=next_id),
            keep=True,
            random_seed=seed
        )
        ts = msprime.sim_mutations(
            ts,
            rate=deep_rate,
            start_time=handoff_ticks,
            model=msprime.SLiMMutationModel(
                type=0, next_id=pyslim.next_slim_mutation_id(ts)),
            keep=True,
            random_seed=seed + 1
        )

    print(f"The tree sequence now has {ts.num_mutations} mutations,\n"
        f"and mean pairwise nucleotide diversity is {ts.diversity():0.3e}.")

    return ts

def get_vars_df(ts, Q_scaling=1.0, times_already_unscaled=False,
                handoff_ticks=None, deep_Q_scaling=None):
    """Tabulate the variants, undoing SLiM's Q rescaling.

    `Q_scaling` is the stage-1 rescaling factor (the config key of the same
    name; cattle 0.01, human 10). BOTH engines write s_slim = Q_scaling *
    s_real into the SLiM script -- raw SLiM via `initializeMutationType(...,
    Q * -4e-4)` in farm_*.slim, stdpopsim via slim_engine.py's
    `"Q * " + <mean>` (dfe.py:Q_scaled_index marks the mean of "e"/"g"
    distributions) -- so the recorded selection_coeff is divided by Q_scaling
    here to recover the real coefficient. This feeds
    generate_phenos_from_pos(), which takes beta = sqrt(|selco|) * scaling, so
    an uncorrected selco rescales every effect size by sqrt(Q_scaling) and the
    per-variant variance explained by Q_scaling.

    `times_already_unscaled` is True for stdpopsim output: its
    _recap_and_rescale does `table.time *= slim_scaling_factor` on the node,
    migration and mutation tables, so times come back in real generations --
    but it never touches selection_coeff. Raw SLiM output (cattle) gets
    neither correction, so its times still need multiplying by Q_scaling.
    Conflating the two is what left the human arm's selco at 10x: the cattle
    call correctly passed the reciprocal of Q_scaling, the human call passed
    nothing.
    """
    # One pass over ts.variants(), not six. The previous version built each
    # column from its own comprehension, so every variant was decoded six times;
    # it also summed genotypes with Python's builtin sum() over a numpy array,
    # which at 208k sites x 16k haplotypes is billions of interpreted-level
    # additions. This is the dominant cost of stage 2.
    #
    # The denominator is ts.num_samples, NOT 2 * ts.num_samples. num_samples
    # counts sample NODES (haplotypes), so the old expression reported exactly
    # half the true allele frequency -- verified against plink's A1_FREQ on
    # round-2 output, ratio 2.0000 in both arms. Consequences of the old bug:
    # every daf/maf/Vs in the vars tables was halved, and the `min_maf` filter
    # was effectively 2x its stated value (0.02 rather than 0.01), which bites
    # the rare-skewed human SFS harder than the bottlenecked cattle one.
    n_hap = ts.num_samples
    ids, types, selcos, dafs, times, positions = [], [], [], [], [], []
    for v in ts.variants():
        mut = v.site.mutations[0]
        md = mut.metadata['mutation_list'][0]
        ids.append(v.site.id)
        types.append(md['mutation_type'])
        selcos.append(md['selection_coeff'])
        # Sites are biallelic here (remove_fixed drops multiallelics upstream),
        # so genotypes are 0/1 and the sum is the derived allele count.
        dafs.append(v.genotypes.sum() / n_hap)
        times.append(mut.time)
        positions.append(v.site.position)

    tree = pd.DataFrame({
        'id': ids,
        'type': types,
        'selco': selcos,
        'daf': dafs,
        'time': times,
        'position': positions,
    })

    tree['selco'] = tree['selco'] / Q_scaling
    if not times_already_unscaled:
        if handoff_ticks is None:
            tree['time'] = tree['time'] * Q_scaling
        else:
            # Piecewise, for tree sequences spanning the Q=1 -> Q=0.01 handoff:
            # ticks younger than the handoff are worth Q_scaling real generations
            # each, older ones deep_Q_scaling. Diagnostic column only -- nothing
            # downstream reads it -- but a single factor here would misdate every
            # pre-handoff mutation by 100x.
            t = tree['time'].to_numpy(dtype=float)
            recent = np.minimum(t, handoff_ticks) * Q_scaling
            deep = np.maximum(t - handoff_ticks, 0.0) * deep_Q_scaling
            tree['time'] = recent + deep

    tree['maf'] = tree['daf'].apply(lambda x: min([x, 1-x]))
    tree['Vs'] = np.abs(tree['selco']) * tree['daf'] * (1-tree['daf'])
    
    return tree

def build_redistribution_map(donor_positions, recipient_positions, seed):
    # Pair each non-zero-selco donor position with a distinct neutral-selco
    # recipient. Sampling is seeded by `seed` for determinism. If fewer
    # recipients exist than donors, donors are randomly truncated so the
    # pairing is one-to-one. Donor keys are int (truncated) to support the
    # caller's `redistribution.get(int(pos))` lookup; recipient values are
    # kept as their original float position so the downstream pandas lookup
    # `vars[vars['position'] == recipient_position]` matches the float dtype
    # of `vars['position']`.
    rng = np.random.default_rng(seed)
    donors = sorted(int(p) for p in donor_positions)
    recipients = np.array(sorted(float(p) for p in recipient_positions), dtype=float)
    n_keep = min(len(donors), len(recipients))
    if n_keep < len(donors):
        kept_idx = rng.choice(len(donors), size=n_keep, replace=False)
        donors_kept = [donors[i] for i in sorted(kept_idx.tolist())]
    else:
        donors_kept = donors
    chosen_recipients = rng.choice(recipients, size=n_keep, replace=False)
    return donors_kept, dict(zip(donors_kept, [float(x) for x in chosen_recipients]))

def subsample_traits(df, n, seed):
    """Randomly keep `n` rows of `df` without replacement, seeded for
    determinism. Returns `df` unchanged if `n` is None or the pool already holds
    <= n rows (the "use all available" shortfall behavior). Kept rows stay in
    their original (position) order so downstream trait ordering is stable."""
    if n is None or df.shape[0] <= int(n):
        return df
    rng = np.random.default_rng(seed)
    keep = sorted(rng.choice(df.shape[0], size=int(n), replace=False).tolist())
    return df.iloc[keep]

def select_flank_gtex(gtex_vars, min_maf, flank_lo, flank_hi, n, seed):
    """GTEx-only trait loci drawn from the two edge buffers: positions
    <= flank_lo (region start .. 500 kb in) or >= flank_hi (500 kb before the
    region end .. end). Same causal eligibility as the central loci
    (selco != 0, MAF >= min_maf) but measured in the GTEx sample, since these
    traits are only phenotyped there. Randomly keeps `n` (seeded); uses all if
    fewer exist. These give GTEx loci >=500 kb from every central GWAS locus,
    for a non-colocalization comparison."""
    flank = gtex_vars.loc[
        (gtex_vars['maf'] >= min_maf) & (gtex_vars['selco'] != 0)
        & ((gtex_vars['position'] <= flank_lo) | (gtex_vars['position'] >= flank_hi))
    ]
    return subsample_traits(flank, n, seed)

def attach_beta(vars_df, key_df):
    # Add a 'beta' column to vars_df from the per-trait key DataFrame.
    # key_df['position'] is the recipient position (it's what tstrait wrote
    # the effect on), so a left merge places the signed beta on the
    # recipient row. Non-recipient rows get 0.
    key_lite = key_df[['position', 'effect_size']].rename(columns={'effect_size': 'beta'})
    key_lite['position'] = key_lite['position'].astype(float)
    out = vars_df.merge(key_lite, on='position', how='left')
    out['beta'] = out['beta'].fillna(0.0)
    return out

def generate_phenos_from_pos(position, ts, vars, scaling=1, recipient_position=None,
                             flip_seed=False, seed=0, site_by_pos=None, selco_by_pos=None):
    # Donor (`position`) determines the magnitude of beta via its selco;
    # recipient determines which site the effect lands on, the trait id, and the
    # seed used for `tstrait`'s random sign draw. When flip_seed=True (used for
    # GTEx) the sign draw is independent of the GWAS draw at the same recipient.
    #
    # `seed` is the run seed (--seed, i.e. stage2_seed). It is folded into both
    # the sign draw and the environmental noise so stage 2 is reproducible: the
    # noise used to come from the unseeded global np.random, which meant re-running
    # stage 2 with identical inputs produced different phenotypes.
    #
    # site_by_pos / selco_by_pos are optional {position: value} maps. Without them
    # each call did two linear scans of the whole vars frame.
    if recipient_position is None:
        recipient_position = position
    if site_by_pos is not None:
        site_id = int(site_by_pos[recipient_position])
        selco = float(selco_by_pos[position])
    else:
        site_id = int(vars[vars['position'] == recipient_position]['id'].iloc[0])
        selco = float(vars[vars['position'] == position]['selco'].iloc[0])
    # beta = selco_to_beta(selco)
    beta = np.sqrt(np.abs(selco)) * scaling
    recipient_position = int(recipient_position)

    # Independent, reproducible streams for the sign and the noise. Keying on
    # (run seed, recipient position, GWAS-vs-GTEx) preserves the two properties
    # the old scheme had -- different loci get different signs, and the same
    # locus gets independent signs in the GWAS and GTEx samples -- while adding
    # the run seed, which it lacked.
    sign_ss, noise_ss = np.random.SeedSequence(
        [int(seed), recipient_position, int(bool(flip_seed))]
    ).spawn(2)
    sign_seed = int(sign_ss.generate_state(1, dtype=np.uint32)[0])

    pheno = tstrait.sim_phenotype(ts, model=tstrait.TraitModelFixed(beta, random_sign=True), causal_sites=[site_id], h2=1, random_seed=sign_seed)
    pheno.phenotype['environmental_noise'] = np.random.default_rng(noise_ss).normal(
        0, 1, pheno.phenotype.shape[0])
    pheno.phenotype['trait_id'] = 'tr' + str(int(recipient_position))
    pheno.trait['trait_id'] = 'tr' + str(int(recipient_position))
    pheno.phenotype['phenotype'] = pheno.phenotype['genetic_value'] + pheno.phenotype['environmental_noise']

    return pheno

def combine_phenos_to_df(positions, ts, vars, scaling, redistribution=None, flip_seed=False, seed=0):
    id_list = ['tsk_' + str(x.id) for x in ts.individuals()]

    # Look-up maps built once instead of two linear scans of `vars` per trait.
    site_by_pos = dict(zip(vars['position'], vars['id']))
    selco_by_pos = dict(zip(vars['position'], vars['selco']))

    # Accumulate and concat once. Concatenating inside the loop rebuilt both
    # frames on every iteration, which is quadratic in the number of traits.
    pheno_cols = []
    trait_rows = []
    for pos in positions['position']:
        recipient = redistribution.get(int(pos)) if redistribution else None
        tr = generate_phenos_from_pos(pos, ts, vars, scaling, recipient_position=recipient,
                                      flip_seed=flip_seed, seed=seed,
                                      site_by_pos=site_by_pos, selco_by_pos=selco_by_pos)
        tr_id = tr.trait.trait_id.iloc[0]
        pheno_cols.append(tr.phenotype.phenotype.rename(tr_id).reset_index(drop=True))
        trait_rows.append(tr.trait)

    phenos = pd.concat([pd.DataFrame({'FID': 0, 'IID': id_list})] + pheno_cols, axis=1)
    snp_key = (pd.concat(trait_rows, ignore_index=True) if trait_rows else
               pd.DataFrame({'position': [], 'site_id': [], 'effect_size': [],
                             'causal_allele': [], 'allele_freq': [], 'trait_id': []}))

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

# create_pca() removed. It wrote 20 `controlled` SBAMS rows from ts.pca(20),
# which tskit only supports in mode="branch" -- i.e. eigenvectors of the genetic
# relatedness matrix of the local genealogy of the very 2 Mb region being
# fine-mapped. Those went to DAP-G as covariates while stage 5 separately
# regressed plink PCs out of the GLM. Neither arm has population structure to
# correct for (single panmictic populations; environmental noise drawn i.i.d.),
# so the PCs could only remove signal -- and in cattle, where 4*Ne*r*L across
# 2 Mb is of order 10-40, a handful of them span most of the genotype space.
# Measured on round-2 output: ~76-81% of the causal variant's genotype variance
# absorbed in cattle vs ~18% in human. See scripts/run_plink_glm.sh.

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
    commands.append(f'mkdir cgtex cgwas hgtex hgwas eo vcf')
    commands.append('mv ../tmp_upload/*vcf.gz vcf')
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
        commands.append(f'plink2 --pfile cgwas --maf 0.01 --glm hide-covar allow-no-covars --pheno cgwas_traits.scaling_{gwas_scaling}.pheno --out glm/cgwas_glm_scaling_{gwas_scaling}')
        commands.append('# cgtex')
        commands.append('plink2 --vcf cgtex.vcf --make-pgen --out cgtex')
        commands.append(f'plink2 --pfile cgtex --maf 0.01 --glm hide-covar allow-no-covars --pheno cgtex_traits.scaling_{gtex_scaling}.pheno --out glm/cgtex_glm_scaling_{gtex_scaling}')
    if run_hums:
        commands.append('# hgwas')
        commands.append('plink2 --vcf hgwas.vcf --make-pgen --out hgwas')
        commands.append(f'plink2 --pfile hgwas --maf 0.01 --glm hide-covar allow-no-covars --pheno hgwas_traits.scaling_{gwas_scaling}.pheno --out glm/hgwas_glm_scaling_{gwas_scaling}')
        commands.append('# hgtex')
        commands.append('plink2 --vcf hgtex.vcf --make-pgen --out hgtex')
        commands.append(f'plink2 --pfile hgtex --maf 0.01 --glm hide-covar allow-no-covars --pheno hgtex_traits.scaling_{gtex_scaling}.pheno --out glm/hgtex_glm_scaling_{gtex_scaling}')
    command_string = '\n'.join(commands)
    return command_string

def get_arguments():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--gwas_scaling', type=int)
    parser.add_argument('--gtex_scaling', type=int)
    parser.add_argument('--r2_value', type=float) # Note to self: I think I no longer use this. Confirm, then delete.
    parser.add_argument('--Q_scaling', type=float, required=True, help='The stage-1 SLiM rescaling factor (config key Q_scaling; cattle 0.01, human 10). Selection coefficients recorded in the tree sequence are Q_scaling * the real value and are divided back by get_vars_df, which sets beta = sqrt(|selco|) * scaling. Applies to whichever single species this invocation runs.')
    parser.add_argument('--handoff_ticks', type=float, required=False, default=None,
                        help='Cattle deep-history handoff: number of ticks between the '
                             'Q=1 -> Q_scaling switch and the end of the simulation '
                             '(epochs 8-12 = 2400 at Q_scaling=0.01). Set together with '
                             '--deep_Q_scaling when stage 1 used the split-Q architecture. '
                             'Such a tree sequence has a piecewise time scale, so the '
                             'neutral overlay and the time column must be applied in two '
                             'segments; a single rate would lose almost all of the deep '
                             'neutral variation. Omit for a single-Q stage 1.')
    parser.add_argument('--deep_Q_scaling', type=float, required=False, default=None,
                        help='Q_scaling the pre-handoff (deep history) phase ran under; 1 '
                             'for the split-Q architecture. Required with --handoff_ticks.')
    parser.add_argument('--min_maf', type=float, help='The minor allele frequency variants must reach to be tested')
    parser.add_argument('--length', type=float, required=False, default=1e7, help='Genomic region length L (bp). Trait-associated variants are restricted to [5e5, L-5e5] (a 500 kb buffer from each edge). Default 1e7 reproduces the legacy hardcoded [5e5, 9.5e6] window.')
    parser.add_argument('--human_ts_file', type=str)
    parser.add_argument('--cattle_ts_file', type=str)
    parser.add_argument('--cattle_m4_file', type=str, help='A file keeping track of the cattle variants under intense agricultural selection')
    parser.add_argument('--out_dir', type=str)
    parser.add_argument('--n_samples', type=float, required=False, default=None, help='The number of individuals in the GTEx sample and GWAS sample combined (default is the full number of simulated individuals)')
    parser.add_argument('--gtex_size', type=float, required=False, default=500, help='The number of individuals in the GTEx sample (default 500; -1 creates outputs of 1k, 500, and 250)')
    parser.add_argument('--already_includes_neutral', type=str2bool, required=False, default=False, help='Only needed if you already added neutral mutations in the input data')
    parser.add_argument('--seed', type=int, default=19930224)
    parser.add_argument('--neutral_trait_vars', type=str2bool, required=False, default=False, help='Redistribute each non-zero effect size from its causal (selco != 0) donor variant to a random eligible neutral (selco == 0) recipient. Trait IDs are named for the recipient position. The redistribution is identical across GWAS and GTEx samples and is seeded by --seed.')
    parser.add_argument('--n_central_traits', type=int, required=False, default=None, help='Number of central trait loci to keep -- these are the GWAS traits AND the shared GTEx traits (same positions). Drawn from the eligible causal pool (selco != 0, MAF >= min_maf, central [5e5, L-5e5] window). Randomly subsampled (seeded by --seed) when the pool is larger; all are used when fewer exist. Default None = use all eligible (legacy behavior).')
    parser.add_argument('--n_flank_gtex_traits', type=int, required=False, default=50, help='Number of GTEx-ONLY trait loci drawn from the two 500 kb edge flanks ([0,5e5] U [L-5e5,L]) with the same causal eligibility (selco != 0, MAF >= min_maf in the GTEx sample). Added to the GTEx trait set on top of the shared central loci, so each GWAS locus can be compared to GTEx loci >=500 kb away. Randomly subsampled (seeded); uses all if fewer exist. Default 50; set 0 to disable.')
    return parser.parse_args()

if __name__ == '__main__':
    print('Getting arguments')
    args = get_arguments()
    gwas_scaling = args.gwas_scaling
    gtex_scaling = args.gtex_scaling
    r2_value = args.r2_value
    min_maf = args.min_maf
    gtex_size = args.gtex_size

    # Trait-associated (causal) variants are restricted to the central region,
    # excluding a 500 kb buffer from each edge. Computed from L so a smaller
    # region gets the correct central window (L=2 Mb -> [5e5, 1.5e6], a 1 Mb
    # center). At the default L=1e7 this is [5e5, 9.5e6] -- identical to the
    # previous hardcoded bounds, so 10 Mb runs are unaffected.
    region_length = args.length
    trait_pos_lo = 5e5
    trait_pos_hi = region_length - 5e5

    if gtex_size < 0:
        multiple_gtex_outputs = True
        gtex_size = 1000
    else:
        multiple_gtex_outputs = False

    run_hums = not args.human_ts_file is None
    run_cows = not args.cattle_ts_file is None

    # A single --Q_scaling cannot serve two species with different stage-1
    # rescaling factors (cattle 0.01, human 10). The Snakemake stage-2 rule
    # always passes exactly one ts file (rules/common.smk: params.ts_flag), so
    # this only fires for a hand-written dual-species invocation; fail loudly
    # rather than silently mis-scaling one arm's effect sizes.
    if run_hums and run_cows:
        raise SystemExit(
            'Refusing to run both species in one invocation: --Q_scaling '
            f'({args.Q_scaling}) applies to a single species. Run cattle and '
            'human separately, each with its own --Q_scaling.'
        )
    q_scaling = args.Q_scaling

    # Clean up the data
    print('Loading data')
    if run_cows:
        # cows = tskit.load('/Users/noah/comparative_colocalization/data/simulations/demographic_simulations/cattle/farm_selection_Q_0.01.L_10000000.seed_20250303.full.ts')
        cows = tskit.load(args.cattle_ts_file)
        if args.cattle_m4_file is not None:
            cows = relabel_ag_variants(cows, args.cattle_m4_file)
        print('Cattle demography')
        # add_neutral's Q is the RECIPROCAL of the config's Q_scaling: it
        # divides the overlay rate (mut_rate = 8.4e-9 / Q), and raw SLiM
        # branch lengths are inflated by 1/Q_scaling, so the rate must shrink
        # by the same factor. 1/0.01 == 100.0 exactly, so this is bit-identical
        # to the literal 100 it replaces.
        cows = remove_fixed(add_neutral(cows, Q=1 / q_scaling, seed=args.seed,
                                        handoff_ticks=args.handoff_ticks,
                                        deep_Q_scaling=args.deep_Q_scaling))
        # generate_nucleotides defaults to seed=None, i.e. a fresh random draw
        # each run. Only the REF/ALT letters depend on it -- genotypes, positions
        # and phenotypes are unaffected -- but leaving it unseeded meant two runs
        # with identical inputs produced byte-different VCFs.
        cows = pyslim.convert_alleles(pyslim.generate_nucleotides(cows, seed=args.seed))
        cows = remove_position_zero(cows)
    if run_hums:
        # hums = tskit.load('/Users/noah/comparative_colocalization/data/simulations/demographic_simulations/human/human_selection_Q_1.L10000000.seed_20250521.full.ts')
        hums = tskit.load(args.human_ts_file)
        print('Human demography')
        if not args.already_includes_neutral:
            # Deliberately left at the default Q=1, i.e. the unscaled 8.4e-9.
            # stdpopsim's _recap_and_rescale has already put times back into
            # real generations, so unlike the cattle branch there is nothing to
            # compensate for here. (Dead path today -- every config sets
            # already_includes_neutral: True. Stage 1's own overlay in
            # human_simulation_o2.py used to carry a `* Q` here, the same
            # conflation that was fixed for selco; it is now plain 8.4e-9 and
            # agrees with this branch.)
            hums = remove_fixed(add_neutral(hums, seed=args.seed))
        # tskit's write_vcf emits REF="" / ALT="<int>" for SLiM-origin mutations
        # (SLiMMutationModel uses integer alleles, no ancestral_state string).
        # plink2 --make-pgen then writes a malformed pvar and crashes. Apply the
        # nucleotide overlay unconditionally — mirrors what the cattle branch
        # does above, and adds nothing when already_includes_neutral=True.
        hums = pyslim.convert_alleles(pyslim.generate_nucleotides(hums, seed=args.seed))
        hums = remove_position_zero(hums)

    # Split into GWAS and GTEx
    print('Splitting into GWAS and GTEx')
    # --n_samples is parsed as float; it feeds range()/stride math below (and
    # when unset defaults to an int num_individuals), so coerce to int here or
    # `range(0, 2*n_samples, ...)` raises TypeError on a float.
    n_samples = args.n_samples if args.n_samples is None else int(args.n_samples)
    hgtex_subsamples = {}
    cgtex_subsamples = {}
    if run_hums:
        if n_samples is None:
            n_samples = hums.num_individuals
        print('Human samples: ' + str(n_samples))
        samps = [i for i in range(0, 2*n_samples, int(n_samples/(gtex_size/2)))] + [i+1 for i in range(0, 2*n_samples, int(n_samples/(gtex_size/2)))]
        samps.sort()
        # Set for the membership test below; `samps` itself stays a sorted list
        # because it indexes ts.samples(). As a list the test was O(N^2)
        # (18,000 x 2,000 comparisons per species).
        samps_set = set(samps)
        hgtex = hums.simplify(hums.samples()[samps])
        hgwas = hums.simplify(hums.samples()[[i for i in range(0, 2*n_samples) if i not in samps_set]])
        if multiple_gtex_outputs:
            # Sub-sample 500 individuals into hgtex_smaller and 250 into hgtex_smallest from
            # the 1000-individual hgtex. Same paired-haplotype stride pattern as the gwas/gtex
            # split (pull haplotypes (i, i+1) at a stride determined by the target sample size)
            # so each retained individual still has both chromosomes.
            for sub_size, sub_name in [(500, 'hgtex_smaller'), (250, 'hgtex_smallest')]:
                step = int(1000 / (sub_size / 2))
                sub_samps = [i for i in range(0, 2*1000, step)] + [i+1 for i in range(0, 2*1000, step)]
                sub_samps.sort()
                hgtex_subsamples[sub_name] = {
                    'ts': hgtex.simplify(hgtex.samples()[sub_samps]),
                    'retained_inds': [s // 2 for s in sub_samps[::2]],
                }
    if run_cows:
        if n_samples is None:
            n_samples = cows.num_individuals
        print('Cattle samples: ' + str(n_samples))
        samps = [i for i in range(0, 2*n_samples, int(n_samples/(gtex_size/2)))] + [i+1 for i in range(0, 2*n_samples, int(n_samples/(gtex_size/2)))]
        samps.sort()
        samps_set = set(samps)  # see the human branch above
        cgtex = cows.simplify(cows.samples()[samps])
        cgwas = cows.simplify(cows.samples()[[i for i in range(0, 2*n_samples) if i not in samps_set]])
        if multiple_gtex_outputs:
            for sub_size, sub_name in [(500, 'cgtex_smaller'), (250, 'cgtex_smallest')]:
                step = int(1000 / (sub_size / 2))
                sub_samps = [i for i in range(0, 2*1000, step)] + [i+1 for i in range(0, 2*1000, step)]
                sub_samps.sort()
                cgtex_subsamples[sub_name] = {
                    'ts': cgtex.simplify(cgtex.samples()[sub_samps]),
                    'retained_inds': [s // 2 for s in sub_samps[::2]],
                }

    # Get the variants
    print('Getting variants')
    # Cattle .ts files come straight from SLiM (sim.treeSeqOutput), so neither
    # the selection coefficients nor the times have been un-rescaled. Human
    # .ts files come from stdpopsim, which un-rescales times but not selection
    # coefficients -- hence times_already_unscaled=True on that side only.
    if run_cows:
        print(f'Undoing Q_scaling={q_scaling} on cattle selection coefficients '
              f'(selco / {q_scaling}) and times (time * {q_scaling})')
        time_kw = dict(handoff_ticks=args.handoff_ticks, deep_Q_scaling=args.deep_Q_scaling)
        cvars = get_vars_df(cows, Q_scaling=q_scaling, **time_kw)
        cgwas_vars = get_vars_df(cgwas, Q_scaling=q_scaling, **time_kw)
        cgtex_vars = get_vars_df(cgtex, Q_scaling=q_scaling, **time_kw)
        print('Number of cattle vars: ' + str(len(cvars)))
        print('Median |selco| at selected sites: '
              f"{np.abs(cvars.loc[cvars['selco'] != 0, 'selco']).median():.4e}")
        for name, data in cgtex_subsamples.items():
            data['vars'] = get_vars_df(data['ts'], Q_scaling=q_scaling, **time_kw)
            print(f'Number of {name} vars: ' + str(len(data['vars'])))
    if run_hums:
        print(f'Undoing Q_scaling={q_scaling} on human selection coefficients '
              f'(selco / {q_scaling}); stdpopsim already un-rescaled the times')
        hvars = get_vars_df(hums, Q_scaling=q_scaling, times_already_unscaled=True)
        hgwas_vars = get_vars_df(hgwas, Q_scaling=q_scaling, times_already_unscaled=True)
        hgtex_vars = get_vars_df(hgtex, Q_scaling=q_scaling, times_already_unscaled=True)
        print('Number of human vars: ' + str(len(hvars)))
        print('Median |selco| at selected sites: '
              f"{np.abs(hvars.loc[hvars['selco'] != 0, 'selco']).median():.4e}")
        for name, data in hgtex_subsamples.items():
            data['vars'] = get_vars_df(data['ts'], Q_scaling=q_scaling,
                                       times_already_unscaled=True)
            print(f'Number of {name} vars: ' + str(len(data['vars'])))

    # Select the relevant alleles
    print('Selecting relevant alleles')
    if run_hums:
        hgtex_maf_dict = {var.position: var.maf for _, var in hgtex_vars.iterrows()}
        hcausative_maf01 = hgwas_vars.loc[(hgwas_vars['maf'] >= min_maf) & (hgwas_vars['selco'] != 0) & (hgwas_vars['position'] > trait_pos_lo) & (hgwas_vars['position'] < trait_pos_hi) & (hgwas_vars.position.isin(hgtex_vars.position))]
        print(f'Number of causative human variants: {hcausative_maf01.shape[0]}')
        hcausative_maf01 = hcausative_maf01[hcausative_maf01.apply(lambda row: hgtex_maf_dict[row.position] >= min_maf, axis=1)]
    if run_cows:
        cgtex_maf_dict = {var.position: var.maf for _, var in cgtex_vars.iterrows()}
        ccausative_maf01 = cgwas_vars.loc[(cgwas_vars['maf'] >= min_maf) & (cgwas_vars['selco'] != 0) & (cgwas_vars['position'] > trait_pos_lo) & (cgwas_vars['position'] < trait_pos_hi) & (cgwas_vars.position.isin(cgtex_vars.position))]
        print(f'Number of causative cattle variants: {ccausative_maf01.shape[0]}')
        ccausative_maf01 = ccausative_maf01[ccausative_maf01.apply(lambda row: cgtex_maf_dict[row.position] >= min_maf, axis=1)]

    # If --neutral_trait_vars is set, redistribute each non-zero effect from
    # its selco != 0 donor to a random selco == 0 recipient. Recipients are
    # filtered identically to donors (MAF in GWAS + position window + GTEx
    # cross-check). With fewer recipients than donors, donors are randomly
    # truncated (seeded by --seed) so the pairing stays one-to-one.
    hredist = credist = None
    if args.neutral_trait_vars:
        if run_hums:
            hneutral_maf01 = hgwas_vars.loc[(hgwas_vars['maf'] >= min_maf) & (hgwas_vars['selco'] == 0) & (hgwas_vars['position'] > trait_pos_lo) & (hgwas_vars['position'] < trait_pos_hi) & (hgwas_vars.position.isin(hgtex_vars.position))]
            hneutral_maf01 = hneutral_maf01[hneutral_maf01.apply(lambda row: hgtex_maf_dict[row.position] >= min_maf, axis=1)]
            hkept, hredist = build_redistribution_map(hcausative_maf01['position'], hneutral_maf01['position'], seed=args.seed)
            print(f'Human donors: {len(hcausative_maf01)}, neutral recipients: {len(hneutral_maf01)}, paired: {len(hkept)}')
            hcausative_maf01 = hcausative_maf01[hcausative_maf01['position'].astype(int).isin(hkept)]
        if run_cows:
            cneutral_maf01 = cgwas_vars.loc[(cgwas_vars['maf'] >= min_maf) & (cgwas_vars['selco'] == 0) & (cgwas_vars['position'] > trait_pos_lo) & (cgwas_vars['position'] < trait_pos_hi) & (cgwas_vars.position.isin(cgtex_vars.position))]
            cneutral_maf01 = cneutral_maf01[cneutral_maf01.apply(lambda row: cgtex_maf_dict[row.position] >= min_maf, axis=1)]
            ckept, credist = build_redistribution_map(ccausative_maf01['position'], cneutral_maf01['position'], seed=args.seed)
            print(f'Cattle donors: {len(ccausative_maf01)}, neutral recipients: {len(cneutral_maf01)}, paired: {len(ckept)}')
            ccausative_maf01 = ccausative_maf01[ccausative_maf01['position'].astype(int).isin(ckept)]

    # Down-select trait loci.
    #  - GWAS traits = central causal loci, subsampled to --n_central_traits
    #    (used for BOTH the GWAS sample and the shared GTEx loci, so the shared
    #    set is identical between the two).
    #  - GTEx traits = those shared central loci PLUS --n_flank_gtex_traits
    #    GTEx-only loci from the 500 kb edge flanks.
    # Distinct seed offsets keep the central and flank draws (and the two
    # species) independent but deterministic given --seed (= stage2_seed).
    print('Selecting trait loci')
    if run_hums:
        hcausative_maf01 = subsample_traits(hcausative_maf01, args.n_central_traits, args.seed)
        hflank_gtex = select_flank_gtex(hgtex_vars, min_maf, trait_pos_lo, trait_pos_hi,
                                        args.n_flank_gtex_traits, args.seed + 10**7)
        hgtex_trait_pos = pd.concat([hcausative_maf01, hflank_gtex])
        print(f'Human trait loci -- GWAS/central: {hcausative_maf01.shape[0]}, '
              f'GTEx flank: {hflank_gtex.shape[0]}, GTEx total: {hgtex_trait_pos.shape[0]}')
    if run_cows:
        ccausative_maf01 = subsample_traits(ccausative_maf01, args.n_central_traits, args.seed + 2*10**7)
        cflank_gtex = select_flank_gtex(cgtex_vars, min_maf, trait_pos_lo, trait_pos_hi,
                                        args.n_flank_gtex_traits, args.seed + 3*10**7)
        cgtex_trait_pos = pd.concat([ccausative_maf01, cflank_gtex])
        print(f'Cattle trait loci -- GWAS/central: {ccausative_maf01.shape[0]}, '
              f'GTEx flank: {cflank_gtex.shape[0]}, GTEx total: {cgtex_trait_pos.shape[0]}')

    # Creating phenotypes
    print('Creating phenotypes')
    if run_cows:
        cgwas_key_maf01, cgwas_traits_maf01 = combine_phenos_to_df(ccausative_maf01, cgwas, cgwas_vars, scaling=gwas_scaling, redistribution=credist, seed=args.seed)
        cgtex_key_maf01, cgtex_traits_maf01 = combine_phenos_to_df(cgtex_trait_pos, cgtex, cgtex_vars, scaling=gtex_scaling, redistribution=credist, flip_seed=True, seed=args.seed)
        for name, data in cgtex_subsamples.items():
            sub_traits = cgtex_traits_maf01.iloc[data['retained_inds']].copy().reset_index(drop=True)
            sub_traits['IID'] = ['tsk_' + str(i) for i in range(len(sub_traits))]
            data['traits'] = sub_traits
    if run_hums:
        hgwas_key_maf01, hgwas_traits_maf01 = combine_phenos_to_df(hcausative_maf01, hgwas, hgwas_vars, scaling=gwas_scaling, redistribution=hredist, seed=args.seed)
        hgtex_key_maf01, hgtex_traits_maf01 = combine_phenos_to_df(hgtex_trait_pos, hgtex, hgtex_vars, scaling=gtex_scaling, redistribution=hredist, flip_seed=True, seed=args.seed)
        # For subsamples, take the simulated trait rows of the retained individuals.
        # Same simulated outcomes, just observed in fewer individuals (power-analysis setup).
        # Renumber IIDs so they match the re-indexed individuals in each subsampled VCF (tsk_0..tsk_N-1).
        for name, data in hgtex_subsamples.items():
            sub_traits = hgtex_traits_maf01.iloc[data['retained_inds']].copy().reset_index(drop=True)
            sub_traits['IID'] = ['tsk_' + str(i) for i in range(len(sub_traits))]
            data['traits'] = sub_traits

    # Write everything out
    print('Writing everything out')
    # sim_dir = f'/Users/noah/comparative_colocalization/data/simulations/gwas_and_eqtl_mapping/gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}'
    sim_dir = f'{args.out_dir}/gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}'
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    if run_cows:
        attach_beta(cgwas_vars, cgwas_key_maf01).to_csv(f'{sim_dir}/cgwas_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)
        attach_beta(cgtex_vars, cgtex_key_maf01).to_csv(f'{sim_dir}/cgtex_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)
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

        # And we're also going to create the output for plink -> sumstats -> coloc
        cgwas_traits_maf01.to_csv(f'{sim_dir}/cgwas_traits.scaling_{gwas_scaling}.pheno', sep='\t', index=False)
        cgtex_traits_maf01.to_csv(f'{sim_dir}/cgtex_traits.scaling_{gtex_scaling}.pheno', sep='\t', index=False)

        for name, data in cgtex_subsamples.items():
            sub_ts, sub_traits, sub_vars = data['ts'], data['traits'], data['vars']
            attach_beta(sub_vars, cgtex_key_maf01).to_csv(f'{sim_dir}/{name}_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)
            sub_traits.to_csv(f'{sim_dir}/{name}_traits_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)
            with open(f'{sim_dir}/{name}.vcf', 'w') as f:
                sub_ts.write_vcf(f)
            write_ts_as_sbams(sub_ts, f'{sim_dir}/{name}_scaling_{gtex_scaling}_geno.sbams')
            write_traits_as_sbams(sub_traits, f'{sim_dir}/{name}_scaling_{gtex_scaling}_pheno.sbams')
            sub_traits.to_csv(f'{sim_dir}/{name}_traits.scaling_{gtex_scaling}.pheno', sep='\t', index=False)

    if run_hums:
        attach_beta(hgwas_vars, hgwas_key_maf01).to_csv(f'{sim_dir}/hgwas_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)
        attach_beta(hgtex_vars, hgtex_key_maf01).to_csv(f'{sim_dir}/hgtex_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)
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

        hgwas_traits_maf01.to_csv(f'{sim_dir}/hgwas_traits.scaling_{gwas_scaling}.pheno', sep='\t', index=False)
        hgtex_traits_maf01.to_csv(f'{sim_dir}/hgtex_traits.scaling_{gtex_scaling}.pheno', sep='\t', index=False)

        for name, data in hgtex_subsamples.items():
            sub_ts, sub_traits, sub_vars = data['ts'], data['traits'], data['vars']
            attach_beta(sub_vars, hgtex_key_maf01).to_csv(f'{sim_dir}/{name}_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)
            sub_traits.to_csv(f'{sim_dir}/{name}_traits_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.tsv', sep='\t', index=False)
            with gzip.open(f'{sim_dir}/{name}.vcf.gz', 'wt') as f:
                sub_ts.write_vcf(f)
            write_ts_as_sbams(sub_ts, f'{sim_dir}/{name}_scaling_{gtex_scaling}_geno.sbams')
            write_traits_as_sbams(sub_traits, f'{sim_dir}/{name}_scaling_{gtex_scaling}_pheno.sbams')
            sub_traits.to_csv(f'{sim_dir}/{name}_traits.scaling_{gtex_scaling}.pheno', sep='\t', index=False)

    # Write out a long string to the file sim_dir/gwas_scaling_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.txt
    with open(f'{sim_dir}/gwas_scaling_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}.txt', 'w') as f:
        f.write(get_commands(gwas_scaling, gtex_scaling, min_maf, sim_dir, run_hums, run_cows))



# Included in final analysis

# Not included in final analysis

# Too soon to tell
# python ../create_gwas_files_and_phenotypes.py --gwas_scaling 35 --gtex_scaling 35 --r2_value 0.2 --min_maf 0.01 --cattle_ts_file '/Users/noah/tmp/selsims/revision_farm_selection_mult_100_gen_23_muts_26_bottlenecked_sd24.full.ts' --cattle_m4_file '/Users/noah/tmp/selsims/revision_farm_selection_mult_100_gen_23_muts_26_bottlenecked_sd24.m4_marks.tsv' --out_dir '/Users/noah/tmp/selsims/outdir/'
# python ../create_gwas_files_and_phenotypes.py --gtex_size 50000 --gwas_scaling 35 --gtex_scaling 35 --r2_value 0.2 --min_maf 0.01 --human_ts_file '/Users/noah/tmp/selsims/largesim/hts_19930224_large.ts' --out_dir '/Users/noah/tmp/selsims/largesim/outputs/' --already_includes_neutral True
