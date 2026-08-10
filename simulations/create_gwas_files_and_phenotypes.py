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

# sys.path[0] is this script's directory, so `helpers` resolves however the script
# is invoked -- by absolute path from a Snakemake rule, or by hand from anywhere.
from helpers import causal_power, params_record, paths, synthetic_dfe


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

def causal_eligible(vars_df, pos_lo, pos_hi, causal_min_maf, synthetic=False):
    """Central-window variants that may be drawn as causal, in EITHER sampling scheme.

    One predicate for both, which is the point: `causal_min_maf` used to gate the
    uniform pool and be ignored by the power pool, on the grounds that the
    detection-power weight subsumed a frequency floor. It does not -- a weight is a
    soft preference and a floor is a hard cut, and having one knob mean two things
    depending on another knob made "the pool at causal_min_maf=0.01" ambiguous.
    They now compose.

    `maf > 0` is not the floor, it is well-formedness, and it survives
    `causal_min_maf = 0` (the default). `remove_fixed` runs on the full tree
    sequence BEFORE the GWAS/GTEx split, so a site can still be monomorphic within
    one panel; there it is not a variant at all, its power is undefined rather than
    small, and an effect assigned to it would produce a phenotype of pure noise.
    The uniform path never had to say this because a non-zero causal_min_maf
    excluded those sites incidentally.

    `synthetic` drops the `selco != 0` requirement, for the neutral simulations
    (category H) where no variant carries a selection coefficient and the causal
    set is drawn from all of them -- see --synthetic_dfe_effects.

    NOTE this is the GWAS-panel pool only. The uniform path additionally
    intersects with the GTEx panel at the call site; the power path deliberately
    does not.
    """
    keep = ((vars_df['maf'] > 0)
            & (vars_df['maf'] >= causal_min_maf)
            & (vars_df['position'] > pos_lo)
            & (vars_df['position'] < pos_hi))
    if not synthetic:
        keep &= (vars_df['selco'] != 0)
    return vars_df.loc[keep]


def select_central_power(pool, n, scaling, sampling_n, sig_p, min_power,
                         min_pool_multiple, seed, maf_by_position=None, label=''):
    """Draw `n` causal loci with inclusion probability proportional to GWAS power.

    Each candidate's weight is its probability of reaching `sig_p` in a GWAS of
    `sampling_n` individuals, given beta = sqrt(|selco|) * `scaling` and its MAF
    in the GWAS panel -- see helpers/causal_power.py for the model. The draw is
    pi-PS (exact first-order inclusion probabilities), NOT numpy's successive
    weighted choice, which would over-represent the high-power tail relative to
    its own weights.

    `maf_by_position` overrides the MAF used for the weight, keyed by the DONOR's
    position. It is how --neutral_trait_vars is handled: there the effect is moved
    to a neutral recipient, so the power a locus actually has is set by the
    recipient's frequency and the donor's beta, and weighting on the donor's own
    MAF would rank the pool by a quantity no downstream test ever sees.

    Raises SystemExit when too little of the pool is detectable at all. A draw
    from a pool where only a handful of variants carry real weight is not a
    weighted sample -- it is those few plus an arbitrary tail -- and silently
    returning it would look exactly like a successful run.

    Returns (chosen_rows, diagnostics) where diagnostics covers the WHOLE pool,
    one row per candidate, so the draw is auditable after the fact."""
    if maf_by_position is None:
        maf = pool['maf'].to_numpy(dtype=float)
    else:
        maf = np.array([maf_by_position[p] for p in pool['position']], dtype=float)
    beta = np.sqrt(np.abs(pool['selco'].to_numpy(dtype=float))) * scaling
    power = causal_power.detection_power(maf, beta, sampling_n, sig_p)

    n_above = int((power >= min_power).sum())
    required = int(np.ceil(min_pool_multiple * n))
    if n_above < required:
        raise SystemExit(
            f'Refusing to draw {n} {label} causal variants by power: only {n_above} of '
            f'{pool.shape[0]} pool variants reach power >= {min_power} at '
            f'sampling_gwas_n={sampling_n:g} and sig_p={sig_p:g}, but '
            f'{required} are required (min_pool_multiple={min_pool_multiple:g}). '
            'A draw this concentrated is not a weighted sample. Lower '
            '--sampling_min_power or --sampling_min_pool_multiple, raise '
            '--sampling_gwas_n or the phenotype scaling, or loosen --sampling_sig_p.'
        )

    pi = causal_power.inclusion_probabilities(power, n)
    idx = causal_power.systematic_sample(pi, np.random.default_rng(seed))

    diagnostics = pd.DataFrame({
        'position': pool['position'].to_numpy(),
        'selco': pool['selco'].to_numpy(),
        'maf': pool['maf'].to_numpy(),
        'maf_for_power': maf,
        'beta': beta,
        'power': power,
        'pi': pi,
        'selected': np.isin(np.arange(pool.shape[0]), idx),
    })
    return pool.iloc[idx], diagnostics


def select_gtex_topup(gtex_vars, chosen, pos_lo, pos_hi, n, seed, redistribution=None,
                      causal_min_maf=0.0, synthetic=False):
    """The GTEx central trait set under --causal_sampling power.

    The GWAS causal variants are no longer required to segregate in the GTEx panel
    -- the uniform path guarantees that by intersecting the eligible pool with
    `gtex_vars.position`, and this one does not -- so the shared set is whichever
    of them do. The remaining slots are filled UNIFORMLY from central selco != 0
    GTEx variants that are not already trait loci, up to `n` in total.

    Every membership test here is on the TRAIT position -- the redistribution
    recipient under --neutral_trait_vars, the variant's own position otherwise --
    never on the row's position, because the trait position is what
    combine_phenos_to_df resolves a site id for and what names the trait. Two
    consequences, both load-bearing:

      * a locus whose trait position the GTEx panel lacks cannot be shared at all
        (the site lookup would KeyError), and the same test gates top-up
        candidates; and
      * a top-up candidate is excluded when its trait position is already taken.
        Under redistribution that is not the same as excluding the drawn rows: a
        drawn locus's DONOR position is still an ordinary GTEx variant, but
        selecting it would remap it onto the recipient the drawn locus already
        occupies and emit two GTEx traits with one name.

    Returns (shared, topup). `shared` stays as rows of `chosen` -- the donor frame
    -- so that the same `redistribution` map remaps it on the GTEx side exactly as
    it does on the GWAS side. Downstream, combine_phenos_to_df reads two columns
    off these rows: `position` (the donor, remapped to the trait position) and
    `selco` (the donor's, which sets beta). Both are donor properties and both are
    correct here; what is NOT safe is resolving either against the GTEx frame,
    since a donor drawn by power need not segregate in the GTEx panel.
    """
    def trait_pos(p):
        return redistribution.get(int(p), p) if redistribution else p

    # "Present in the GTEx panel" means SEGREGATING there. A site the panel carries
    # but is monomorphic for cannot be an eQTL: the effect would land on a column of
    # identical genotypes and the trait would be pure noise. The uniform path got
    # this for free from its causal_min_maf cross-check.
    gtex_pos = set(gtex_vars.loc[gtex_vars['maf'] > 0, 'position'])
    chosen_trait_pos = {trait_pos(p) for p in chosen['position']}
    shared = chosen.loc[[trait_pos(p) in gtex_pos for p in chosen['position']]]

    n_missing = max(int(n) - shared.shape[0], 0)
    # Same eligibility rule the GWAS pool used, measured in this panel. It gained
    # causal_min_maf when the floor stopped being scheme-dependent: a top-up locus
    # is as causative as a drawn one, so a floor that excluded the latter and not
    # the former would put variants in the GTEx causal set that the config forbade.
    candidates = causal_eligible(gtex_vars, pos_lo, pos_hi, causal_min_maf,
                                 synthetic=synthetic)
    cand_trait_pos = [trait_pos(p) for p in candidates['position']]
    eligible = [t in gtex_pos and t not in chosen_trait_pos for t in cand_trait_pos]
    topup_pool = candidates.loc[eligible]
    topup = subsample_traits(topup_pool, n_missing, seed) if n_missing else topup_pool.iloc[:0]
    return shared, topup


def trait_partner_table(gwas_chosen, gtex_shared, redistribution=None):
    """One row per GWAS causal locus: its trait id, and its GTEx partner if it has one.

    Under uniform sampling every GWAS locus has a partner by construction, so this
    is an invariant check. Under power sampling it is the answer key: both scorers
    -- helpers/summarize_coloc.py and figures_and_tables/figure2_revision2.ipynb --
    currently infer the pairing from trait-name equality alone, which cannot tell
    "this GWAS locus has no eQTL to find" apart from "colocalization failed", so
    this table is what a corrected scorer needs to read instead."""
    def trait_id(p):
        eff = redistribution.get(int(p), p) if redistribution else p
        return 'tr' + str(int(eff))

    shared_pos = set(gtex_shared['position'])
    positions = [int(p) for p in gwas_chosen['position']]
    shared = [p in shared_pos for p in gwas_chosen['position']]
    ids = [trait_id(p) for p in gwas_chosen['position']]
    return pd.DataFrame({
        'gwas_trait': ids,
        'gtex_trait': [i if s else '' for i, s in zip(ids, shared)],
        'position': positions,
        'shared': shared,
    })


def flank_eligible(gtex_vars, causal_min_maf, flank_lo, flank_hi, synthetic=False):
    """The GTEx-only flank pool, before the `n` subsample. Split out from
    select_flank_gtex so the synthetic-DFE draw can cover the same candidate set
    the subsample will pick from -- see helpers/synthetic_dfe.assign_by_position,
    whose determinism depends on being handed the whole universe of positions
    rather than the ones that happened to be drawn."""
    keep = ((gtex_vars['maf'] > 0)
            & (gtex_vars['maf'] >= causal_min_maf)
            & ((gtex_vars['position'] <= flank_lo) | (gtex_vars['position'] >= flank_hi)))
    if not synthetic:
        keep &= (gtex_vars['selco'] != 0)
    return gtex_vars.loc[keep]


def select_flank_gtex(gtex_vars, causal_min_maf, flank_lo, flank_hi, n, seed,
                      synthetic=False):
    """GTEx-only trait loci drawn from the two edge buffers: positions
    <= flank_lo (region start .. 500 kb in) or >= flank_hi (500 kb before the
    region end .. end). Same causal eligibility as the central loci
    (selco != 0, MAF >= causal_min_maf) but measured in the GTEx sample, since
    these traits are only phenotyped there. Randomly keeps `n` (seeded); uses all
    if fewer exist. These give GTEx loci >=500 kb from every central GWAS locus,
    for a non-colocalization comparison.

    The floor is causal_min_maf, not min_maf: these loci ARE causative, so they
    are drawn from the same pool as the central ones. See --causal_min_maf.
    `maf > 0` accompanies the floor for the same reason it does in
    causal_eligible: at causal_min_maf = 0 the floor alone would admit sites
    monomorphic in this panel.

    `synthetic` drops the `selco != 0` requirement -- see causal_eligible."""
    return subsample_traits(
        flank_eligible(gtex_vars, causal_min_maf, flank_lo, flank_hi, synthetic),
        n, seed)

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
                             flip_seed=False, seed=0, site_by_pos=None, selco=None):
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
    # The two quantities are resolved from DIFFERENT sources on purpose:
    #
    #   site id  is a property of THIS PANEL -- which row of this tree sequence the
    #            effect lands on -- so it is looked up in `vars`, keyed by the
    #            recipient. select_gtex_topup guarantees the recipient is there.
    #   selco    is a property of the MUTATION, identical in every panel that
    #            carries it, so it is passed in from the donor's own row. It used
    #            to be looked up in `vars` too, which is wrong whenever `positions`
    #            and `vars` come from different panels -- see combine_phenos_to_df.
    #
    # `site_by_pos` is an optional {position: id} map; without it each call did a
    # linear scan of the whole vars frame.
    if recipient_position is None:
        recipient_position = position
    if site_by_pos is not None:
        site_id = int(site_by_pos[recipient_position])
    else:
        site_id = int(vars[vars['position'] == recipient_position]['id'].iloc[0])
    if selco is None:
        selco = vars[vars['position'] == position]['selco'].iloc[0]
    selco = float(selco)
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

    # Look-up map built once instead of a linear scan of `vars` per trait.
    site_by_pos = dict(zip(vars['position'], vars['id']))

    # `positions` and `vars` are NOT always frames from the same panel, and the
    # donor's selection coefficient is read off the row rather than out of `vars`
    # for exactly that reason. Under causal_sampling=power the shared central rows
    # are DONOR rows from the GWAS frame -- select_gtex_topup returns `chosen`
    # verbatim -- and with neutral_trait_vars the donor need not segregate in the
    # GTEx panel at all; three quarters of them do not, so the old
    # selco_by_pos[donor] lookup raised KeyError on essentially every B run under
    # power sampling. selco is mutation metadata and identical in both panels, so
    # the row carries the same value the lookup found wherever it succeeded.
    #
    # Only the recipient has to be a row of THIS tree sequence, since its site id
    # is where tstrait puts the effect.
    trait_pos = [redistribution.get(int(p), p) if redistribution else p
                 for p in positions['position']]
    absent = sorted({p for p in trait_pos if p not in site_by_pos})
    if absent:
        raise SystemExit(
            f'{len(absent)} of {len(trait_pos)} trait positions are not in this panel '
            f'(first few: {absent[:5]}). select_gtex_topup is supposed to guarantee '
            'every trait position segregates here; a caller has changed.')

    # Accumulate and concat once. Concatenating inside the loop rebuilt both
    # frames on every iteration, which is quadratic in the number of traits.
    pheno_cols = []
    trait_rows = []
    for pos, selco, recipient in zip(positions['position'], positions['selco'], trait_pos):
        tr = generate_phenos_from_pos(pos, ts, vars, scaling, recipient_position=recipient,
                                      flip_seed=flip_seed, seed=seed,
                                      site_by_pos=site_by_pos, selco=selco)
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

# get_commands() removed. It wrote a copy-paste shell crib to
# gwas_scaling_*_gtex_*_maf_*.txt for a decommissioned qsub/muhee workflow:
# hardcoded `plink2 --maf 0.01` (which no longer even tracks the config),
# `--pca` steps that went away with create_pca(), and scp lines to a host
# nobody uses. It documented what someone once typed, not what ran. The
# stage2_params.txt written at the end of __main__ replaces it with an actual
# record of the parameters used and the pool sizes they produced.

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
    parser.add_argument('--min_maf', type=float, help='The minor allele frequency variants must reach to be tested in the GWAS (the pipeline passes this same value to plink2 --maf in stage 5). Does NOT gate causative-variant eligibility -- see --causal_min_maf.')
    parser.add_argument('--causal_min_maf', type=float, required=False, default=0.0,
                        help='MAF floor a variant must clear to be ELIGIBLE AS CAUSATIVE: the '
                             'central causal pool (in both the GWAS and the GTEx sample, and '
                             'under BOTH --causal_sampling schemes), the central neutral '
                             'recipient pool under --neutral_trait_vars, and the flanking '
                             'GTEx-only loci. Also names the output directory, since it '
                             'is the only MAF floor that changes what this script produces. '
                             'Independent of --min_maf (which gates the GWAS) and of the '
                             "pipeline's fm_min_maf (which gates the SBAMS handed to DAP-G): set "
                             'it below those two and the true causal variant falls out of the '
                             'tested set, reachable only through a variant that tags it. Default '
                             '0 = no floor, which is NOT the same as no filter: MAF > 0 '
                             'well-formedness applies regardless, since a site monomorphic within '
                             'a panel is not a variant there. Deliberately not defaulted to '
                             '--min_maf, so a hand invocation cannot couple the two floors.')
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
    parser.add_argument('--synthetic_dfe_effects', type=str2bool, required=False,
                        default=False,
                        help='For tree sequences with no selection in them at all (the '
                             'category-H neutral simulations). Every polymorphic variant '
                             'in the eligible window is a causal candidate, and its '
                             'effect-size parameter is DRAWN from the DFE with the neutral '
                             'class removed (helpers/synthetic_dfe.py) rather than read '
                             'out of the tree sequence. One draw per position, shared '
                             'across the GWAS and GTEx panels, seeded by --seed; beta is '
                             'then sqrt(|s|) * the panel multiplier exactly as in the '
                             'selected categories. The GTEx central set is built by '
                             'top-up in both sampling schemes. Mutually exclusive with '
                             '--neutral_trait_vars, which is a different neutral model: '
                             'there the genome is still the output of a run under '
                             'selection and only the assignment of effects is neutral.')
    parser.add_argument('--n_central_traits', type=int, required=False, default=None, help='Number of central trait loci to keep -- these are the GWAS traits AND the shared GTEx traits (same positions). Drawn from the eligible causal pool (selco != 0, MAF >= --causal_min_maf, central [5e5, L-5e5] window). Randomly subsampled (seeded by --seed) when the pool is larger; all are used when fewer exist. Default None = use all eligible (legacy behavior).')
    parser.add_argument('--causal_sampling', choices=('uniform', 'power'), required=False,
                        default='uniform',
                        help="How the central causal loci are drawn. 'uniform' (default, and "
                             'what every run through round 3 used) takes a uniform sample of the '
                             'pool that clears --causal_min_maf in both panels; the GWAS and '
                             'shared-GTEx causal sets are then identical by construction. '
                             "'power' instead weights every eligible central variant by its "
                             'probability of being detected in a GWAS of --sampling_gwas_n '
                             'individuals and draws with inclusion probability proportional to '
                             'that power. --causal_min_maf gates the central pool under both '
                             'schemes. What still differs is the GTEx side: in power mode the '
                             'GWAS causal set is not intersected with the GTEx panel, so it may '
                             'contain positions absent there -- the GTEx central set is then '
                             'topped up uniformly, and the pairing is written to '
                             '{prefix}trait_partners_*.tsv.')
    parser.add_argument('--sampling_gwas_n', type=float, required=False, default=None,
                        help='GWAS size (individuals) the detection power in --causal_sampling '
                             'power is computed for. The pipeline passes the config\'s gwas_size '
                             '(8000 for categories A and E). Defaults to the realized GWAS panel '
                             'size for a hand invocation. Larger values shift the causal set '
                             'toward rarer and smaller-effect variants.')
    parser.add_argument('--sampling_sig_p', type=float, required=False, default=5e-8,
                        help='Two-sided p-value threshold defining "detected" for the power '
                             'weight (default 5e-8, matching GWAS_SIG_P in '
                             'figures_and_tables/figure2_revision2.ipynb). Only used under '
                             '--causal_sampling power.')
    parser.add_argument('--sampling_min_power', type=float, required=False, default=0.05,
                        help='Power a pool variant must reach to count toward the sanity check '
                             'below (default 0.05). Variants under it stay in the pool and keep '
                             'their (tiny) weight; this is a diagnostic floor, not a filter.')
    parser.add_argument('--sampling_min_pool_multiple', type=float, required=False, default=2.0,
                        help='Refuse to run unless at least this multiple of --n_central_traits '
                             'pool variants reach --sampling_min_power (default 2). Guards '
                             'against a "weighted" draw that is really a handful of certainties '
                             'plus an arbitrary tail.')
    parser.add_argument('--n_flank_gtex_traits', type=int, required=False, default=50, help='Number of GTEx-ONLY trait loci drawn from the two 500 kb edge flanks ([0,5e5] U [L-5e5,L]) with the same causal eligibility (selco != 0, MAF >= --causal_min_maf in the GTEx sample). Added to the GTEx trait set on top of the shared central loci, so each GWAS locus can be compared to GTEx loci >=500 kb away. Randomly subsampled (seeded); uses all if fewer exist. Default 50; set 0 to disable.')
    return parser.parse_args()

if __name__ == '__main__':
    print('Getting arguments')
    args = get_arguments()
    gwas_scaling = args.gwas_scaling
    gtex_scaling = args.gtex_scaling
    r2_value = args.r2_value
    min_maf = args.min_maf
    # The floor for causative-variant eligibility, and the one that names the
    # output dir. min_maf is only passed through to the GWAS downstream; it has
    # no effect on anything this script decides.
    causal_min_maf = args.causal_min_maf
    gtex_size = args.gtex_size

    # Two neutral models, and they are not composable. --neutral_trait_vars keeps
    # a genome shaped by selection and moves each donor's coefficient onto a
    # neutral recipient; --synthetic_dfe_effects starts from a genome with no
    # selection in it and draws the coefficient outright. Running both would
    # redistribute drawn effects from donors that were themselves arbitrary.
    synthetic_effects = args.synthetic_dfe_effects
    if synthetic_effects and args.neutral_trait_vars:
        raise SystemExit(
            '--synthetic_dfe_effects and --neutral_trait_vars are mutually exclusive: '
            'they are two different neutral models. Pick the one the category calls for '
            '(H uses --synthetic_dfe_effects; B and D use --neutral_trait_vars).'
        )

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

    # Select the relevant alleles.
    #
    # Eligibility is MAF >= causal_min_maf in the GWAS sample, a position inside
    # the central window, and (unless --synthetic_dfe_effects) selco != 0. Under
    # the uniform scheme the pool is ALSO intersected with the GTEx panel and
    # re-tested against causal_min_maf there, which is what makes the GWAS and
    # shared-GTEx causal sets identical by construction. Note the two MAF tests
    # are applied in sequence, so the counts are logged after BOTH -- the pool
    # size is exactly the quantity that moves when causal_min_maf moves, so it
    # must not be reported one filter early.
    print('Selecting relevant alleles')
    pool_counts = {}
    power_sampling = args.causal_sampling == 'power'
    if power_sampling and args.n_central_traits is None:
        raise SystemExit(
            '--causal_sampling power needs an explicit --n_central_traits: at the '
            'default --causal_min_maf 0 the pool it draws from is every polymorphic '
            'central variant, so "use all eligible" would make every one of them '
            'causative.'
        )

    # Whether the central GTEx trait set is built by select_gtex_topup (the GWAS
    # loci that GTEx happens to carry, plus a uniform fill) rather than being the
    # GWAS set verbatim. Power sampling has always worked that way, because it does
    # not intersect its pool with the GTEx panel. The synthetic-DFE categories do
    # too, in BOTH schemes, because that is the model as specified: draw 50 in the
    # GWAS panel, whichever of them segregate in GTEx become causal there, top the
    # rest up.
    #
    # Consequence, and it is not a bug: a GWAS locus the GTEx panel lacks has no
    # partner and cannot colocalize, while under uniform+non-synthetic every locus
    # has one by construction. Both scorers currently count a partnerless locus as
    # a colocalization failure, so read {panel}_trait_partners_*.tsv before
    # comparing a rate from here against one from a uniform non-synthetic arm.
    topup_gtex = power_sampling or synthetic_effects

    # {position: drawn selection coefficient}, per species, under
    # --synthetic_dfe_effects. Built over every pool that can yield a trait locus
    # -- central GWAS, central GTEx, GTEx flanks -- so a locus causal in more than
    # one of them carries one coefficient, and so the assignment does not depend on
    # which loci the draw happens to pick. Distinct seed offsets keep the two
    # species independent; 6e7/7e7 stay clear of the 1e7-5e7 offsets the trait
    # draws already use.
    hsynth = csynth = (None, None)

    def apply_synth(frame, synth):
        """Put the drawn coefficients on a GTEx-side trait frame, or pass it through.

        The GWAS-side pool gets this inside central_pool, before the power weight is
        computed. The GTEx top-up and flank frames are drawn from the GTEx panel and
        never pass through there, so they get it here -- combine_phenos_to_df reads
        `selco` off each row to set beta, and without this they would carry the
        tree sequence's 0 and produce traits with no genetic component at all."""
        if not synthetic_effects:
            return frame
        return synthetic_dfe.apply_to(frame, *synth)

    def synthetic_maps(gwas_vars, gtex_vars, seed):
        pools = [
            causal_eligible(gwas_vars, trait_pos_lo, trait_pos_hi, causal_min_maf,
                            synthetic=True),
            causal_eligible(gtex_vars, trait_pos_lo, trait_pos_hi, causal_min_maf,
                            synthetic=True),
            flank_eligible(gtex_vars, causal_min_maf, trait_pos_lo, trait_pos_hi,
                           synthetic=True),
        ]
        smap, tmap = synthetic_dfe.assign_by_position(pools, seed)
        print(f'Synthetic DFE: drew {len(smap)} coefficients (one per candidate '
              f'position across the central GWAS, central GTEx and flank pools); '
              f'median |s| {np.median(np.abs(list(smap.values()))):.4e}')
        return smap, tmap

    def central_pool(gwas_vars, gtex_vars, gtex_maf, species_label, synth=(None, None)):
        """The central causal pool for one species, with its log line and counts."""
        pool = causal_eligible(gwas_vars, trait_pos_lo, trait_pos_hi,
                               causal_min_maf, synthetic=synthetic_effects)
        if synthetic_effects:
            # Before anything reads selco -- select_central_power weights on
            # sqrt(|selco|) * scaling, so the draw has to be in place first. This
            # is the "draw s for the whole pool, then compute power" order.
            pool = synthetic_dfe.apply_to(pool, *synth)
        pool_counts['causal_eligible_gwas_sample'] = int(pool.shape[0])
        selco_clause = '' if synthetic_effects else 'selco != 0, '
        if topup_gtex:
            pool_counts['causal_eligible'] = int(pool.shape[0])
            print(f'Number of causative {species_label} variants ({selco_clause}central '
                  f'window, MAF >= {causal_min_maf} in the GWAS sample): '
                  f'{pool_counts["causal_eligible"]} (not intersected with the GTEx '
                  'panel -- the GTEx central set is topped up)')
        else:
            pool = pool.loc[pool['position'].isin(gtex_vars['position'])]
            # `> 0` as well as the floor: at causal_min_maf = 0 the floor alone
            # would keep a site the GTEx panel carries but is monomorphic for, and
            # an effect there is not an eQTL. Implied at any non-zero floor, so
            # this is inert for every run that predates the default change.
            pool = pool[pool.apply(
                lambda row: gtex_maf[row.position] > 0
                and gtex_maf[row.position] >= causal_min_maf, axis=1)]
            pool_counts['causal_eligible'] = int(pool.shape[0])
            print(f'Number of causative {species_label} variants ({selco_clause}MAF >= '
                  f'{causal_min_maf} in the GWAS sample): '
                  f'{pool_counts["causal_eligible_gwas_sample"]}; also >= '
                  f'{causal_min_maf} in the GTEx sample: {pool_counts["causal_eligible"]}')
        return pool
    # Power is computed for a GWAS of this size, which is the config's gwas_size when
    # the pipeline invokes us. A hand invocation with no --sampling_gwas_n falls back to
    # the panel that was actually built, so the weight describes this run's own GWAS.
    realized_gwas_n = (hgwas.num_individuals if run_hums else
                       cgwas.num_individuals if run_cows else 0)
    sampling_n = (float(args.sampling_gwas_n) if args.sampling_gwas_n is not None
                  else float(realized_gwas_n))
    if power_sampling:
        print(f'Causal sampling: power-weighted (pi-PS), sampling_gwas_n={sampling_n:g} '
              f'(realized GWAS panel {realized_gwas_n}), sig_p={args.sampling_sig_p:g}')
    if run_hums:
        # Only the GTEx-intersecting path needs this, and building it costs a full
        # iterrows() pass over ~200k variants.
        hgtex_maf_dict = ({} if topup_gtex else
                          {var.position: var.maf for _, var in hgtex_vars.iterrows()})
        if synthetic_effects:
            hsynth = synthetic_maps(hgwas_vars, hgtex_vars, args.seed + 6*10**7)
        hcausative_maf01 = central_pool(hgwas_vars, hgtex_vars, hgtex_maf_dict,
                                        'human', hsynth)
    if run_cows:
        cgtex_maf_dict = ({} if topup_gtex else
                          {var.position: var.maf for _, var in cgtex_vars.iterrows()})
        if synthetic_effects:
            csynth = synthetic_maps(cgwas_vars, cgtex_vars, args.seed + 7*10**7)
        ccausative_maf01 = central_pool(cgwas_vars, cgtex_vars, cgtex_maf_dict,
                                        'cattle', csynth)

    # If --neutral_trait_vars is set, redistribute each non-zero effect from
    # its selco != 0 donor to a random selco == 0 recipient. Recipients are
    # filtered identically to donors (MAF in GWAS + position window + GTEx
    # cross-check where the donors get one) -- including the causal_min_maf floor,
    # since in these categories the recipients ARE the causative variants. With
    # fewer recipients than donors, donors are randomly truncated (seeded by
    # --seed) so the pairing stays one-to-one.
    #
    # Under --causal_sampling power the recipient pool loses the GTEx cross-check,
    # matching the donor pool, and the power weight is computed from the
    # RECIPIENT's frequency with the donor's beta -- the recipient is where the
    # effect lands, so it is the only frequency any downstream test sees. It keeps
    # the causal_min_maf floor: that used to be dropped too, which meant B's
    # recipients under power came from a pool the config's floor did not describe.
    hredist = credist = None
    hpower_maf = cpower_maf = None
    recipient_filter = (f'MAF >= {causal_min_maf} in the GWAS panel' if topup_gtex
                        else f'MAF >= {causal_min_maf} in both panels')

    def neutral_recipients(gwas_vars, gtex_vars, gtex_maf):
        """The selco == 0 recipient pool -- the donor predicate with selco inverted."""
        # maf > 0 for the same reason as causal_eligible: a recipient monomorphic
        # in this panel would carry the effect nowhere.
        pool = gwas_vars.loc[(gwas_vars['selco'] == 0)
                             & (gwas_vars['maf'] > 0)
                             & (gwas_vars['maf'] >= causal_min_maf)
                             & (gwas_vars['position'] > trait_pos_lo)
                             & (gwas_vars['position'] < trait_pos_hi)]
        if not topup_gtex:
            pool = pool.loc[pool['position'].isin(gtex_vars['position'])]
            pool = pool[pool.apply(
                lambda row: gtex_maf[row.position] > 0
                and gtex_maf[row.position] >= causal_min_maf, axis=1)]
        return pool

    if args.neutral_trait_vars:
        if run_hums:
            hneutral_maf01 = neutral_recipients(hgwas_vars, hgtex_vars, hgtex_maf_dict)
            hkept, hredist = build_redistribution_map(hcausative_maf01['position'], hneutral_maf01['position'], seed=args.seed)
            print(f'Human donors: {len(hcausative_maf01)}, neutral recipients '
                  f'({recipient_filter}): {len(hneutral_maf01)}, paired: {len(hkept)}')
            pool_counts['neutral_recipients'] = int(len(hneutral_maf01))
            pool_counts['donors_paired'] = int(len(hkept))
            hcausative_maf01 = hcausative_maf01[hcausative_maf01['position'].astype(int).isin(hkept)]
            if power_sampling:
                _hmaf = dict(zip(hgwas_vars['position'], hgwas_vars['maf']))
                hpower_maf = {p: _hmaf[hredist[int(p)]] for p in hcausative_maf01['position']}
        if run_cows:
            cneutral_maf01 = neutral_recipients(cgwas_vars, cgtex_vars, cgtex_maf_dict)
            ckept, credist = build_redistribution_map(ccausative_maf01['position'], cneutral_maf01['position'], seed=args.seed)
            print(f'Cattle donors: {len(ccausative_maf01)}, neutral recipients '
                  f'({recipient_filter}): {len(cneutral_maf01)}, paired: {len(ckept)}')
            pool_counts['neutral_recipients'] = int(len(cneutral_maf01))
            pool_counts['donors_paired'] = int(len(ckept))
            ccausative_maf01 = ccausative_maf01[ccausative_maf01['position'].astype(int).isin(ckept)]
            if power_sampling:
                _cmaf = dict(zip(cgwas_vars['position'], cgwas_vars['maf']))
                cpower_maf = {p: _cmaf[credist[int(p)]] for p in ccausative_maf01['position']}

    # Down-select trait loci.
    #  - GWAS traits = central causal loci, cut to --n_central_traits (uniformly
    #    under the default sampling scheme, where they are ALSO the shared GTEx
    #    loci because the eligible pool was intersected with the GTEx panel).
    #  - GTEx traits = the shared central loci PLUS --n_flank_gtex_traits
    #    GTEx-only loci from the 500 kb edge flanks.
    # Distinct seed offsets keep the central and flank draws (and the two
    # species) independent but deterministic given --seed (= stage2_seed).
    #
    # Under --causal_sampling power the central GWAS draw is pi-PS by detection
    # power and is NOT intersected with the GTEx panel, so the shared GTEx set is
    # only whichever of the drawn loci the GTEx panel carries; the remaining
    # central GTEx slots are filled uniformly (offsets 4e7 / 5e7). That is the one
    # place the two panels' causal sets can differ -- see the partner table written
    # below.
    print('Selecting trait loci')
    power_diag = {}
    partners = {}
    if run_hums:
        if power_sampling:
            hcausative_maf01, power_diag['hgwas'] = select_central_power(
                hcausative_maf01, args.n_central_traits, gwas_scaling, sampling_n,
                args.sampling_sig_p, args.sampling_min_power, args.sampling_min_pool_multiple,
                args.seed, maf_by_position=hpower_maf, label='human')
        else:
            hcausative_maf01 = subsample_traits(hcausative_maf01, args.n_central_traits, args.seed)
        if topup_gtex:
            hshared, htopup = select_gtex_topup(hgtex_vars, hcausative_maf01, trait_pos_lo,
                                                trait_pos_hi, args.n_central_traits,
                                                args.seed + 4*10**7, redistribution=hredist,
                                                causal_min_maf=causal_min_maf,
                                                synthetic=synthetic_effects)
            htopup = apply_synth(htopup, hsynth)
            hgtex_central = [hshared, htopup]
        else:
            hshared, htopup = hcausative_maf01, hcausative_maf01.iloc[:0]
            # Concatenated as a single frame, exactly as before this option existed --
            # an empty frame in the list would be a needless chance to shift a dtype.
            hgtex_central = [hcausative_maf01]
        hflank_gtex = apply_synth(
            select_flank_gtex(hgtex_vars, causal_min_maf, trait_pos_lo, trait_pos_hi,
                              args.n_flank_gtex_traits, args.seed + 10**7,
                              synthetic=synthetic_effects), hsynth)
        hgtex_trait_pos = pd.concat(hgtex_central + [hflank_gtex])
        partners['hgwas'] = trait_partner_table(hcausative_maf01, hshared, hredist)
        print(f'Human trait loci -- GWAS/central: {hcausative_maf01.shape[0]}, '
              f'GTEx shared: {hshared.shape[0]}, GTEx central top-up: {htopup.shape[0]}, '
              f'GTEx flank: {hflank_gtex.shape[0]}, GTEx total: {hgtex_trait_pos.shape[0]}')
        pool_counts.update(central_kept=int(hcausative_maf01.shape[0]),
                           gtex_shared_kept=int(hshared.shape[0]),
                           gtex_topup_kept=int(htopup.shape[0]),
                           flank_kept=int(hflank_gtex.shape[0]),
                           gtex_traits_total=int(hgtex_trait_pos.shape[0]))
    if run_cows:
        if power_sampling:
            ccausative_maf01, power_diag['cgwas'] = select_central_power(
                ccausative_maf01, args.n_central_traits, gwas_scaling, sampling_n,
                args.sampling_sig_p, args.sampling_min_power, args.sampling_min_pool_multiple,
                args.seed + 2*10**7, maf_by_position=cpower_maf, label='cattle')
        else:
            ccausative_maf01 = subsample_traits(ccausative_maf01, args.n_central_traits, args.seed + 2*10**7)
        if topup_gtex:
            cshared, ctopup = select_gtex_topup(cgtex_vars, ccausative_maf01, trait_pos_lo,
                                                trait_pos_hi, args.n_central_traits,
                                                args.seed + 5*10**7, redistribution=credist,
                                                causal_min_maf=causal_min_maf,
                                                synthetic=synthetic_effects)
            ctopup = apply_synth(ctopup, csynth)
            cgtex_central = [cshared, ctopup]
        else:
            cshared, ctopup = ccausative_maf01, ccausative_maf01.iloc[:0]
            cgtex_central = [ccausative_maf01]   # see the human branch
        cflank_gtex = apply_synth(
            select_flank_gtex(cgtex_vars, causal_min_maf, trait_pos_lo, trait_pos_hi,
                              args.n_flank_gtex_traits, args.seed + 3*10**7,
                              synthetic=synthetic_effects), csynth)
        cgtex_trait_pos = pd.concat(cgtex_central + [cflank_gtex])
        partners['cgwas'] = trait_partner_table(ccausative_maf01, cshared, credist)
        print(f'Cattle trait loci -- GWAS/central: {ccausative_maf01.shape[0]}, '
              f'GTEx shared: {cshared.shape[0]}, GTEx central top-up: {ctopup.shape[0]}, '
              f'GTEx flank: {cflank_gtex.shape[0]}, GTEx total: {cgtex_trait_pos.shape[0]}')
        pool_counts.update(central_kept=int(ccausative_maf01.shape[0]),
                           gtex_shared_kept=int(cshared.shape[0]),
                           gtex_topup_kept=int(ctopup.shape[0]),
                           flank_kept=int(cflank_gtex.shape[0]),
                           gtex_traits_total=int(cgtex_trait_pos.shape[0]))

    # How many GWAS causal loci actually have a same-named GTEx partner. Under
    # uniform sampling this is all of them by construction; under power sampling it
    # is not, and both scorers (helpers/summarize_coloc.py and the figure-2 notebook)
    # currently define a true positive as trait-name equality and keep every GWAS
    # trait in the denominator -- so a partnerless locus reads as a colocalization
    # failure, and as a false positive if it hits a neighbour. Fixing that is
    # deliberately out of scope here; this number is how we know how much it matters.
    for tbl in partners.values():
        n_shared = int(tbl['shared'].sum())
        print(f'GWAS causal loci with a GTEx partner: {n_shared} of {tbl.shape[0]}'
              + ('' if n_shared == tbl.shape[0] else
                 ' -- the rest cannot colocalize and will score as failures'))
    if power_sampling:
        for name, diag in power_diag.items():
            sel = diag.loc[diag['selected']]
            pool_counts[f'{name}_power_above_floor'] = int(
                (diag['power'] >= args.sampling_min_power).sum())
            pool_counts[f'{name}_power_sum'] = float(diag['power'].sum())
            pool_counts[f'{name}_power_median_pool'] = float(diag['power'].median())
            pool_counts[f'{name}_power_median_selected'] = float(sel['power'].median())
            print(f'{name} power: pool median {diag["power"].median():.3g}, selected median '
                  f'{sel["power"].median():.3g}, {pool_counts[f"{name}_power_above_floor"]} of '
                  f'{diag.shape[0]} at or above {args.sampling_min_power:g}')

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
    # The maf_ component is causal_min_maf, NOT min_maf: it is the only one of the
    # three MAF floors that changes what lands in this directory.
    #
    # The value is rendered through paths.causal_maf_token rather than
    # interpolated directly, because argparse hands us a FLOAT and the config
    # hands helpers/paths.py:stage2_inner() whatever YAML parsed. At 0.01 and
    # 0.001 those format the same and the two sides agreed by accident; at the
    # new default of 0 they do not -- int 0 renders "0" and float 0.0 renders
    # "0.0", so Snakemake would declare maf_0/ as the rule's output while this
    # script wrote maf_0.0/, and the rule would fail with its outputs missing
    # after stage 2 had already done all the work. One function, both callers.
    maf_token = paths.causal_maf_token(causal_min_maf)
    sim_dir = f'{args.out_dir}/gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{maf_token}'
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    # `_causal_power` exists only under power sampling -- there is no power column
    # to report otherwise. It covers the WHOLE pool, not just the drawn loci,
    # which is what makes the draw checkable after the fact: the selected set's
    # power distribution means nothing without the pool's.
    suffix = f'gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{maf_token}.tsv'
    if power_sampling:
        for name, diag in power_diag.items():
            diag.to_csv(f'{sim_dir}/{name}_causal_power_{suffix}', sep='\t', index=False)

    # The partner table is gated on TOP-UP, not on power sampling. It used to be
    # gated on the latter because the two coincided; they stopped coinciding when
    # the drawn-DFE arms (H, I) began topping up under uniform sampling too. Those
    # are precisely the runs where a GWAS locus can have no GTEx trait at its causal
    # variant, and where both scorers would otherwise read "no eQTL to find" as
    # "colocalization failed" -- opposite results. Where the GTEx set is intersected
    # instead, every locus has a partner by construction and nothing is written, so
    # a plain uniform run's output directory is unchanged.
    if topup_gtex:
        for name, tbl in partners.items():
            tbl.to_csv(f'{sim_dir}/{name}_trait_partners_{suffix}', sep='\t', index=False)

    # The drawn coefficients, for the WHOLE candidate pool rather than the drawn
    # loci -- same reasoning as the power sidecar above: the selected set's
    # distribution says nothing without the pool's. The `*_vars_*.tsv` files
    # deliberately keep selco = 0, because these variants really do have no
    # selection coefficient; this is where the assigned effect-size parameter is
    # recorded, and it is what plot_selco_af() would have to read for category H.
    if synthetic_effects:
        for name, (smap, tmap) in (('hgwas', hsynth), ('cgwas', csynth)):
            if not smap:
                continue
            pd.DataFrame({
                'position': list(smap.keys()),
                'selco': list(smap.values()),
                'mutation_type': [tmap[p] for p in smap],
            }).to_csv(f'{sim_dir}/{name}_synthetic_dfe_{suffix}', sep='\t', index=False)

    if run_cows:
        attach_beta(cgwas_vars, cgwas_key_maf01).to_csv(f'{sim_dir}/cgwas_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{maf_token}.tsv', sep='\t', index=False)
        attach_beta(cgtex_vars, cgtex_key_maf01).to_csv(f'{sim_dir}/cgtex_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{maf_token}.tsv', sep='\t', index=False)
        cgwas_traits_maf01.to_csv(f'{sim_dir}/cgwas_traits_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{maf_token}.tsv', sep='\t', index=False)
        cgtex_traits_maf01.to_csv(f'{sim_dir}/cgtex_traits_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{maf_token}.tsv', sep='\t', index=False)

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
            attach_beta(sub_vars, cgtex_key_maf01).to_csv(f'{sim_dir}/{name}_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{maf_token}.tsv', sep='\t', index=False)
            sub_traits.to_csv(f'{sim_dir}/{name}_traits_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{maf_token}.tsv', sep='\t', index=False)
            with open(f'{sim_dir}/{name}.vcf', 'w') as f:
                sub_ts.write_vcf(f)
            write_ts_as_sbams(sub_ts, f'{sim_dir}/{name}_scaling_{gtex_scaling}_geno.sbams')
            write_traits_as_sbams(sub_traits, f'{sim_dir}/{name}_scaling_{gtex_scaling}_pheno.sbams')
            sub_traits.to_csv(f'{sim_dir}/{name}_traits.scaling_{gtex_scaling}.pheno', sep='\t', index=False)

    if run_hums:
        attach_beta(hgwas_vars, hgwas_key_maf01).to_csv(f'{sim_dir}/hgwas_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{maf_token}.tsv', sep='\t', index=False)
        attach_beta(hgtex_vars, hgtex_key_maf01).to_csv(f'{sim_dir}/hgtex_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{maf_token}.tsv', sep='\t', index=False)
        hgwas_traits_maf01.to_csv(f'{sim_dir}/hgwas_traits_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{maf_token}.tsv', sep='\t', index=False)
        hgtex_traits_maf01.to_csv(f'{sim_dir}/hgtex_traits_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{maf_token}.tsv', sep='\t', index=False)

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
            attach_beta(sub_vars, hgtex_key_maf01).to_csv(f'{sim_dir}/{name}_vars_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{maf_token}.tsv', sep='\t', index=False)
            sub_traits.to_csv(f'{sim_dir}/{name}_traits_gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{maf_token}.tsv', sep='\t', index=False)
            with gzip.open(f'{sim_dir}/{name}.vcf.gz', 'wt') as f:
                sub_ts.write_vcf(f)
            write_ts_as_sbams(sub_ts, f'{sim_dir}/{name}_scaling_{gtex_scaling}_geno.sbams')
            write_traits_as_sbams(sub_traits, f'{sim_dir}/{name}_scaling_{gtex_scaling}_pheno.sbams')
            sub_traits.to_csv(f'{sim_dir}/{name}_traits.scaling_{gtex_scaling}.pheno', sep='\t', index=False)

    # Record what this stage-2 run actually used, and what it produced.
    #
    # Replaces the old gwas_scaling_*_gtex_*_maf_*.txt, which was a copy-paste crib
    # for a decommissioned qsub/muhee workflow -- hardcoded `--maf 0.01`, plink
    # --pca steps the pipeline no longer runs -- rather than a record of anything.
    #
    # This file travels with the directory when it is symlink-adopted into another
    # output root, which is what lets the Snakefile verify that adopted phenotypes
    # were built under the causal MAF floor the adopting run actually wants.
    ts_file = args.human_ts_file if run_hums else args.cattle_ts_file
    stage1 = {'file': os.path.basename(ts_file), 'path': ts_file}
    try:
        st = os.stat(ts_file)
        stage1.update(size=st.st_size, mtime=int(st.st_mtime))
    except OSError:
        pass
    # Keyed by CONFIG name rather than CLI name, so stage2_uid here and in the
    # run-level record are computed over the same dict and can be compared.
    params_record.write_stage2_params(
        f'{sim_dir}/stage2_params.txt',
        values={
            'causal_min_maf': causal_min_maf,
            'min_maf': min_maf,
            'gwas_scaling': gwas_scaling,
            'gtex_scaling': gtex_scaling,
            'gtex_size': args.gtex_size,
            'n_samples': args.n_samples,
            'n_central_traits': args.n_central_traits,
            'n_flank_gtex_traits': args.n_flank_gtex_traits,
            'neutral_trait_vars': args.neutral_trait_vars,
            'synthetic_dfe_effects': args.synthetic_dfe_effects,
            'already_includes_neutral': args.already_includes_neutral,
            'L': args.length,
            'Q_scaling': q_scaling,
            'handoff_ticks': args.handoff_ticks,
            'deep_Q_scaling': args.deep_Q_scaling,
            'stage2_seed': args.seed,
            'causal_sampling': args.causal_sampling,
            # The RESOLVED size, not args.sampling_gwas_n, so a hand invocation that
            # let it default records the number the weights were actually built from.
            'sampling_gwas_n': sampling_n,
            'sampling_sig_p': args.sampling_sig_p,
            'sampling_min_power': args.sampling_min_power,
            'sampling_min_pool_multiple': args.sampling_min_pool_multiple,
            'species': 'human' if run_hums else 'cattle',
        },
        pools=pool_counts,
        stage1=stage1,
    )
    print(f'Wrote {sim_dir}/stage2_params.txt')



# Included in final analysis

# Not included in final analysis

# Too soon to tell
# python ../create_gwas_files_and_phenotypes.py --gwas_scaling 35 --gtex_scaling 35 --r2_value 0.2 --min_maf 0.01 --cattle_ts_file '/Users/noah/tmp/selsims/revision_farm_selection_mult_100_gen_23_muts_26_bottlenecked_sd24.full.ts' --cattle_m4_file '/Users/noah/tmp/selsims/revision_farm_selection_mult_100_gen_23_muts_26_bottlenecked_sd24.m4_marks.tsv' --out_dir '/Users/noah/tmp/selsims/outdir/'
# python ../create_gwas_files_and_phenotypes.py --gtex_size 50000 --gwas_scaling 35 --gtex_scaling 35 --r2_value 0.2 --min_maf 0.01 --human_ts_file '/Users/noah/tmp/selsims/largesim/hts_19930224_large.ts' --out_dir '/Users/noah/tmp/selsims/largesim/outputs/' --already_includes_neutral True
