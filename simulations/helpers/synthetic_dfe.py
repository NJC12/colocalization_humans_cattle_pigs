"""The distribution of fitness effects, with the neutral class removed.

Categories A and E-G get their effect sizes from selection coefficients that a
forward simulation actually assigned: SLiM draws each new mutation from a
three-component DFE, and stage 2 reads the coefficient back out of the tree
sequence. The neutral simulations (category H) have no forward phase and no
selection coefficients at all, so the coefficient has to be drawn here instead.

WHAT "REMOVE THE NEUTRAL CLASS AND RENORMALIZE" COMES OUT AS

Every arm mutates at 1.4e-8 per bp per generation, split into a DFE-applied
5.6e-9 and a neutral 8.4e-9 overlaid afterwards by msprime. So the full
distribution over new mutations is 60% neutral and 40% spread across the three
DFE components, and dropping the neutral class renormalizes the remainder back
to exactly the proportions SLiM was configured with. That is why PROPORTIONS
below is the g1 vector verbatim rather than something recomputed: the neutral
class is not a DFE component to begin with, it is a separate mutation rate.

The parameters are the same ones in four other places, and the five must not
drift:

    human_simulation_o2.py:70-91      stdpopsim MutationType / DFE
    farm_burn_in_e2.slim:35-38        initializeMutationType / g1
    farm_selection.slim:40-43         same
    farm_selection_from_ep8.slim:182-186   same, plus the m4 selection target

Coefficients here are UNSCALED (real, Q = 1). The SLiM definitions carry a
`Q *` factor and `get_vars_df` divides it back out; there is no forward phase to
rescale for, so there is nothing to undo.
"""

import numpy as np

# m1 exponential deleterious, m2 gamma deleterious, m3 exponential positive.
PROPORTIONS = (0.975, 0.024, 0.001)

M1_MEAN = 4e-4            # exponential scale
M2_MEAN = 0.03            # gamma mean; SLiM's "g" is parameterised (mean, shape)
M2_SHAPE = 0.206
M3_MEAN = 0.005           # exponential scale

# SLiM mutation-type ids, so a drawn coefficient can be reported alongside the
# class it came from. Matches the `type` column get_vars_df reads out of the
# mutation metadata (0 = neutral, 1/2/3 = the DFE components).
MUTATION_TYPES = (1, 2, 3)


def draw(n, rng):
    """`n` independent draws from the truncated DFE.

    Returns ``(selco, mutation_type)``, both length-`n` arrays. `selco` is signed
    the way SLiM signs it: negative for the two deleterious components, positive
    for m3.

    The component is chosen per element rather than by a multinomial count plus a
    permutation (which is how ``farm_create_orig_pop_e2.py:selection_coeff_bulk``
    does it). Same distribution; this form keeps element i's draw a function of
    the stream position alone, so growing the candidate pool does not reshuffle
    the coefficients assigned to the variants that were already in it.
    """
    n = int(n)
    if n == 0:
        return np.empty(0, dtype=float), np.empty(0, dtype=int)

    comp = rng.choice(len(PROPORTIONS), size=n, p=np.asarray(PROPORTIONS, dtype=float))

    # Draw every component for every element and select, rather than drawing
    # per-component counts. Costs three variates per element and buys the same
    # position-stability the per-element `choice` above is for.
    m1 = -rng.exponential(scale=M1_MEAN, size=n)
    m2 = -rng.gamma(shape=M2_SHAPE, scale=M2_MEAN / M2_SHAPE, size=n)
    m3 = rng.exponential(scale=M3_MEAN, size=n)

    selco = np.select([comp == 0, comp == 1, comp == 2], [m1, m2, m3])
    mtype = np.asarray(MUTATION_TYPES, dtype=int)[comp]
    return selco, mtype


def assign_by_position(frames, seed):
    """One draw per POSITION, shared across every panel that carries it.

    `frames` is any iterable of ``*_vars_*``-shaped DataFrames. The union of
    their positions is sorted and drawn for in that order, so the map is a
    function of (seed, position set) alone -- not of the order the frames were
    passed, nor of which panel a position happened to appear in first.

    Sharing is the point. In categories A/E-G a selection coefficient is mutation
    metadata: the same variant carries the same value in the GWAS panel and the
    GTEx panel, and only the multiplier differs between them. A per-panel draw
    would break that, and a locus causal in both panels would get two unrelated
    effect sizes.

    Returns ``(selco_by_position, mutation_type_by_position)``.
    """
    positions = sorted({float(p) for frame in frames for p in frame['position']})
    selco, mtype = draw(len(positions), np.random.default_rng(seed))
    return (dict(zip(positions, (float(s) for s in selco))),
            dict(zip(positions, (int(t) for t in mtype))))


def apply_to(frame, selco_by_position, mutation_type_by_position=None):
    """Return a copy of `frame` with `selco` (and `type`) taken from the maps.

    Overwriting `selco` rather than adding a column is deliberate: it is what
    lets every downstream consumer run unchanged. ``select_central_power``
    computes the detection-power weight as ``sqrt(|selco|) * scaling`` and
    ``generate_phenos_from_pos`` computes beta the same way, both off the row
    they are handed.

    Raises KeyError on a position the map does not cover, which is the failure
    you want: it means a trait locus was drawn from a frame that was not passed
    to ``assign_by_position``, and the effect size would otherwise be silently
    whatever the tree sequence happened to carry.
    """
    out = frame.copy()
    out['selco'] = [selco_by_position[float(p)] for p in out['position']]
    if mutation_type_by_position is not None:
        out['type'] = [mutation_type_by_position[float(p)] for p in out['position']]
    return out
