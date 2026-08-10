"""Power-weighted (pi-PS) selection of causal variants for the stage-2 simulation.

Why this exists
---------------
Stage 2 has always drawn its 50 central causal variants UNIFORMLY from the pool of
variants with a non-zero selection coefficient inside the central window.  In the
human arm that pool is dominated by variants no GWAS of realistic size could ever
find -- 52-64% of it sits below MAF 0.01 -- so the simulated "GWAS loci" are not
the kind of loci a real study reports, and the colocalization rate is measured on a
locus set that no real analysis would ever assemble.

This module implements the alternative: weight each candidate by its PROBABILITY OF
BEING DETECTED in a GWAS of a stated size, then draw the causal set with inclusion
probability proportional to that power.

The power model
---------------
The simulated phenotype is

    y = +/- beta * g + N(0, 1),      g in {0, 1, 2},

built by `create_gwas_files_and_phenotypes.py:generate_phenos_from_pos` with
beta = sqrt(|selco|) * scaling.  It is never standardized, so with allele frequency
f the genotype contributes variance 2f(1-f)beta^2 on top of the unit environmental
variance:

    Var(y)   = 1 + 2f(1-f) beta^2
    vexp     = 2f(1-f) beta^2
    E[chi^2] = ncp = N * vexp / (1 + vexp)

and the Wald statistic is chi2_1(ncp).  Power at a two-sided threshold `sig_p` is
the non-central tail P(chi2_1(ncp) > qchisq(sig_p, 1, lower=FALSE)).

This is the same identity used by `add_gwas_significance()` in
`figures_and_tables/figure2_revision2.ipynb` and stated in
`helpers/pval_rescale.py:27`; the three must not drift.  Note the factor of 2 in
vexp -- an older version of `figure2.ipynb` omitted it in one of its two power
functions, which understates the ncp by half.

Why pi-PS and not `rng.choice(p=weights)`
-----------------------------------------
numpy's weighted draw without replacement is SUCCESSIVE sampling: the first pick is
proportional to the weights, but each later pick renormalises over what is left, so
a variant's probability of appearing ANYWHERE in the final sample is not
proportional to its weight (high-weight units are over-represented relative to
their weight, and the distortion grows with the sample fraction).

`inclusion_probabilities` + `systematic_sample` give exact first-order inclusion
probabilities instead: P(variant i is in the drawn set) == pi_i, with pi_i
proportional to power_i wherever the cap does not bind.  That is what "sample
weighted by the probability of detection" means when stated precisely, and it is
what makes the resulting causal set interpretable as a sample from the
detection-weighted distribution.

The cap matters and is not a technicality.  With 50 draws from a pool where 40
variants have power ~1, a raw proportional rule would ask for pi > 1 on those 40;
the standard iterative fix pins them at 1 (they are certainties) and redistributes
the remaining budget over the rest.

Second-order (joint) inclusion probabilities of systematic sampling are not
uniform, and some pairs have zero joint probability.  That would matter for
variance estimation of a survey total; it does not matter here, where the drawn set
is the object of interest rather than an estimator.
"""

import numpy as np
from scipy import stats

__all__ = [
    "detection_power",
    "inclusion_probabilities",
    "systematic_sample",
    "significance_threshold",
]


def significance_threshold(sig_p):
    """Chi-square (1 df) critical value for a two-sided p-value threshold.

    29.716... at the genome-wide 5e-8.
    """
    sig_p = float(sig_p)
    if not (0.0 < sig_p < 1.0):
        raise ValueError(f"sig_p must lie in (0, 1), got {sig_p!r}")
    return float(stats.chi2.isf(sig_p, df=1))


def detection_power(maf, beta, n, sig_p):
    """P(the additive test rejects at `sig_p`) for the simulated phenotype model.

    `maf` and `beta` broadcast against each other; `n` is the GWAS sample size in
    INDIVIDUALS (the genotype g is 0/1/2, so no factor of two here).  See the
    module docstring for the ncp identity.

    Returns `sig_p` exactly when the ncp is zero (beta == 0): that is the size of
    the test, and `scipy.stats.ncx2` is not reliable at nc == 0.
    """
    maf = np.asarray(maf, dtype=float)
    beta = np.asarray(beta, dtype=float)
    n = float(n)

    if n <= 0:
        raise ValueError(f"n must be positive, got {n!r}")
    if np.any(~np.isfinite(maf)) or np.any((maf <= 0.0) | (maf >= 1.0)):
        raise ValueError("maf must be finite and lie strictly in (0, 1)")
    if np.any(~np.isfinite(beta)):
        raise ValueError("beta must be finite")

    thr = significance_threshold(sig_p)

    vexp = 2.0 * maf * (1.0 - maf) * beta**2
    ncp = np.asarray(n * vexp / (1.0 + vexp), dtype=float)

    power = np.full(ncp.shape, float(sig_p), dtype=float)
    nz = ncp > 0.0
    if np.any(nz):
        power[nz] = stats.ncx2.sf(thr, df=1, nc=ncp[nz])
    return power if ncp.ndim else float(power)


def inclusion_probabilities(weights, n):
    """First-order inclusion probabilities for a pi-PS sample of size `n`.

    Returns `pi` with ``sum(pi) == n`` (to float64), ``0 <= pi <= 1``, and pi
    proportional to `weights` on the units where the cap does not bind.  Units
    whose proportional share would exceed 1 are pinned at 1 and the budget is
    redistributed over the remainder, iterating until nothing is over -- the
    standard construction, and the only one that can satisfy both the total and
    the pi <= 1 constraint simultaneously.

    Zero-weight units get pi == 0 and can never be drawn; `n` must therefore not
    exceed the number of strictly positive weights.
    """
    w = np.asarray(weights, dtype=float)
    if w.ndim != 1:
        raise ValueError(f"weights must be 1-D, got shape {w.shape}")
    n = int(n)
    if n < 0:
        raise ValueError(f"n must be non-negative, got {n!r}")
    if np.any(~np.isfinite(w)) or np.any(w < 0.0):
        raise ValueError("weights must be finite and non-negative")
    if n > w.size:
        raise ValueError(f"cannot draw {n} units from a pool of {w.size}")

    pi = np.zeros_like(w)
    if n == 0:
        return pi

    positive = w > 0.0
    if int(positive.sum()) < n:
        raise ValueError(
            f"cannot draw {n} units: only {int(positive.sum())} of {w.size} have a "
            "strictly positive weight"
        )

    free = positive.copy()
    budget = n
    while True:
        share = budget * w[free] / w[free].sum()
        over = share >= 1.0
        if not over.any():
            pi[free] = share
            return pi
        pinned = np.flatnonzero(free)[over]
        pi[pinned] = 1.0
        free[pinned] = False
        budget -= pinned.size
        if budget == 0:
            return pi


def systematic_sample(pi, rng):
    """Randomized systematic pi-PS draw. Returns exactly ``round(sum(pi))`` indices.

    P(index i is returned) == pi[i] exactly.  The pool is permuted first: run on
    the natural (position-sorted) order, systematic sampling would space the picks
    evenly along the chromosome and induce spatial structure in which loci can
    co-occur.

    `rng` is a `numpy.random.Generator` supplied by the caller -- this module never
    creates one, matching the convention in `helpers/pval_rescale.py` and
    `create_gwas_files_and_phenotypes.py`.  Returned indices are sorted ascending,
    so a caller slicing a position-ordered frame keeps it in position order (the
    same contract `subsample_traits` has).
    """
    pi = np.asarray(pi, dtype=float)
    if pi.ndim != 1:
        raise ValueError(f"pi must be 1-D, got shape {pi.shape}")
    if np.any(~np.isfinite(pi)) or np.any((pi < 0.0) | (pi > 1.0 + 1e-12)):
        raise ValueError("pi must be finite and lie in [0, 1]")

    total = float(pi.sum())
    n = int(round(total))
    if abs(total - n) > 1e-8:
        raise ValueError(
            f"sum(pi) = {total!r} is not an integer sample size; build pi with "
            "inclusion_probabilities()"
        )
    if n == 0:
        return np.empty(0, dtype=int)

    order = rng.permutation(pi.size)
    cum = np.cumsum(pi[order])
    # One point per draw on a random offset grid. sum(pi) == n and every pi <= 1,
    # so each unit's interval contains at most one point: exactly n distinct units.
    points = rng.random() + np.arange(n, dtype=float)
    hit = np.searchsorted(cum, points, side="right")
    # Guard against a float32-scale overshoot at the very top of the grid.
    hit = np.minimum(hit, pi.size - 1)
    return np.sort(order[hit])
