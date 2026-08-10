"""Tests for helpers/causal_power.py.

These check the algebra and the guard rails: the ncp identity, the constraints
pi-PS is defined by (sum(pi) == n, pi <= 1, proportionality where uncapped), and
every documented error path.

One test is a CALIBRATION check rather than an algebraic one --
`test_empirical_inclusion_frequencies_match_pi` runs the sampler many times and
compares realized frequencies to their nominal pi.  It is seeded and its tolerance
is a Monte Carlo band, not an identity.
"""

import numpy as np
import pytest
from scipy import stats

from helpers.causal_power import (
    detection_power,
    inclusion_probabilities,
    significance_threshold,
    systematic_sample,
)


# --------------------------------------------------------------- the power model

def test_genome_wide_threshold_matches_qchisq():
    # qchisq(5e-8, df = 1, lower.tail = FALSE) in figure2_revision2.ipynb
    assert significance_threshold(5e-8) == pytest.approx(29.716756, rel=1e-6)


@pytest.mark.parametrize("sig_p", [5e-8, 1e-5, 0.05])
def test_power_equals_the_size_of_the_test_when_beta_is_zero(sig_p):
    assert detection_power(0.25, 0.0, 8000, sig_p) == pytest.approx(sig_p, rel=1e-12)


def test_power_matches_the_noncentral_chisquare_tail_directly():
    maf, beta, n, sig_p = 0.2, 0.15, 8000, 5e-8
    vexp = 2 * maf * (1 - maf) * beta**2
    ncp = n * vexp / (1 + vexp)
    expected = stats.ncx2.sf(stats.chi2.isf(sig_p, df=1), df=1, nc=ncp)
    assert detection_power(maf, beta, n, sig_p) == pytest.approx(expected, rel=1e-12)


def test_power_uses_the_factor_of_two_in_variance_explained():
    # The half-size ncp (the missing-2 bug in the older figure2.ipynb) gives a
    # materially different answer, so the test would catch a regression to it.
    maf, beta, n, sig_p = 0.05, 0.2, 8000, 5e-8
    vexp = 2 * maf * (1 - maf) * beta**2
    half = maf * (1 - maf) * beta**2
    thr = stats.chi2.isf(sig_p, df=1)
    got = detection_power(maf, beta, n, sig_p)
    assert got == pytest.approx(stats.ncx2.sf(thr, df=1, nc=n * vexp / (1 + vexp)), rel=1e-12)
    assert got != pytest.approx(stats.ncx2.sf(thr, df=1, nc=n * half / (1 + half)), rel=1e-3)


def test_power_is_monotone_in_effect_size_sample_size_and_maf():
    assert detection_power(0.2, 0.05, 8000, 5e-8) < detection_power(0.2, 0.10, 8000, 5e-8)
    assert detection_power(0.2, 0.10, 4000, 5e-8) < detection_power(0.2, 0.10, 8000, 5e-8)
    assert detection_power(0.02, 0.10, 8000, 5e-8) < detection_power(0.20, 0.10, 8000, 5e-8)


def test_power_is_symmetric_in_the_sign_of_beta():
    assert detection_power(0.3, -0.12, 8000, 5e-8) == detection_power(0.3, 0.12, 8000, 5e-8)


def test_power_saturates_for_a_large_effect():
    assert detection_power(0.4, 2.0, 8000, 5e-8) == pytest.approx(1.0, abs=1e-12)


def test_power_broadcasts_over_arrays():
    maf = np.array([0.01, 0.1, 0.4])
    beta = np.array([0.5, 0.2, 0.05])
    got = detection_power(maf, beta, 8000, 5e-8)
    assert got.shape == (3,)
    for i in range(3):
        assert got[i] == pytest.approx(detection_power(maf[i], beta[i], 8000, 5e-8), rel=1e-12)


@pytest.mark.parametrize("bad", [0.0, 1.0, -0.1, 1.5, np.nan])
def test_power_rejects_maf_outside_the_open_unit_interval(bad):
    with pytest.raises(ValueError, match="maf must be finite"):
        detection_power(bad, 0.1, 8000, 5e-8)


def test_power_rejects_non_positive_sample_size():
    with pytest.raises(ValueError, match="n must be positive"):
        detection_power(0.2, 0.1, 0, 5e-8)


@pytest.mark.parametrize("bad", [0.0, 1.0, -1.0, 2.0])
def test_power_rejects_a_threshold_outside_the_open_unit_interval(bad):
    with pytest.raises(ValueError, match=r"sig_p must lie in \(0, 1\)"):
        detection_power(0.2, 0.1, 8000, bad)


# ------------------------------------------------------- inclusion probabilities

def test_inclusion_probabilities_sum_to_the_sample_size():
    w = np.array([0.9, 0.5, 0.2, 0.05, 0.01, 1e-6])
    pi = inclusion_probabilities(w, 3)
    assert pi.sum() == pytest.approx(3.0, abs=1e-12)


def test_inclusion_probabilities_never_exceed_one():
    # A pool where a proportional rule alone would demand pi > 1.
    w = np.array([1.0, 1.0, 1.0, 1e-6, 1e-6])
    pi = inclusion_probabilities(w, 4)
    assert np.all(pi <= 1.0 + 1e-12)
    assert pi.sum() == pytest.approx(4.0, abs=1e-12)


def test_saturated_units_are_pinned_at_certainty():
    # w[0] would take 2 * 100/103 = 1.94 of the budget under a bare proportional
    # rule, so it is pinned at 1 and the remaining 1 is split over the rest.
    w = np.array([100.0, 1.0, 1.0, 1.0])
    pi = inclusion_probabilities(w, 2)
    assert pi[0] == pytest.approx(1.0, abs=1e-12)
    assert np.allclose(pi[1:], 1.0 / 3.0)


def test_equal_weights_give_a_uniform_sample_fraction():
    pi = inclusion_probabilities(np.ones(10), 4)
    assert np.allclose(pi, 0.4)


def test_pi_is_proportional_to_weight_where_the_cap_does_not_bind():
    w = np.array([0.4, 0.2, 0.1, 0.05])
    pi = inclusion_probabilities(w, 1)
    assert np.all(pi < 1.0)
    assert np.allclose(pi / w, pi[0] / w[0])


def test_pi_preserves_the_ordering_of_the_weights():
    rng = np.random.default_rng(20260801)
    w = np.sort(rng.random(40))[::-1] + 1e-9
    pi = inclusion_probabilities(w, 12)
    assert np.all(np.diff(pi) <= 1e-12)


def test_capping_is_idempotent():
    w = np.array([5.0, 4.0, 0.3, 0.2, 0.1])
    pi = inclusion_probabilities(w, 3)
    assert np.allclose(inclusion_probabilities(pi, 3), pi)


def test_zero_weight_units_can_never_be_drawn():
    pi = inclusion_probabilities(np.array([1.0, 0.0, 1.0, 0.0]), 2)
    assert pi[1] == 0.0 and pi[3] == 0.0


def test_drawing_the_whole_pool_makes_every_unit_certain():
    pi = inclusion_probabilities(np.array([0.9, 0.1, 0.5]), 3)
    assert np.allclose(pi, 1.0)


def test_a_zero_size_sample_gives_all_zero_pi():
    assert np.all(inclusion_probabilities(np.array([1.0, 2.0]), 0) == 0.0)


def test_inclusion_probabilities_reject_a_sample_larger_than_the_pool():
    with pytest.raises(ValueError, match="cannot draw 5 units from a pool of 3"):
        inclusion_probabilities(np.ones(3), 5)


def test_inclusion_probabilities_reject_too_few_positive_weights():
    with pytest.raises(ValueError, match="only 2 of 5 have a strictly positive weight"):
        inclusion_probabilities(np.array([1.0, 0.0, 2.0, 0.0, 0.0]), 3)


@pytest.mark.parametrize("bad", [np.array([1.0, -1.0]), np.array([1.0, np.nan]),
                                 np.array([1.0, np.inf])])
def test_inclusion_probabilities_reject_negative_or_non_finite_weights(bad):
    with pytest.raises(ValueError, match="finite and non-negative"):
        inclusion_probabilities(bad, 1)


# ------------------------------------------------------------- systematic sample

def test_systematic_sample_returns_exactly_n_distinct_sorted_indices():
    rng = np.random.default_rng(11)
    pi = inclusion_probabilities(rng.random(200) + 1e-6, 50)
    idx = systematic_sample(pi, rng)
    assert idx.size == 50
    assert np.unique(idx).size == 50
    assert np.all(np.diff(idx) > 0)
    assert idx.min() >= 0 and idx.max() < 200


def test_certain_units_are_always_drawn():
    rng = np.random.default_rng(12)
    w = np.array([10.0, 10.0, 10.0] + [0.001] * 50)
    pi = inclusion_probabilities(w, 5)
    assert np.allclose(pi[:3], 1.0)
    for _ in range(200):
        assert set([0, 1, 2]).issubset(set(systematic_sample(pi, rng).tolist()))


def test_zero_weight_units_are_never_drawn():
    rng = np.random.default_rng(13)
    w = np.array([1.0, 0.0, 1.0, 0.0, 1.0, 1.0])
    pi = inclusion_probabilities(w, 2)
    for _ in range(200):
        assert not ({1, 3} & set(systematic_sample(pi, rng).tolist()))


def test_the_draw_is_reproducible_from_the_seed():
    pi = inclusion_probabilities(np.linspace(0.01, 1.0, 60), 20)
    a = systematic_sample(pi, np.random.default_rng(2026))
    b = systematic_sample(pi, np.random.default_rng(2026))
    assert np.array_equal(a, b)


def test_empirical_inclusion_frequencies_match_pi():
    """Calibration, not algebra: P(i drawn) should equal pi_i.

    20k replicates gives SE <= sqrt(0.25/20000) = 0.0035 per unit; the 4 SE band
    below is loose enough not to flake and tight enough that the successive-draw
    alternative (`rng.choice(p=w, replace=False)`), which over-weights the top
    units by several percentage points here, would fail it.
    """
    rng = np.random.default_rng(20260801)
    w = np.array([0.9, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01])
    pi = inclusion_probabilities(w, 4)
    reps = 20000
    counts = np.zeros(w.size)
    for _ in range(reps):
        counts[systematic_sample(pi, rng)] += 1
    freq = counts / reps
    se = np.sqrt(pi * (1 - pi) / reps)
    assert np.all(np.abs(freq - pi) <= 4 * se + 1e-12)


def test_systematic_sample_rejects_pi_that_is_not_an_integer_total():
    with pytest.raises(ValueError, match="is not an integer sample size"):
        systematic_sample(np.array([0.3, 0.3, 0.3]), np.random.default_rng(0))


@pytest.mark.parametrize("bad", [np.array([1.5, 0.5]), np.array([-0.5, 1.0, 0.5]),
                                 np.array([np.nan, 1.0])])
def test_systematic_sample_rejects_pi_outside_the_unit_interval(bad):
    with pytest.raises(ValueError, match=r"pi must be finite and lie in \[0, 1\]"):
        systematic_sample(bad, np.random.default_rng(0))


def test_an_empty_draw_returns_an_empty_index_array():
    idx = systematic_sample(np.zeros(5), np.random.default_rng(0))
    assert idx.size == 0


# ----------------------------------------------------- the two composed together

def test_power_weighting_prefers_detectable_variants():
    """The point of the whole module: the drawn set is enriched for power.

    The pool is deliberately shaped like the real human one -- MAF and |selco|
    both log-uniform over several decades, so most of it is undetectable at
    n = 8000 and only a tail is findable.
    """
    rng = np.random.default_rng(99)
    maf = 10.0 ** rng.uniform(-4.0, np.log10(0.5), 5000)
    selco = 10.0 ** rng.uniform(-6.0, -2.0, 5000)
    power = detection_power(maf, np.sqrt(selco) * 10.0, 8000, 5e-8)
    # Strongly bimodal, like the real pool: the typical variant is invisible even
    # though a quarter of them are findable.
    assert np.median(power) < 0.01
    assert (power >= 0.05).mean() < 0.3

    idx = systematic_sample(inclusion_probabilities(power, 50), rng)
    assert np.median(power[idx]) > 0.99
    assert np.mean(power[idx]) > 4 * np.mean(power)
