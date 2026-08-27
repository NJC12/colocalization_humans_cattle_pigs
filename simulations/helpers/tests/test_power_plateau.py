"""Saturating power sampling: weight = min(power, plateau).

The premise is about real data. Past some detection probability a study finds the
variant regardless, so relative power above that point should not decide which
loci become causal -- the raw proportional rule disagrees, making a power-0.99
variant 1.24x likelier to be drawn than a power-0.80 one.

What has to hold, and what silently would not:

  * the default must be a no-op, or every power-sampled arm already on disk
    becomes unreproducible;
  * below the plateau the draw must stay PROPORTIONAL to power, not merely
    monotone -- `inclusion_probabilities` is scale-invariant, so this is a real
    property and not a tautology;
  * at or above the plateau every variant must carry EQUAL weight, which is the
    whole claim;
  * the pool-eligibility guard must keep using raw power. Flattening the top
    first would let a pool of near-certainties look like a pool of marginal
    variants and pass a guard designed to reject exactly that;
  * a saturated run must not share a directory with its unsaturated twin. They
    differ in which loci are causal, so a shared path silently mixes two
    simulations.
"""

import numpy as np
import pytest

from helpers import causal_power as cp
from helpers import params_record, paths


# ------------------------------------------------------------- the transform

def test_the_default_is_a_no_op():
    """Power is already bounded by 1, so plateau=1 must return it untouched.
    Every arm predating this knob relies on that."""
    p = np.array([1e-9, 0.05, 0.5, 0.8, 1.0])
    assert np.array_equal(cp.saturated_weights(p), p)
    assert np.array_equal(cp.saturated_weights(p, cp.NO_PLATEAU), p)
    assert cp.NO_PLATEAU == paths.NO_POWER_PLATEAU == 1.0


def test_below_the_plateau_the_weight_is_untouched():
    p = np.array([0.001, 0.1, 0.25, 0.4999])
    assert np.array_equal(cp.saturated_weights(p, 0.5), p)


def test_at_and_above_the_plateau_every_weight_is_equal():
    p = np.array([0.5, 0.6, 0.8, 0.99, 1.0])
    w = cp.saturated_weights(p, 0.5)
    assert np.all(w == 0.5)


@pytest.mark.parametrize("bad", [0, -0.1, 1.5, np.inf])
def test_a_plateau_outside_zero_to_one_is_refused(bad):
    with pytest.raises(ValueError):
        cp.saturated_weights(np.array([0.5]), bad)


def test_negative_or_nonfinite_power_is_refused():
    for bad in (np.array([-0.1]), np.array([np.nan])):
        with pytest.raises(ValueError):
            cp.saturated_weights(bad, 0.5)


# --------------------------------------------- what it does to the actual draw

def test_the_sub_plateau_draw_stays_proportional_to_power():
    """The claim is "proportional below, flat above" -- not just "monotone".

    inclusion_probabilities is scale-invariant, so ratios among the unsaturated
    units must survive the transform exactly.
    """
    power = np.array([0.02, 0.04, 0.08, 0.9, 0.95, 1.0])
    pi = cp.inclusion_probabilities(cp.saturated_weights(power, 0.5), 2)
    lo = pi[:3]
    assert lo[1] / lo[0] == pytest.approx(2.0)
    assert lo[2] / lo[0] == pytest.approx(4.0)


def test_the_high_power_tail_stops_being_favoured_among_itself():
    power = np.array([0.02, 0.04, 0.08, 0.9, 0.95, 1.0])
    plain = cp.inclusion_probabilities(power, 2)
    sat = cp.inclusion_probabilities(cp.saturated_weights(power, 0.5), 2)
    # Unsaturated: the three near-certain variants are ranked against each other.
    assert plain[3] < plain[4] < plain[5]
    # Saturated: they are indistinguishable.
    assert sat[3] == pytest.approx(sat[4]) == pytest.approx(sat[5])
    # And the low-power tail gains, because the budget is no longer skewed.
    assert sat[0] > plain[0]


def test_inclusion_probabilities_still_sum_to_the_sample_size():
    """The pi <= 1 clamp is a different mechanism and must survive the change."""
    power = np.array([0.01, 0.2, 0.5, 0.9, 1.0, 1.0, 1.0])
    for plateau in (1.0, 0.5, 0.1):
        pi = cp.inclusion_probabilities(cp.saturated_weights(power, plateau), 3)
        assert pi.sum() == pytest.approx(3.0)
        assert np.all(pi <= 1.0 + 1e-12)


def test_the_pool_guard_reads_raw_power_not_the_saturated_weight():
    """select_central_power counts variants at `>= min_power` BEFORE saturating.

    A plateau below min_power would otherwise lift the whole pool over the bar,
    turning the guard against degenerate draws into a rubber stamp.
    """
    import inspect
    from pathlib import Path
    src = Path(paths.__file__).parent.parent / "create_gwas_files_and_phenotypes.py"
    body = src.read_text()
    i_guard = body.index("n_above = int((power >= min_power).sum())")
    i_sat = body.index("weights = causal_power.saturated_weights(")
    assert i_guard < i_sat, "the guard must be evaluated on raw power, before saturation"
    assert "(power >= min_power)" in body


# ----------------------------------------------------------------- the paths

def test_a_saturated_run_gets_its_own_directory():
    base = dict(causal_sampling="power", sampling_gwas_n=50000,
                stage1_seed=11, pipeline="human", causal_min_maf=0)
    plain = paths.causal_sampling_segment(base)
    sat = paths.causal_sampling_segment({**base, "sampling_power_plateau": 0.5})
    assert plain == ".psamp_50000"
    assert sat == ".psamp_50000_sat05"
    assert plain != sat
    assert paths.stage2_run_tag({**base, "sampling_power_plateau": 0.5}).endswith("_sat05")


def test_the_default_plateau_leaves_every_existing_path_unchanged():
    base = dict(causal_sampling="power", sampling_gwas_n=8000,
                stage1_seed=11, pipeline="human", causal_min_maf=0)
    assert paths.causal_sampling_segment(base) == ".psamp_8000"
    assert paths.causal_sampling_segment({**base, "sampling_power_plateau": 1.0}) == ".psamp_8000"


def test_uniform_sampling_emits_no_segment_whatever_the_plateau():
    cfg = dict(causal_sampling="uniform", sampling_gwas_n=8000,
               sampling_power_plateau=0.5, stage1_seed=11)
    assert paths.causal_sampling_segment(cfg) == ""


@pytest.mark.parametrize("value,token", [(0.5, "05"), (0.25, "025"),
                                         (0.75, "075"), (0.9, "09")])
def test_the_token_never_contains_a_dot(value, token):
    """stage2_run_tag is dot-joined and summarize_coloc prefix-globs it."""
    got = paths.power_plateau_token(value)
    assert got == token and "." not in got


# ------------------------------------------------------------- the provenance

def test_the_plateau_is_a_stage2_determining_key():
    """It changes which variants are causal, so two runs differing only in it
    must hash to different stage2_uids -- otherwise the Snakefile's provenance
    guard would let one adopt the other's phenotypes."""
    assert "sampling_power_plateau" in params_record.STAGE2_KEYS
    base = {k: 1 for k in params_record.STAGE2_KEYS}
    a = dict(base, sampling_power_plateau=1.0)
    b = dict(base, sampling_power_plateau=0.5)
    assert params_record.uid(a, params_record.STAGE2_KEYS) != \
           params_record.uid(b, params_record.STAGE2_KEYS)
