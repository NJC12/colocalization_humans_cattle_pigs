"""Unit tests for helpers/pval_rescale.py.

These check the algebra and the guard rails.  The statistical CLAIM -- that the
predicted distribution is calibrated against real regressions -- is measured
separately by helpers/validate_pval_rescale.py, not asserted here.
"""

import warnings

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from helpers import pval_rescale as pr


# --------------------------------------------------------------------------
# p <-> z round trips, including the underflow regime
# --------------------------------------------------------------------------

@pytest.mark.parametrize("z", [0.5, 1.96, 5.45, 30.0, 100.0, 400.0])
def test_z_neglog10p_round_trip(z):
    back = pr.neglog10p_to_z(pr.z_to_neglog10p(z))
    assert back == pytest.approx(z, rel=1e-9)


def test_p_to_z_matches_scipy_in_the_ordinary_range():
    p = np.array([0.5, 0.05, 5e-8, 1e-100])
    assert pr.p_to_z(p) == pytest.approx(stats.norm.isf(p / 2), rel=1e-9)


def test_p_to_z_survives_subnormal_pvalues():
    """A p-value at the float64 floor must still give a finite statistic."""
    z = pr.p_to_z(1e-320)
    assert np.isfinite(z) and z > 38.0


def test_z_to_neglog10p_never_underflows():
    """The whole point of working in -log10 space."""
    v = pr.z_to_neglog10p(500.0)
    assert np.isfinite(v) and v > 50_000


def test_p_to_z_with_df_uses_the_t_distribution():
    z_norm = pr.p_to_z(0.001)
    t_small = pr.p_to_z(0.001, df=10)
    assert t_small > z_norm  # heavier tails need a bigger statistic
    assert pr.p_to_z(0.001, df=1e6) == pytest.approx(z_norm, rel=1e-3)


def test_p_to_z_rejects_out_of_range():
    with pytest.raises(ValueError):
        pr.p_to_z(1.5)


# --------------------------------------------------------------------------
# The point estimate
# --------------------------------------------------------------------------

def test_same_sample_size_is_the_identity():
    z1 = -4.2
    assert pr.rescale_z(z1, 8000, 8000) == pytest.approx(z1)
    out = pr.predict_interval(z1, 8000, 8000)
    assert out["z2"] == pytest.approx(z1)
    assert out["z2_lo"] == pytest.approx(z1)   # zero-width: t == 1 => sd == 0
    assert out["z2_hi"] == pytest.approx(z1)


def test_sign_is_preserved():
    assert pr.rescale_z(-4.0, 8000, 2000) < 0
    assert pr.rescale_z(+4.0, 8000, 2000) > 0


def test_scaling_follows_sqrt_of_the_information_fraction():
    z1, n1, n2 = 6.0, 8000, 2000
    t = pr.information_fraction(n1, n2)
    assert pr.rescale_z(z1, n1, n2) == pytest.approx(z1 * np.sqrt(t))


def test_r2_and_ncp_agree_exactly_when_maf_is_unchanged():
    """Not merely for small effects -- the (1 - r2) factor cancels identically.

    This is the claim in the module docstring, checked across four orders of
    magnitude of z^2 / N including a variant explaining ~25% of variance.
    """
    n1, n2 = 5000, 40000
    for z1 in (1.0, 5.0, 20.0, 40.0):
        a = pr.rescale_z(z1, n1, n2, method="r2")
        b = pr.rescale_z(z1, n1, n2, method="ncp")
        assert a == pytest.approx(b, rel=1e-12)
    # and with MAF supplied but identical on both sides
    assert pr.rescale_z(40.0, n1, n2, method="r2", maf1=0.2, maf2=0.2) == (
        pytest.approx(pr.rescale_z(40.0, n1, n2, method="ncp"), rel=1e-12)
    )


def test_r2_route_diverges_only_when_maf_changes():
    """Higher heterozygosity at the same N means a larger statistic."""
    z2 = pr.rescale_z(6.0, 5000, 5000, method="r2", maf1=0.05, maf2=0.5)
    assert z2 > 6.0
    # and the direction reverses
    z2 = pr.rescale_z(6.0, 5000, 5000, method="r2", maf1=0.5, maf2=0.05)
    assert z2 < 6.0


def test_ncp_method_refuses_a_maf_change():
    with pytest.raises(ValueError, match="cannot represent a change in allele"):
        pr.rescale_z(6.0, 5000, 5000, method="ncp", maf1=0.1, maf2=0.4)


def test_maf_rescale_refuses_impossible_variance_explained():
    with pytest.raises(ValueError, match="variance explained"):
        pr.rescale_z(90.0, 100, 100, method="r2", maf1=0.01, maf2=0.5)


def test_monotone_in_target_sample_size():
    z = [abs(pr.rescale_z(4.0, 8000, n2)) for n2 in (500, 1000, 4000, 8000, 50_000)]
    assert all(np.diff(z) > 0)
    lp = [pr.predict_interval(4.0, 8000, n2)["neglog10_p2"]
          for n2 in (500, 1000, 4000, 8000, 50_000)]
    assert all(np.diff(lp) > 0)


def test_no_underflow_scaling_up_to_human_sample_sizes():
    out = pr.predict_interval(9.0, 8000, 500_000, design="independent")
    assert np.isfinite(out["neglog10_p2"]) and out["neglog10_p2"] > 1000


def test_requires_both_or_neither_maf():
    with pytest.raises(ValueError, match="both maf1 and maf2"):
        pr.rescale_z(4.0, 8000, 2000, maf1=0.2)


def test_rejects_sample_size_below_df_floor():
    with pytest.raises(ValueError, match="non-positive residual df"):
        pr.rescale_z(4.0, 8000, 10, k_covar=20)


# --------------------------------------------------------------------------
# The interval: which design, and the drift-free bridge
# --------------------------------------------------------------------------

def test_nested_variance_matches_the_brownian_bridge_downward():
    """Var = 1 - t for t < 1."""
    n1, n2 = 8000, 2000
    t = pr.information_fraction(n1, n2)
    out = pr.predict_interval(5.0, n1, n2, design="nested", level=0.95)
    half = (out["z2_hi"] - out["z2_lo"]) / 2
    expected = stats.norm.isf(0.025) * np.sqrt(1 - t)
    assert half == pytest.approx(expected, rel=1e-12)


def test_nested_variance_matches_the_predictive_increment_upward():
    """Var = t - 1 for t > 1: the increment PLUS the plug-in drift's own variance.

    Not (t-1)/t, which assumes the drift is known.  Measured empirically by
    helpers/validate_pval_rescale.py.
    """
    n1, n2 = 2000, 8000
    t = pr.information_fraction(n1, n2)
    out = pr.predict_interval(5.0, n1, n2, design="nested", level=0.95)
    half = (out["z2_hi"] - out["z2_lo"]) / 2
    assert half == pytest.approx(stats.norm.isf(0.025) * np.sqrt(t - 1), rel=1e-12)


def test_nested_variance_is_abs_t_minus_one_in_both_directions():
    for n2 in (500, 2000, 4000, 16_000, 64_000):
        t = pr.information_fraction(8000, n2)
        out = pr.predict_interval(5.0, 8000, n2, design="nested")
        half = (out["z2_hi"] - out["z2_lo"]) / 2
        expected = stats.norm.isf(0.025) * np.sqrt(abs(t - 1))
        assert half == pytest.approx(expected, rel=1e-12)


def test_conditional_design_matches_nested_downward_and_is_narrower_upward():
    """The bridge needs no drift, so the two coincide for t <= 1."""
    def half(design, n2):
        o = pr.predict_interval(5.0, 8000, n2, design=design)
        return (o["z2_hi"] - o["z2_lo"]) / 2

    assert half("nested_conditional", 2000) == pytest.approx(half("nested", 2000))
    assert half("nested_conditional", 32_000) < half("nested", 32_000)


def test_nested_variance_vanishes_at_t_equals_one_and_grows_away_from_it():
    widths = [pr.predict_interval(5.0, 8000, n2)["z2_hi"]
              - pr.predict_interval(5.0, 8000, n2)["z2_lo"]
              for n2 in (8000, 6000, 2000, 500)]
    assert widths[0] == pytest.approx(0.0, abs=1e-12)
    assert all(np.diff(widths) > 0)


def test_independent_design_is_wider_than_nested():
    kw = dict(z1=5.0, n1=8000, n2=2000)
    def width(design):
        o = pr.predict_interval(kw["z1"], kw["n1"], kw["n2"], design=design)
        return o["z2_hi"] - o["z2_lo"]
    assert width("deterministic") == pytest.approx(0.0)
    assert width("nested") < width("independent_conditional") < width("independent")


def test_unknown_design_rejected():
    with pytest.raises(ValueError, match="design must be one of"):
        pr.predict_interval(5.0, 8000, 2000, design="bootstrap")


def test_interval_straddling_zero_reports_no_association_at_the_low_end():
    out = pr.predict_interval(0.3, 8000, 200, design="nested")
    assert out["z2_lo"] < 0 < out["z2_hi"]
    assert out["neglog10_p2_lo"] == pytest.approx(0.0)


def test_neglog10p_ordering_within_the_interval():
    out = pr.predict_interval(6.0, 8000, 3000, design="nested")
    assert out["neglog10_p2_lo"] <= out["neglog10_p2"] <= out["neglog10_p2_hi"]


def test_level_rejected_outside_open_unit_interval():
    with pytest.raises(ValueError, match="level must lie"):
        pr.predict_interval(5.0, 8000, 2000, level=1.0)


# --------------------------------------------------------------------------
# Winner's curse
# --------------------------------------------------------------------------

def _z_vector(seed=0, n=10_000, n_signal=50):
    rng = np.random.default_rng(seed)
    z = rng.normal(size=n)
    z[:n_signal] += rng.normal(6.0, 1.0, size=n_signal)
    return z


def test_fiqt_shrinks_toward_zero_and_is_monotone():
    z = _z_vector()
    lam = pr.fiqt(z)
    assert np.all(np.abs(lam) <= np.abs(z) + 1e-9)
    assert np.all(np.sign(lam[lam != 0]) == np.sign(z[lam != 0]))
    order = np.argsort(z)
    assert np.all(np.diff(lam[order]) >= -1e-9)  # monotone in z


def test_fiqt_shrinks_null_statistics_hardest():
    """Relative, not absolute: a null at |z| ~ 0.8 cannot shrink by more than 0.8.

    Nulls should collapse to essentially nothing while true signal keeps most of
    its magnitude -- that asymmetry is what makes the correction useful.
    """
    z = _z_vector()
    lam = pr.fiqt(z)
    null_frac = np.mean(1 - np.abs(lam[1000:]) / np.abs(z[1000:]))
    signal_frac = np.mean(1 - np.abs(lam[:50]) / np.abs(z[:50]))
    assert null_frac > 0.9
    assert signal_frac < 0.35
    assert null_frac > signal_frac


def test_fiqt_refuses_a_single_variant():
    with pytest.raises(ValueError, match="genome-wide z vector"):
        pr.fiqt(np.array([5.0]))


def test_winners_curse_rejected_for_nested_down_sampling():
    """The bridge is drift-free, so shrinking would bias the prediction."""
    frame = pd.DataFrame({"BETA": [0.1, 0.2, 0.3], "SE": [0.02, 0.02, 0.02]})
    with pytest.raises(ValueError, match="drift-free"):
        pr.rescale_summary_stats(frame, n1=8000, n2=2000,
                                 design="nested", winners_curse=True)


def test_up_scaling_warns_that_selection_is_not_modelled():
    frame = pd.DataFrame({"BETA": [0.1, 0.2], "SE": [0.02, 0.02]})
    with pytest.warns(RuntimeWarning, match="extrapolating UP"):
        pr.rescale_summary_stats(frame, n1=2000, n2=8000, design="nested")


def test_down_scaling_does_not_warn():
    frame = pd.DataFrame({"BETA": [0.1, 0.2], "SE": [0.02, 0.02]})
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        pr.rescale_summary_stats(frame, n1=8000, n2=2000, design="nested")


@pytest.mark.filterwarnings("ignore:extrapolating UP")
def test_winners_curse_allowed_and_lowers_the_prediction_when_scaling_up():
    z = _z_vector()
    frame = pd.DataFrame({"BETA": z * 0.02, "SE": np.full(z.size, 0.02)})
    plain = pr.rescale_summary_stats(frame, n1=2000, n2=8000, design="nested")
    shrunk = pr.rescale_summary_stats(frame, n1=2000, n2=8000, design="nested",
                                      winners_curse=True)
    top = np.argsort(-np.abs(plain["z1"].to_numpy()))[:50]
    assert np.all(np.abs(shrunk["z2"].to_numpy()[top])
                  <= np.abs(plain["z2"].to_numpy()[top]) + 1e-9)
    assert np.mean(np.abs(shrunk["z2"].to_numpy()[top])) < np.mean(
        np.abs(plain["z2"].to_numpy()[top]))


# --------------------------------------------------------------------------
# Frame-level behaviour: input precedence and the MAC guard
# --------------------------------------------------------------------------

def _frame():
    return pd.DataFrame({
        "ID": ["a", "b", "c", "d"],
        "BETA": [0.20, -0.15, 0.30, 0.05],
        "SE": [0.02, 0.02, 0.02, 0.02],
        "A1_FREQ": [0.40, 0.25, 0.002, 0.10],
        "P": [1e-30, 1e-14, 1e-50, 0.012],
    })


def test_beta_se_is_preferred_over_p_and_keeps_the_sign():
    out = pr.rescale_summary_stats(_frame(), n1=8000, n2=2000, p_col="P")
    assert (out["z1_source"] == "beta/se").all()
    assert out["z1"].iloc[1] < 0


def test_p_column_fallback_warns_about_the_lost_sign():
    frame = _frame().drop(columns=["BETA", "SE"])
    with pytest.warns(RuntimeWarning, match="sign of the effect is unknown"):
        out = pr.rescale_summary_stats(frame, n1=8000, n2=2000, p_col="P")
    assert (out["z1"] > 0).all()
    assert (out["z1_source"] == "p").all()


def test_missing_all_input_columns_is_an_error():
    with pytest.raises(ValueError, match="need either beta"):
        pr.rescale_summary_stats(pd.DataFrame({"ID": ["a"]}), n1=8000, n2=2000)


def test_low_mac_rows_are_withheld_not_guessed():
    """maf=0.002 at n2=2000 gives an expected MAC of 8 -- below the floor."""
    out = pr.rescale_summary_stats(_frame(), n1=8000, n2=2000,
                                   maf_col="A1_FREQ", min_mac=20)
    assert out["expected_mac_n2"].iloc[2] == pytest.approx(8.0)
    assert not out["rescalable"].iloc[2]
    assert np.isnan(out["neglog10_p2"].iloc[2])
    assert out["rescalable"].drop(index=2).all()
    assert out.loc[out["rescalable"], "neglog10_p2"].notna().all()


def test_mac_guard_uses_the_target_sample_size():
    """The same variant clears the floor at a larger n2."""
    out = pr.rescale_summary_stats(_frame(), n1=8000, n2=8000,
                                   maf_col="A1_FREQ", min_mac=20)
    assert out["rescalable"].all()


def test_no_maf_column_means_no_guard():
    out = pr.rescale_summary_stats(_frame(), n1=8000, n2=2000)
    assert out["rescalable"].all()
    assert out["expected_mac_n2"].isna().all()


def test_zero_standard_error_rejected():
    frame = _frame()
    frame.loc[0, "SE"] = 0.0
    with pytest.raises(ValueError, match="non-positive standard errors"):
        pr.rescale_summary_stats(frame, n1=8000, n2=2000)
