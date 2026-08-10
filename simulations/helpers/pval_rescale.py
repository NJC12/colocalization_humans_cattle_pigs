#!/usr/bin/env python
"""Rescale a continuous-trait GWAS association to a different sample size.

Given a per-variant test statistic measured at sample size N1, predict what the
statistic -- and hence the p-value -- would be at sample size N2.

Definitions
-----------
z             signed Wald statistic, beta / SE.  Preferred input.  When only a
              two-sided p-value is available, |z| is recovered from it and the
              sign is unknown (reported as +).
df            residual degrees of freedom, n - k_covar - 2, where `k_covar`
              counts covariates EXCLUDING the intercept and the tested variant.
              Everything below is expressed in df rather than n so that
              `n2 == n1` is an exact identity and small-n cases stay honest.
t             information fraction, df2 / df1.  Optionally multiplied by a MAF
              factor -- see `rescale_z`.
r2            partial variance explained by the variant, t_stat^2 / (t_stat^2 + df).
              The identity `lambda = n * r2 / (1 - r2)` is what makes the
              rescaling exact rather than a small-effect approximation.

The scaling law
---------------
For a quantitative trait under linear regression, the Wald statistic satisfies
z^2 ~ chi2_1(lambda) with non-centrality

    lambda = (beta / SE)^2 ~= N * 2f(1-f) beta^2 / Var(y) = N * r2 / (1 - r2)

so lambda is LINEAR IN N and |z| scales as sqrt(N).  Sham & Purcell 2014
(Nat Rev Genet 15:335-346, https://www.nature.com/articles/nrg3706);
Wang & Xu 2019 (Heredity 123:287-306,
https://www.nature.com/articles/s41437-019-0205-3).  The same linearity is what
LD-score regression exploits genome-wide, E[chi2_j] = 1 + N*a + (N h2 / M) l_j
(Bulik-Sullivan et al. 2015, Nat Genet 47:291-295,
https://www.nature.com/articles/ng.3211).

A NOTE ON "EXACTNESS".  It is tempting to say the sqrt(N) rule is a small-effect
approximation that an r2-based route repairs.  That is not right when MAF is
held fixed: substituting the plug-in r2_hat = t1^2 / (t1^2 + df1) into
lambda = df * r2 / (1 - r2) makes the (1 - r2) factor cancel, leaving

    t2 = t1 * sqrt(df2 / df1)

for ANY effect size.  The r2 route earns its keep only when the ALLELE FREQUENCY
differs between the two settings, because then r2 itself changes
(r2 ~ 2f(1-f) beta^2) and (1 - r2) no longer cancels.  `method="r2"` is
therefore required for maf1 != maf2 and is otherwise identical to
`method="ncp"`.  The residual (1 - r2) bookkeeping is the same algebra behind
the summary-statistic conversion b = z / sqrt(2f(1-f)(n + z^2)) used by SMR
(Zhu et al. 2016, Nat Genet 48:481-487, https://www.nature.com/articles/ng.3538)
-- note the `n + z^2` rather than `n`.

Uncertainty: which design are you in?
-------------------------------------
The point estimate is z1 * sqrt(t) in every case; the designs differ only in the
SPREAD, and the spread is what decides whether "this locus would be significant
at N2" is a claim or a coin flip.  Model the statistic as a Brownian motion in
information time, B(t) = z(t) sqrt(t) -- the "B-value" of Lan & Wittes 1988
(Biometrics 44:579-585, https://pubmed.ncbi.nlm.nih.gov/3390511/); textbook
treatment in Proschan, Lan & Wittes 2006
(https://link.springer.com/book/10.1007/978-0-387-44970-8).  This is the
conditional-power machinery of group-sequential trial monitoring, which asks
exactly our question.

  design="nested"        N2 individuals are drawn from (t<1) or contain (t>1)
                         the same cohort.
                             z2 | z1 ~ N(z1 sqrt(t), |t - 1|)
                         The two directions reach the same variance by different
                         routes:
                           t < 1  Brownian BRIDGE.  The drift cancels EXACTLY --
                                  Var = 1 - t holds whatever the true effect is.
                                  No prior, no assumed effect size, and -- see
                                  below -- no winner's curse.
                           t > 1  the increment is independent of z1 but DOES
                                  carry the drift, so the plug-in theta_hat = z1
                                  contributes its own uncertainty:
                                  (t-1)^2/t + (t-1)/t = t - 1.  This is the
                                  predictive variance under a flat prior on the
                                  drift; it is unbounded as t grows, whereas the
                                  drift-known value (t-1)/t is not.  Getting this
                                  wrong understates the spread badly -- at t = 4
                                  it is 3 versus 0.75.

  design="nested_conditional"
                         classical conditional power: the same nested setup but
                         treating the drift as GIVEN rather than estimated,
                             z2 | z1 ~ N(z1 sqrt(t), 1 - min(t, 1/t))
                         Identical to "nested" for t <= 1 (the bridge needs no
                         drift).  For t > 1 it is the narrower, optimistic
                         interval -- use it only when the effect size is genuinely
                         known a priori rather than read off this same data.

  design="independent"   a hypothetical separate study of size N2.
                             z2 | z1 ~ N(z1 sqrt(t), 1 + t)     (predictive,
                         flat prior on the true non-centrality) -- the
                         predictive-power / assurance version, O'Hagan, Stevens &
                         Campbell 2005, Pharm Stat 4:187-201,
                         https://doi.org/10.1002/pst.175.  Pass
                         design="independent_conditional" for the variance-1
                         version that conditions on the true effect being known.

  design="deterministic" no uncertainty: the point estimate alone.  Honest only
                         as "what p-value does this ESTIMATED effect give at N2",
                         never as a prediction.

Winner's curse applies in ONE direction only
--------------------------------------------
If the variant was selected for being significant, |z1| is upward-biased for the
true non-centrality.  That bias propagates only when the extrapolation depends on
an ESTIMATE of the true effect:

  * t < 1, design="nested": NO correction.  The bridge is drift-free, so nothing
    about the true effect enters.  Shrinking z1 here would actively bias the
    prediction downward.  `rescale_summary_stats` refuses it.
  * t > 1, or design="independent": the uncorrected extrapolation IS
    anti-conservative.  `fiqt` implements the FDR inverse quantile transform of
    Bigdeli et al. 2016 (Bioinformatics 32:2598-2603,
    https://academic.oup.com/bioinformatics/article/32/17/2598/2450747), whose
    stated purpose is predicting which signals appear in much larger scans.
    Alternatives: conditional likelihood (Zhong & Prentice 2008, Biostatistics
    9:621-634, https://academic.oup.com/biostatistics/article/9/4/621/258822) and
    the survey + R package of Forde, Hemani & Ferguson 2023 (PLoS Genet
    19:e1010546, https://amandaforde.github.io/winnerscurse/).

    BUT: it is OFF by default, because in this project's regime it makes things
    worse, not better.  See the measured numbers below -- with sparse signal (a
    handful of causal variants in a regional scan) the FDR step is very
    conservative and FIQT OVER-shrinks, turning a +1.1 z over-prediction into a
    -1.7 z under-prediction and lowering coverage further.  Bigdeli et al. report
    FIQT at its best with many signals and large samples; a 2 Mb window with a
    dozen causal variants is the opposite corner.  Measure before enabling it.

Measured calibration
--------------------
From helpers/validate_pval_rescale.py (8000 individuals, 3000 variants with LD,
5 covariates, 40 nested draws per target; regenerate with `--out`):

    leg   t        subset               slope   resid_sd  theory_sd  coverage
    down  0.4996   pass_mac_guard       1.004   0.7053    0.7074     0.9502
    down  0.2493   pass_mac_guard       1.009   0.8666    0.8664     0.9497
    down  0.1242   pass_mac_guard       0.9931  0.9413    0.9358     0.9486
    down  0.0617   pass_mac_guard       0.9974  0.9693    0.9687     0.9497
    up    4.011    all_variants         0.4031  1.736     1.735      0.9498
    up    4.011    selected, no corr.   --      1.980     1.735      0.8456
    up    4.011    selected, FIQT       --      2.041     1.735      0.7702

So: DOWN-sampling is calibrated to three decimal places at every size, which is
the design this project uses.  UP-sampling is calibrated on unselected variants
and degrades once variants are chosen for being significant -- the interval
models sampling noise, not the selection, and no correction tried here repairs
it.  Treat up-scaling of a hand-picked hit as indicative only.  (The `up` slope
of 0.40 is regression dilution, not miscalibration: the PREDICTOR z1 is itself
noisy at t = 4.  Coverage is the metric to read there.)

What this does NOT fix
----------------------
The transform holds beta, allele frequency, Var(y) and the model constant and
changes only N.  It does not repair differing LD (a lead-SNP p-value is a tagging
statistic), differing trait definitions, or confounding -- rescale
genomic-control-corrected statistics if lambda_GC > 1.  And "sample size" must
mean EFFECTIVE sample size: livestock GWAS on de-regressed EBVs or daughter yield
deviations carry an effective record count, not the animal count (Garrick, Taylor
& Fernando 2009, GSE 41:55,
https://gsejournal.biomedcentral.com/articles/10.1186/1297-9686-41-55).  Feeding
a raw animal count in overstates the information available.

Everything is returned in -log10(p).  At human sample sizes a bare p-value
underflows float64 (p < 1e-308) almost immediately, so no function here returns
one.

CLI
---
    pval_rescale.py --sumstats stats.tsv --n1 400000 --n2 8000 \
        --beta-col BETA --se-col SE --maf-col A1_FREQ --design nested \
        --out rescaled.tsv
"""

import argparse
import sys
import warnings

import numpy as np
from scipy import special, stats

__all__ = [
    "p_to_z",
    "neglog10p_to_z",
    "z_to_neglog10p",
    "rescale_z",
    "predict_interval",
    "fiqt",
    "rescale_summary_stats",
]

# Expected minor-allele count below which the normal approximation to the Wald
# statistic is not trustworthy and a rescaled p-value should not be emitted.
DEFAULT_MIN_MAC = 20.0

DESIGNS = ("nested", "nested_conditional", "independent",
           "independent_conditional", "deterministic")

_LOG10 = np.log(10.0)


def _resid_df(n, k_covar):
    """Residual df for testing one variant with `k_covar` extra covariates.

    The fitted model is y ~ 1 + g + (k_covar covariates), so df = n - k_covar - 2.
    """
    df = np.asarray(n, dtype=float) - float(k_covar) - 2.0
    if np.any(df <= 0):
        raise ValueError(
            f"non-positive residual df: n={n}, k_covar={k_covar} leaves "
            f"df={df}. Sample size must exceed k_covar + 2."
        )
    return df


def p_to_z(p, df=None):
    """Two-sided p-value -> |z| (or |t| when `df` is given).

    Uses `scipy.special.ndtri_exp` on the log scale so that p-values far below
    the float64 subnormal limit still map to a finite statistic.  Returns the
    MAGNITUDE only -- a p-value carries no sign.  Prefer `beta / SE` when it is
    available.
    """
    p = np.asarray(p, dtype=float)
    if np.any((p < 0) | (p > 1)):
        raise ValueError("p-values must lie in [0, 1]")
    if df is None:
        # Phi(-z) = p/2  =>  z = -ndtri_exp(log(p/2))
        with np.errstate(divide="ignore"):
            return -special.ndtri_exp(np.log(p) - np.log(2.0))
    return stats.t.isf(p / 2.0, _as_df(df))


def neglog10p_to_z(neglog10_p, df=None):
    """-log10(two-sided p) -> |z| (or |t| when `df` is given).

    The entry point to use when a summary-statistic file stores -log10 p rather
    than p (as Pan-UKB and plink2 `--glm ... cols=+nlog10p` do), since the
    round-trip through a bare p-value would underflow.
    """
    neglog10_p = np.asarray(neglog10_p, dtype=float)
    if np.any(neglog10_p < 0):
        raise ValueError("-log10(p) must be non-negative")
    log_half_p = -neglog10_p * _LOG10 - np.log(2.0)
    if df is None:
        return -special.ndtri_exp(log_half_p)
    # No log-scale inverse survival for the t distribution; fall back, and warn
    # where it will have underflowed.
    p = np.exp(log_half_p) * 2.0
    if np.any(p == 0.0):
        warnings.warn(
            "t-distribution inverse underflowed for very small p; "
            "pass df=None to use the normal approximation on the log scale",
            RuntimeWarning,
            stacklevel=2,
        )
    return stats.t.isf(p / 2.0, _as_df(df))


def z_to_neglog10p(z, df=None):
    """Signed z (or t) -> -log10(two-sided p), computed on the log scale."""
    z = np.abs(np.asarray(z, dtype=float))
    if df is None:
        log_sf = stats.norm.logsf(z)
    else:
        log_sf = stats.t.logsf(z, _as_df(df))
    return -(log_sf + np.log(2.0)) / _LOG10


def _as_df(df):
    df = np.asarray(df, dtype=float)
    if np.any(df <= 0):
        raise ValueError("df must be positive")
    return df


def _maf_var_factor(maf1, maf2):
    """2 f2 (1 - f2) / (2 f1 (1 - f1)), the genotype-variance ratio."""
    f1 = np.asarray(maf1, dtype=float)
    f2 = np.asarray(maf2, dtype=float)
    if np.any((f1 <= 0) | (f1 >= 1) | (f2 <= 0) | (f2 >= 1)):
        raise ValueError("allele frequencies must lie strictly in (0, 1)")
    return (f2 * (1.0 - f2)) / (f1 * (1.0 - f1))


def information_fraction(n1, n2, k_covar=0):
    """df2 / df1 -- the information fraction `t` used throughout."""
    return _resid_df(n2, k_covar) / _resid_df(n1, k_covar)


def rescale_z(z1, n1, n2, k_covar=0, method="r2", maf1=None, maf2=None):
    """Point estimate of the statistic at sample size `n2`.

    Parameters
    ----------
    z1 : signed Wald statistic (beta / SE) observed at `n1`.  Sign is preserved.
    n1, n2 : sample sizes.
    k_covar : covariates excluding the intercept and the tested variant.
    method : "r2" (default) or "ncp".
        "ncp" applies z2 = z1 * sqrt(df2 / df1) directly.
        "r2" routes through the partial variance explained,
        r2 = t^2 / (t^2 + df), rescales it by the genotype-variance ratio when
        MAF changes, and inverts.  With maf1 == maf2 (or both None) the two are
        ALGEBRAICALLY IDENTICAL -- see the module docstring.  Only "r2" accepts
        a MAF change.
    maf1, maf2 : allele frequencies in the two settings.  Both or neither.

    Returns
    -------
    z2 : float or ndarray, signed.
    """
    if method not in ("r2", "ncp"):
        raise ValueError(f"method must be 'r2' or 'ncp', got {method!r}")

    z1 = np.asarray(z1, dtype=float)
    df1 = _resid_df(n1, k_covar)
    df2 = _resid_df(n2, k_covar)

    if (maf1 is None) != (maf2 is None):
        raise ValueError("supply both maf1 and maf2, or neither")
    maf_changes = maf1 is not None and not np.allclose(
        np.asarray(maf1, dtype=float), np.asarray(maf2, dtype=float)
    )
    if maf_changes and method == "ncp":
        raise ValueError(
            "method='ncp' cannot represent a change in allele frequency; "
            "use method='r2' (the (1 - r2) factor no longer cancels)"
        )

    if method == "ncp" or not maf_changes:
        return z1 * np.sqrt(df2 / df1)

    # r2 route with a genuine MAF change.  r2 ~ 2f(1-f) beta^2, so rescaling r2
    # by the genotype-variance ratio holds the standardized effect beta fixed.
    r2_1 = z1**2 / (z1**2 + df1)
    r2_2 = r2_1 * _maf_var_factor(maf1, maf2)
    if np.any(r2_2 >= 1.0):
        raise ValueError(
            "MAF rescaling implies variance explained >= 1; the observed effect "
            "is not compatible with the target allele frequency"
        )
    z2_sq = df2 * r2_2 / (1.0 - r2_2)
    return np.sign(z1) * np.sqrt(z2_sq)


def _design_sd(t, design):
    """SD of z2 | z1 under each design.  `t` is the information fraction."""
    t = np.asarray(t, dtype=float)
    if design == "nested":
        # t<1: Brownian bridge, Var = 1 - t, exact and drift-free.
        # t>1: increment (t-1)/t plus the variance of the plug-in drift
        #      ((t-1)/sqrt(t))^2 * 1, which sums to t - 1.
        return np.sqrt(np.abs(t - 1.0))
    if design == "nested_conditional":
        # Classical conditional power: drift treated as known.
        return np.sqrt(1.0 - np.minimum(t, 1.0 / t))
    if design == "independent":
        return np.sqrt(1.0 + t)  # predictive, flat prior on the true effect
    if design == "independent_conditional":
        return np.ones_like(t)
    if design == "deterministic":
        return np.zeros_like(t)
    raise ValueError(f"design must be one of {DESIGNS}, got {design!r}")


def predict_interval(z1, n1, n2, design="nested", level=0.95, k_covar=0,
                     method="r2", maf1=None, maf2=None):
    """Predicted statistic at `n2` with an interval, as a dict of arrays.

    Returns keys `z2_lo`, `z2`, `z2_hi` (the interval on the signed z scale) and
    `neglog10_p2_lo`, `neglog10_p2`, `neglog10_p2_hi` (on the significance
    scale, where "lo" is the LEAST significant end).  If the z interval straddles
    zero, `neglog10_p2_lo` is 0.0 -- the prediction admits no association at all.
    """
    if not 0.0 < level < 1.0:
        raise ValueError("level must lie strictly in (0, 1)")

    center = rescale_z(z1, n1, n2, k_covar=k_covar, method=method,
                       maf1=maf1, maf2=maf2)
    t = information_fraction(n1, n2, k_covar)
    sd = _design_sd(t, design)
    half = stats.norm.isf((1.0 - level) / 2.0) * sd

    z_lo, z_hi = center - half, center + half
    straddles_zero = (z_lo <= 0.0) & (z_hi >= 0.0)
    least = np.where(straddles_zero, 0.0, np.minimum(np.abs(z_lo), np.abs(z_hi)))
    most = np.maximum(np.abs(z_lo), np.abs(z_hi))

    return {
        "z2_lo": z_lo,
        "z2": center,
        "z2_hi": z_hi,
        "neglog10_p2_lo": z_to_neglog10p(least),
        "neglog10_p2": z_to_neglog10p(center),
        "neglog10_p2_hi": z_to_neglog10p(most),
    }


def fiqt(z, min_p=1e-300):
    """FDR inverse quantile transform: shrink z toward the true non-centrality.

    A faithful port of the reference R implementation in Bigdeli et al. 2016
    (Bioinformatics 32:2598-2603,
    https://academic.oup.com/bioinformatics/article/32/17/2598/2450747): convert
    to two-sided p, Benjamini-Hochberg adjust, and map back through the Gaussian
    quantile.  Statistics more extreme than `min_p` are passed through unchanged,
    as in the original.

    Requires the GENOME-WIDE vector, not a single variant -- the correction is
    borrowed across the whole set of tests.  Guarantees |lambda_hat| <= |z|.
    """
    z = np.asarray(z, dtype=float)
    if z.ndim == 0 or z.size < 2:
        raise ValueError(
            "fiqt needs the genome-wide z vector; the shrinkage is borrowed "
            "across all tests and is undefined for a single variant"
        )
    p = 2.0 * stats.norm.sf(np.abs(z))
    p = np.maximum(p, min_p)
    adj = stats.false_discovery_control(p, method="bh")
    lam = np.sign(z) * stats.norm.isf(adj / 2.0)
    # Preserve the extreme tail exactly, as the reference implementation does.
    extreme = np.abs(z) > stats.norm.isf(min_p / 2.0)
    lam = np.where(extreme, z, lam)
    return lam


def rescale_summary_stats(
    df,
    n1,
    n2,
    beta_col="BETA",
    se_col="SE",
    p_col=None,
    neglog10p_col=None,
    maf_col=None,
    maf2_col=None,
    k_covar=0,
    method="r2",
    design="nested",
    level=0.95,
    winners_curse=False,
    min_mac=DEFAULT_MIN_MAC,
):
    """Rescale a whole summary-statistic frame.  Returns a new DataFrame.

    Added columns: `z1`, `z1_source`, `z2`, `z2_lo`, `z2_hi`, `neglog10_p2`,
    `neglog10_p2_lo`, `neglog10_p2_hi`, `expected_mac_n2`, `rescalable`.
    Rows failing the MAC guard get NaN predictions and `rescalable=False` rather
    than a number that the normal approximation cannot support.
    """
    import pandas as pd  # local: keeps the numeric API pandas-free

    out = df.copy()

    # --- z1, preferring beta/SE over any p-value round-trip -----------------
    if beta_col in out.columns and se_col in out.columns:
        se = out[se_col].to_numpy(dtype=float)
        if np.any(se <= 0):
            raise ValueError(f"non-positive standard errors in {se_col!r}")
        z1 = out[beta_col].to_numpy(dtype=float) / se
        source = "beta/se"
    elif neglog10p_col is not None and neglog10p_col in out.columns:
        z1 = neglog10p_to_z(out[neglog10p_col].to_numpy(dtype=float))
        source = "neglog10p"
    elif p_col is not None and p_col in out.columns:
        z1 = p_to_z(out[p_col].to_numpy(dtype=float))
        source = "p"
    else:
        raise ValueError(
            "need either beta+se columns, or a p / -log10 p column; "
            f"frame has {list(out.columns)}"
        )
    if source != "beta/se":
        warnings.warn(
            f"recovered |z| from {source}; the sign of the effect is unknown "
            "and is reported as positive",
            RuntimeWarning,
            stacklevel=2,
        )

    t = float(information_fraction(n1, n2, k_covar))

    # --- winner's curse: legal only when the target needs an effect estimate -
    if t > 1.0 and design.startswith("nested"):
        warnings.warn(
            f"extrapolating UP (t={t:.3g}): the interval covers sampling noise "
            "only. If these variants were selected for being significant, "
            "coverage falls below nominal and no correction here repairs it -- "
            "see the measured numbers in the module docstring",
            RuntimeWarning,
            stacklevel=2,
        )

    if winners_curse:
        if design == "nested" and t < 1.0:
            raise ValueError(
                "winner's-curse correction is invalid for nested DOWN-sampling: "
                "the Brownian bridge is drift-free, so the true effect never "
                "enters and shrinking z1 would bias the prediction downward"
            )
        z1 = fiqt(z1)

    # --- MAF handling and the low-MAC guard ---------------------------------
    maf1 = maf2 = None
    if maf_col is not None:
        maf1 = out[maf_col].to_numpy(dtype=float)
        maf2 = (out[maf2_col].to_numpy(dtype=float)
                if maf2_col is not None else maf1)
        minor = np.minimum(maf1, 1.0 - maf1)
        expected_mac = 2.0 * n2 * np.minimum(maf2, 1.0 - maf2)
        del minor
    else:
        expected_mac = np.full(len(out), np.nan)

    pred = predict_interval(z1, n1, n2, design=design, level=level,
                            k_covar=k_covar, method=method,
                            maf1=maf1, maf2=maf2)

    rescalable = (np.isnan(expected_mac) | (expected_mac >= min_mac)) & np.isfinite(z1)

    out["z1"] = z1
    out["z1_source"] = source
    for key, values in pred.items():
        out[key] = np.where(rescalable, values, np.nan)
    out["expected_mac_n2"] = expected_mac
    out["rescalable"] = rescalable
    return out


def _parse_args(argv):
    p = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="See the module docstring for the method and its citations.",
    )
    p.add_argument("--sumstats", required=True, help="input TSV")
    p.add_argument("--out", required=True, help="output TSV")
    p.add_argument("--n1", type=float, required=True, help="observed sample size")
    p.add_argument("--n2", type=float, required=True, help="target sample size")
    p.add_argument("--beta-col", default="BETA")
    p.add_argument("--se-col", default="SE")
    p.add_argument("--p-col", default=None)
    p.add_argument("--neglog10p-col", default=None)
    p.add_argument("--maf-col", default=None)
    p.add_argument("--maf2-col", default=None,
                   help="target-setting allele frequency, if it differs")
    p.add_argument("--k-covar", type=int, default=0,
                   help="covariates excluding intercept and tested variant")
    p.add_argument("--method", default="r2", choices=("r2", "ncp"))
    p.add_argument("--design", default="nested", choices=DESIGNS)
    p.add_argument("--level", type=float, default=0.95)
    p.add_argument("--winners-curse", action="store_true",
                   help="FIQT-shrink z1 first; invalid for nested down-sampling")
    p.add_argument("--min-mac", type=float, default=DEFAULT_MIN_MAC)
    return p.parse_args(argv)


def main(argv=None):
    import pandas as pd

    args = _parse_args(argv)
    frame = pd.read_csv(args.sumstats, sep="\t")
    out = rescale_summary_stats(
        frame,
        n1=args.n1,
        n2=args.n2,
        beta_col=args.beta_col,
        se_col=args.se_col,
        p_col=args.p_col,
        neglog10p_col=args.neglog10p_col,
        maf_col=args.maf_col,
        maf2_col=args.maf2_col,
        k_covar=args.k_covar,
        method=args.method,
        design=args.design,
        level=args.level,
        winners_curse=args.winners_curse,
        min_mac=args.min_mac,
    )
    out.to_csv(args.out, sep="\t", index=False)
    n_drop = int((~out["rescalable"]).sum())
    print(
        f"rescaled {len(out) - n_drop}/{len(out)} variants "
        f"n={args.n1:g} -> n={args.n2:g} (t={float(information_fraction(args.n1, args.n2, args.k_covar)):.4g}, "
        f"design={args.design})",
        file=sys.stderr,
    )
    if n_drop:
        print(
            f"  {n_drop} withheld: expected minor-allele count at n2 below "
            f"{args.min_mac:g}",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
