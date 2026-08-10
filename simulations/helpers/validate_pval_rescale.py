#!/usr/bin/env python
"""Measure whether helpers/pval_rescale.py is actually calibrated.

The module asserts a distribution: for a NESTED subsample at information
fraction t, z2 | z1 ~ N(z1 sqrt(t), 1 - min(t, 1/t)).  This script tests that
against real least-squares refits rather than taking it on faith.

Design
------
Simulate a cohort with the features that could break the Brownian-bridge
argument -- LD between variants, covariates that get re-fitted in every
subsample, a re-estimated residual variance, and a MAF spectrum reaching down to
where the normal approximation fails:

    haplotypes   first-order autoregressive latent Gaussian thresholded at each
                 variant's frequency, giving exponential LD decay
    genotypes    sum of two independent haplotypes (HWE)
    phenotype    a few causal variants + covariate effects + Gaussian noise
    GWAS         per-variant OLS of y on 1 + g + covariates, computed via
                 Frisch-Waugh-Lovell so it is exact, not approximate

Then, for each target size, draw many nested subsamples, re-run the GWAS, and
compare the observed statistic against the predicted distribution.

Reported
--------
  slope, intercept   regression of observed z2 on predicted centre.  Target 1, 0
                     when DOWN-sampling.  When UP-sampling the predictor is
                     itself noisy (z1 carries SD 1 into a centre inflated by
                     sqrt(t)), so regression dilution attenuates the slope
                     towards ~1/t -- read `coverage`, not `slope`, on that leg.
  resid_sd           SD of (observed - predicted centre).  Target sqrt(1 - t)
                     for down-sampling; the `theory_sd` column is that target.
  coverage           fraction inside the nominal 95% interval.  Target 0.95.
  low-MAC rows       the same numbers restricted to variants the MAC guard
                     withholds, to confirm the guard sits where the model breaks.

The up-sampling leg additionally selects variants on significance in a small
subsample and predicts the full cohort, with and without the FIQT winner's-curse
correction -- the asymmetry claimed in the module docstring (correction helps
going up, is invalid going down).

Usage
-----
    python helpers/validate_pval_rescale.py --out calibration.tsv
"""

import argparse
import sys

import numpy as np
from scipy import stats

try:
    from helpers import pval_rescale as pr
except ImportError:  # running the file directly from helpers/
    import os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from helpers import pval_rescale as pr


def simulate_cohort(n, m, k_covar, rng, ld_rho=0.85, n_causal=12, h2=0.15,
                    maf_min=0.002, maf_max=0.5):
    """Genotypes with LD, covariates, and a quantitative phenotype.

    Returns (G, y, C, maf) with G float32 (n x m), C the n x (k_covar + 1)
    design matrix including the intercept.
    """
    maf = rng.uniform(maf_min, maf_max, size=m)
    thresh = stats.norm.isf(maf)  # latent > thresh  =>  carries the minor allele

    n_hap = 2 * n
    G = np.empty((n, m), dtype=np.float32)
    x = rng.standard_normal(n_hap)
    root = np.sqrt(1.0 - ld_rho**2)
    for j in range(m):
        if j:
            x = ld_rho * x + root * rng.standard_normal(n_hap)
        hap = (x > thresh[j]).astype(np.float32)
        G[:, j] = hap[:n] + hap[n:]

    C = np.column_stack([np.ones(n), rng.standard_normal((n, k_covar))])

    causal = rng.choice(m, size=n_causal, replace=False)
    # equal per-variant contribution on the standardized scale
    gvar = 2.0 * maf[causal] * (1.0 - maf[causal])
    beta = np.sqrt(h2 / n_causal / gvar) * rng.choice([-1.0, 1.0], size=n_causal)
    g = (h2 > 0) * (G[:, causal].astype(np.float64) @ beta)

    cov_effect = C[:, 1:] @ rng.normal(0.3, 0.1, size=k_covar)
    y = g + cov_effect + rng.standard_normal(n) * np.sqrt(max(1.0 - h2, 1e-6))
    return G, y, C, maf, causal


def gwas(G, y, C, rows=None, chunk=512):
    """Per-variant OLS t-statistic for y ~ 1 + g + C, via Frisch-Waugh-Lovell.

    Residualizing y and each g on C gives exactly the full-model coefficient and
    residual sum of squares, so this is the same number plink2 would report, not
    an approximation.
    """
    if rows is not None:
        G, y, C = G[rows], y[rows], C[rows]
    n, m = G.shape
    k = C.shape[1] - 1                      # covariates, excluding intercept
    df = n - k - 2                          # intercept + covariates + the variant
    if df <= 0:
        raise ValueError(f"n={n} too small for k_covar={k}")

    q, _ = np.linalg.qr(C)                  # orthonormal basis for the covariates
    yr = y - q @ (q.T @ y)
    yy = float(yr @ yr)

    tstat = np.empty(m)
    beta_out = np.empty(m)
    mac = np.empty(m)
    for s in range(0, m, chunk):
        block = G[:, s:s + chunk].astype(np.float64)
        mac[s:s + chunk] = block.sum(axis=0)
        block -= q @ (q.T @ block)
        gg = np.einsum("ij,ij->j", block, block)
        gy = block.T @ yr
        with np.errstate(divide="ignore", invalid="ignore"):
            beta = gy / gg
            rss = yy - beta * gy
            se = np.sqrt(rss / df / gg)
            tstat[s:s + chunk] = beta / se
        beta_out[s:s + chunk] = beta
    # variants monomorphic in this subsample carry no information
    tstat[~np.isfinite(tstat)] = np.nan
    # report the MINOR allele count actually observed
    mac = np.minimum(mac, 2 * n - mac)
    return tstat, beta_out, mac


def _fit_line(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 10:
        return np.nan, np.nan
    slope, intercept = np.polyfit(x[ok], y[ok], 1)
    return slope, intercept


def down_sampling_leg(G, y, C, z1, maf, n1, targets, n_draws, k_covar, rng,
                      min_mac, level=0.95):
    rows_out = []
    for n2 in targets:
        pred = pr.predict_interval(z1, n1, n2, design="nested", level=level,
                                   k_covar=k_covar)
        t = float(pr.information_fraction(n1, n2, k_covar))
        theory_sd = np.sqrt(1.0 - min(t, 1.0 / t))
        expected_mac = 2.0 * n2 * np.minimum(maf, 1.0 - maf)
        keep = expected_mac >= min_mac

        obs, ctr, inside, mask_keep = [], [], [], []
        for _ in range(n_draws):
            idx = rng.choice(len(y), size=n2, replace=False)
            z2, _, _ = gwas(G, y, C, rows=idx)
            obs.append(z2)
            ctr.append(pred["z2"])
            inside.append((z2 >= pred["z2_lo"]) & (z2 <= pred["z2_hi"]))
            mask_keep.append(keep)

        obs = np.concatenate(obs)
        ctr = np.concatenate(ctr)
        inside = np.concatenate(inside)
        mask_keep = np.concatenate(mask_keep)

        for label, sel in (("pass_mac_guard", mask_keep),
                           ("withheld_low_mac", ~mask_keep)):
            if sel.sum() == 0:
                continue
            o, c, i = obs[sel], ctr[sel], inside[sel]
            ok = np.isfinite(o)
            slope, intercept = _fit_line(c, o)
            rows_out.append(dict(
                leg="down", n1=n1, n2=n2, t=round(t, 5), subset=label,
                n_variants=int(sel.sum() / n_draws), n_obs=int(ok.sum()),
                slope=slope, intercept=intercept,
                resid_sd=float(np.nanstd(o - c)), theory_sd=float(theory_sd),
                coverage=float(np.mean(i[ok])),
            ))
    return rows_out


def up_sampling_leg(G, y, C, z_full, n1, n_small, n_draws, k_covar, rng,
                    select_p=1e-4, level=0.95):
    """Discover in a nested subsample, predict the full cohort.

    Here the drift does NOT cancel, the plug-in uses z_small as the effect
    estimate, and selecting on significance biases it upward -- so this is where
    the winner's-curse correction belongs.
    """
    z_crit = stats.norm.isf(select_p / 2)
    t = float(pr.information_fraction(n_small, n1, k_covar))
    theory_sd = float(np.sqrt(abs(t - 1.0)))
    rows_out = []

    # ---- unselected: tests the variance formula without winner's curse ------
    resid, inside_all, ctr_all, obs_all = [], [], [], []
    draws = []
    for _ in range(n_draws):
        idx = rng.choice(len(y), size=n_small, replace=False)
        z_small, _, _ = gwas(G, y, C, rows=idx)
        draws.append(z_small)
        finite = np.isfinite(z_small)
        pred = pr.predict_interval(z_small[finite], n_small, n1,
                                   design="nested", level=level, k_covar=k_covar)
        obs_all.append(z_full[finite])
        ctr_all.append(pred["z2"])
        resid.append(z_full[finite] - pred["z2"])
        inside_all.append((z_full[finite] >= pred["z2_lo"])
                          & (z_full[finite] <= pred["z2_hi"]))
    obs_all = np.concatenate(obs_all)
    ctr_all = np.concatenate(ctr_all)
    slope, intercept = _fit_line(ctr_all, obs_all)
    rows_out.append(dict(
        leg="up", n1=n_small, n2=n1, t=round(t, 5), subset="all_variants",
        n_variants=obs_all.size // max(n_draws, 1), n_obs=int(obs_all.size),
        slope=slope, intercept=intercept,
        resid_sd=float(np.nanstd(np.concatenate(resid))), theory_sd=theory_sd,
        coverage=float(np.mean(np.concatenate(inside_all))),
        mean_overprediction=float(np.nanmean(-np.concatenate(resid))),
    ))

    # ---- selected on significance: this is where winner's curse bites -------
    for correction in ("none", "fiqt"):
        bias, inside, n_sel = [], [], 0
        for z_small in draws:
            finite = np.isfinite(z_small)
            lam = z_small.copy()
            if correction == "fiqt":
                lam[finite] = pr.fiqt(z_small[finite])
            sel = finite & (np.abs(z_small) > z_crit)
            if sel.sum() == 0:
                continue
            n_sel += int(sel.sum())
            pred = pr.predict_interval(lam[sel], n_small, n1, design="nested",
                                       level=level, k_covar=k_covar)
            bias.append(np.abs(pred["z2"]) - np.abs(z_full[sel]))
            inside.append((z_full[sel] >= pred["z2_lo"])
                          & (z_full[sel] <= pred["z2_hi"]))
        if not bias:
            continue
        bias = np.concatenate(bias)
        rows_out.append(dict(
            leg="up", n1=n_small, n2=n1, t=round(t, 5),
            subset=f"selected_p<{select_p:g}_wc={correction}",
            n_variants=n_sel // max(n_draws, 1), n_obs=int(bias.size),
            slope=np.nan, intercept=np.nan,
            resid_sd=float(np.nanstd(bias)), theory_sd=theory_sd,
            coverage=float(np.mean(np.concatenate(inside))),
            mean_overprediction=float(np.nanmean(bias)),
        ))
    return rows_out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n1", type=int, default=8000, help="full cohort size")
    p.add_argument("--m", type=int, default=3000, help="variants")
    p.add_argument("--k-covar", type=int, default=5)
    p.add_argument("--targets", type=int, nargs="+",
                   default=[4000, 2000, 1000, 500])
    p.add_argument("--n-draws", type=int, default=40)
    p.add_argument("--n-small", type=int, default=2000,
                   help="discovery size for the up-sampling leg")
    p.add_argument("--min-mac", type=float, default=pr.DEFAULT_MIN_MAC)
    p.add_argument("--seed", type=int, default=20260731)
    p.add_argument("--out", default=None, help="write the table as TSV")
    args = p.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    print(f"simulating {args.n1} individuals x {args.m} variants "
          f"({args.k_covar} covariates, LD rho=0.85) ...", file=sys.stderr)
    G, y, C, maf, causal = simulate_cohort(args.n1, args.m, args.k_covar, rng)

    z1, _, mac_full = gwas(G, y, C)
    print(f"full-cohort GWAS: {np.sum(np.abs(z1) > 5.45)} variants at p < 5e-8 "
          f"({len(causal)} causal simulated)", file=sys.stderr)

    rows = down_sampling_leg(G, y, C, z1, maf, args.n1, args.targets,
                             args.n_draws, args.k_covar, rng, args.min_mac)
    rows += up_sampling_leg(G, y, C, z1, args.n1, args.n_small,
                            args.n_draws, args.k_covar, rng)

    header = ["leg", "n1", "n2", "t", "subset", "n_variants", "n_obs",
              "slope", "intercept", "resid_sd", "theory_sd", "coverage",
              "mean_overprediction"]
    lines = ["\t".join(header)]
    for r in rows:
        lines.append("\t".join(
            "" if r.get(h) is None or (isinstance(r.get(h), float) and np.isnan(r.get(h)))
            else (f"{r[h]:.4g}" if isinstance(r.get(h), float) else str(r.get(h, "")))
            for h in header
        ))
    table = "\n".join(lines)
    print(table)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(table + "\n")
        print(f"\nwrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
