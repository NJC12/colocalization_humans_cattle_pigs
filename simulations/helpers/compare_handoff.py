#!/usr/bin/env python
"""Compare the split-Q hybrid cattle history against a pure single-Q run.

Checklist item 9. validate_handoff.sbatch produces, for each seed, a `pure` arm
(everything at Q=QS) and a `hybrid` arm (deep history at Q=1, rescaled, then
epochs 8-12 at Q=QS). If the split-Q architecture is sound the two arms are
draws from the same distribution, so the test is a paired comparison across
seeds rather than a single-run eyeball.

Statistics compared, all on the final sample of each arm:
  - folded SFS in log2 minor-allele-count bins (TVD between arms)
  - nucleotide diversity pi
  - Watterson theta
  - Tajima's D
  - LD decay: mean r^2 in distance bins

The interesting quantity is not whether the arms differ -- with finite seeds
they always will -- but whether the between-arm difference is larger than the
between-seed spread WITHIN an arm. That is what the paired t-statistic and the
"within-arm SD" column report. A |t| below ~2 with the between-arm difference
comfortably inside the within-arm spread means the handoff is not detectable.

Usage:
    python helpers/compare_handoff.py \
        --root /n/scratch/.../handoff_validation \
        --out-tsv figures_and_tables/burnin_diagnostics/handoff_validation.tsv \
        --out-plot figures_and_tables/burnin_diagnostics/handoff_validation.png
"""
import argparse
import glob
import os
import re

import numpy as np

N_SFS_BINS = 16
# Upper edges in bp. The 2 Mb production region is 4x wider, but LD structure at
# the short end is what the handoff could plausibly disturb -- one extra
# generation of drift at 1500 parents perturbs haplotype frequencies, and that
# shows up first in tight LD.
LD_EDGES = np.array([1e3, 5e3, 1e4, 5e4, 1e5, 2.5e5, 5e5])


def folded_log2_sfs(ts, n_bins=N_SFS_BINS):
    """Normalized folded SFS binned by log2(minor allele count)."""
    if ts.num_sites == 0:
        return np.zeros(n_bins)
    afs = ts.allele_frequency_spectrum(polarised=True, span_normalise=False, mode="site")
    n = len(afs) - 1
    counts = np.arange(len(afs))
    minor = np.minimum(counts, n - counts)
    keep = (minor > 0) & (afs > 0)
    if not keep.any():
        return np.zeros(n_bins)
    idx = np.minimum(np.floor(np.log2(minor[keep])).astype(int), n_bins - 1)
    binned = np.bincount(idx, weights=afs[keep], minlength=n_bins)[:n_bins]
    tot = binned.sum()
    return binned / tot if tot > 0 else binned


def tvd(p, q):
    if p.sum() == 0 or q.sum() == 0:
        return np.nan
    return 0.5 * np.abs(p - q).sum()


def ld_decay(ts, max_variants=1500, min_maf=0.05, seed=0):
    """Mean r^2 by distance bin, over a random subsample of common variants.

    Subsampled because r^2 is an all-pairs computation: 1500 variants is ~1.1M
    pairs, which is instant, while every common variant in the region would be
    tens of millions. MAF-filtered because r^2 between rare variants is
    dominated by sampling noise and would swamp the comparison.
    """
    rng = np.random.default_rng(seed)
    pos, geno = [], []
    n_hap = ts.num_samples
    for v in ts.variants():
        if v.num_alleles != 2:
            continue
        g = v.genotypes
        p = g.mean()
        if min(p, 1 - p) < min_maf:
            continue
        pos.append(v.site.position)
        geno.append(g)
    if len(pos) < 10:
        return np.full(len(LD_EDGES), np.nan), 0
    pos = np.asarray(pos, dtype=float)
    G = np.asarray(geno, dtype=np.float64)
    if len(pos) > max_variants:
        sel = np.sort(rng.choice(len(pos), max_variants, replace=False))
        pos, G = pos[sel], G[sel]
    Gc = G - G.mean(axis=1, keepdims=True)
    sd = Gc.std(axis=1)
    ok = sd > 0
    pos, Gc, sd = pos[ok], Gc[ok], sd[ok]
    R = (Gc @ Gc.T) / (len(Gc[0]) * np.outer(sd, sd))
    R2 = R ** 2
    iu = np.triu_indices(len(pos), k=1)
    d = np.abs(pos[iu[0]] - pos[iu[1]])
    r2 = R2[iu]
    out = np.full(len(LD_EDGES), np.nan)
    lo = 0.0
    for i, hi in enumerate(LD_EDGES):
        m = (d >= lo) & (d < hi)
        if m.sum() >= 20:
            out[i] = r2[m].mean()
        lo = hi
    return out, len(pos)


def apply_overlay(ts, arm, qs, seed):
    """Add neutral mutations exactly as stage 2 would for this arm.

    This is what makes the comparison production-faithful rather than a bare
    demographic check. The two arms differ in TWO ways, not one: the hybrid took
    the split-Q route through the demography, and its tree sequence therefore
    carries a piecewise time scale that stage 2 must overlay in two segments,
    while the pure arm's is uniform and takes a single rate. Comparing the arms
    only on SLiM's selected mutations would test the first difference and skip
    the second, even though a botched overlay is the failure mode with the
    larger effect on downstream results.

    Both arms should come out at 8.4e-9 per bp per real generation.
    """
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from create_gwas_files_and_phenotypes import add_neutral
    if arm == "hybrid":
        # Epochs 8-12 = 24 real generations, which is 24/qs ticks at the shallow Q.
        # Everything older ran at Q=1, hence deep_Q_scaling=1.
        return add_neutral(ts, Q=1.0 / qs, seed=seed,
                           handoff_ticks=round(24 / qs), deep_Q_scaling=1)
    return add_neutral(ts, Q=1.0 / qs, seed=seed)


def analyse(path, seed, arm=None, qs=None, overlay=False):
    import tskit
    ts = tskit.load(path)
    if overlay:
        ts = apply_overlay(ts, arm, qs, seed)
    ld, n_ld = ld_decay(ts, seed=seed)
    return {
        "n_indiv": ts.num_individuals,
        "n_sites": ts.num_sites,
        "pi": float(ts.diversity(mode="site")),
        "theta_w": float(ts.segregating_sites(mode="site")),
        "tajimas_d": float(ts.Tajimas_D(mode="site")) if ts.num_sites > 1 else np.nan,
        "sfs": folded_log2_sfs(ts),
        "ld": ld,
        "n_ld_vars": n_ld,
    }


def paired_report(name, pure, hyb, fh):
    """Paired across seeds: is the arm difference bigger than the seed spread?"""
    p, h = np.asarray(pure, float), np.asarray(hyb, float)
    ok = np.isfinite(p) & np.isfinite(h)
    p, h = p[ok], h[ok]
    if len(p) < 2:
        line = f"  {name:14s} insufficient finite pairs ({len(p)})"
        print(line); fh.write(line + "\n"); return
    d = h - p
    sd = d.std(ddof=1)
    se = sd / np.sqrt(len(d))
    t = d.mean() / se if se > 0 else np.nan
    within = np.std(np.concatenate([p, h]), ddof=1)
    rel = 100 * d.mean() / p.mean() if p.mean() != 0 else np.nan
    line = (f"  {name:14s} pure={p.mean():>11.5g}  hybrid={h.mean():>11.5g}  "
            f"diff={d.mean():>+11.4g} ({rel:>+6.1f}%)  t={t:>6.2f}  "
            f"within-arm SD={within:.4g}")
    print(line); fh.write(line + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, help="handoff_validation dir holding pure/ and hybrid/")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-plot", default=None)
    ap.add_argument("--overlay", action="store_true",
                    help="Add neutral mutations the way stage 2 would (uniform rate for "
                         "the pure arm, two-segment for the hybrid's piecewise time scale) "
                         "before computing statistics. Recommended: without it the "
                         "comparison tests the demographic handoff but not the overlay "
                         "correction that accompanies it.")
    ap.add_argument("--QS", type=float, default=0.1,
                    help="Shallow Q the arms ran at; must match submit_handoff_validation.sh")
    args = ap.parse_args()

    arms = {}
    for arm in ("pure", "hybrid"):
        for d in sorted(glob.glob(os.path.join(args.root, arm, "seed_*"))):
            seed = int(re.search(r"seed_(\d+)$", d).group(1))
            ts = os.path.join(d, "sel.full.ts")
            if not os.path.exists(ts):
                print(f"  (skip: no sel.full.ts in {d})")
                continue
            arms.setdefault(arm, {})[seed] = analyse(
                ts, seed, arm=arm, qs=args.QS, overlay=args.overlay)
            r = arms[arm][seed]
            print(f"  {arm:6s} seed {seed}: indiv={r['n_indiv']} sites={r['n_sites']} "
                  f"pi={r['pi']:.4e} D={r['tajimas_d']:.3f} ld_vars={r['n_ld_vars']}")

    seeds = sorted(set(arms.get("pure", {})) & set(arms.get("hybrid", {})))
    if not seeds:
        raise SystemExit("no seed has BOTH arms complete; nothing to compare")
    print(f"\n{len(seeds)} paired seeds: {seeds}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)), exist_ok=True)
    with open(args.out_tsv, "w") as fh:
        fh.write("# paired handoff validation: hybrid (split-Q) vs pure (single-Q)\n")
        fh.write(f"# seeds: {seeds}\n\n")
        fh.write("arm\tseed\tn_indiv\tn_sites\tpi\ttheta_w\ttajimas_d\n")
        for arm in ("pure", "hybrid"):
            for s in seeds:
                r = arms[arm][s]
                fh.write(f"{arm}\t{s}\t{r['n_indiv']}\t{r['n_sites']}\t{r['pi']:.6e}\t"
                         f"{r['theta_w']:.6e}\t{r['tajimas_d']:.4f}\n")

        fh.write("\n# paired comparison (hybrid - pure)\n")
        print("\n=== paired comparison (hybrid - pure) ===")
        for key in ("n_sites", "pi", "theta_w", "tajimas_d"):
            paired_report(key, [arms["pure"][s][key] for s in seeds],
                          [arms["hybrid"][s][key] for s in seeds], fh)

        print("\n=== LD decay, mean r^2 by distance bin ===")
        fh.write("\n# LD decay, mean r^2 by distance bin\n")
        lo = 0
        for i, hi in enumerate(LD_EDGES):
            paired_report(f"r2 {int(lo/1000)}-{int(hi/1000)}kb",
                          [arms["pure"][s]["ld"][i] for s in seeds],
                          [arms["hybrid"][s]["ld"][i] for s in seeds], fh)
            lo = hi

        print("\n=== folded SFS: TVD between arms, vs TVD between seeds within an arm ===")
        fh.write("\n# folded SFS TVD\n")
        cross = [tvd(arms["pure"][s]["sfs"], arms["hybrid"][s]["sfs"]) for s in seeds]
        within = []
        for arm in ("pure", "hybrid"):
            for i in range(len(seeds)):
                for j in range(i + 1, len(seeds)):
                    within.append(tvd(arms[arm][seeds[i]]["sfs"], arms[arm][seeds[j]]["sfs"]))
        for lbl, v in (("cross-arm (paired)", cross), ("within-arm (seed pairs)", within)):
            line = f"  {lbl:26s} mean={np.nanmean(v):.4f}  max={np.nanmax(v):.4f}  n={len(v)}"
            print(line); fh.write(line + "\n")
        verdict = ("HANDOFF NOT DETECTABLE: cross-arm SFS distance is within the "
                   "seed-to-seed noise of a single arm."
                   if np.nanmean(cross) <= np.nanmean(within) else
                   "CROSS-ARM SFS DISTANCE EXCEEDS WITHIN-ARM NOISE -- inspect before trusting the handoff.")
        print(f"\n  => {verdict}")
        fh.write(f"\n# {verdict}\n")
    print(f"\nWrote {args.out_tsv}")

    if args.out_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        x = np.arange(N_SFS_BINS)
        for arm, c in (("pure", "C0"), ("hybrid", "C1")):
            m = np.mean([arms[arm][s]["sfs"] for s in seeds], axis=0)
            sd = np.std([arms[arm][s]["sfs"] for s in seeds], axis=0)
            ax[0].plot(x, m, marker="o", ms=3, color=c, label=arm)
            ax[0].fill_between(x, m - sd, m + sd, color=c, alpha=0.2)
        ax[0].set_xlabel("log2(minor allele count) bin"); ax[0].set_ylabel("proportion")
        ax[0].set_title("Folded SFS (mean +/- SD over seeds)"); ax[0].legend(frameon=False)

        centers = np.concatenate([[LD_EDGES[0] / 2], (LD_EDGES[:-1] + LD_EDGES[1:]) / 2])
        for arm, c in (("pure", "C0"), ("hybrid", "C1")):
            m = np.nanmean([arms[arm][s]["ld"] for s in seeds], axis=0)
            sd = np.nanstd([arms[arm][s]["ld"] for s in seeds], axis=0)
            ax[1].plot(centers, m, marker="o", ms=3, color=c, label=arm)
            ax[1].fill_between(centers, m - sd, m + sd, color=c, alpha=0.2)
        ax[1].set_xscale("log"); ax[1].set_xlabel("distance (bp)"); ax[1].set_ylabel("mean r^2")
        ax[1].set_title("LD decay"); ax[1].legend(frameon=False)
        fig.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(args.out_plot)), exist_ok=True)
        fig.savefig(args.out_plot, dpi=150)
        print(f"Wrote {args.out_plot}")


if __name__ == "__main__":
    main()
