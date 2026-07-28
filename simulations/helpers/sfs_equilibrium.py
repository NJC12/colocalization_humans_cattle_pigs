#!/usr/bin/env python
"""Per-|s|-bin site-frequency-spectrum convergence check for the cattle burn-in.

The burn-in starts from an msprime coalescent whose mutations are at NEUTRAL
equilibrium frequencies; farm_create_orig_pop_e2.py then paints selection
coefficients onto them after the fact. The forward burn-in has to run long
enough for each selection class to relax from that neutral starting SFS to its
mutation-selection-drift equilibrium, and the classes do that on very different
timescales -- roughly 1/|s| generations, so ~33 generations for the strong tail
of m2 but ~2,500 for the m1 bulk (97.5% of the DFE, mean |s| = 4e-4).

That is why this reports TVD per |s| bin rather than pooled. A pooled SFS is
dominated by whichever class is most numerous and can look stationary while the
weakly-selected bulk -- the class that actually matters here, because it supplies
most of the trait-associated variants -- is still drifting. The round-2 analysis
(notebook `cattle-burnin-sfs-equilibrium`) used the pooled SFS and concluded
equilibrium by cycle 16000-18000; this exists to check that per class.

Everything hot is a C-level tskit call. Do NOT reach for ts.genotype_matrix() or
per-variant .genotypes here: the Q=0.01 burn-in has 1.7M individuals, i.e. 3.4M
sample nodes, so materializing genotypes for ~46k sites is ~10^11 entries. The
approach below instead partitions sites by |s| with delete_sites() (cheap, O(n))
and calls allele_frequency_spectrum() on each subset.

Usage:
    python helpers/sfs_equilibrium.py \
        --checkpoints '/path/to/stage1/farm_burn_in.Q_0.01.*.cycle_*.ts' \
        --Q_scaling 0.01 \
        --out-tsv  figures_and_tables/burnin_diagnostics/burnin_sfs_by_s.tsv \
        --out-plot figures_and_tables/burnin_diagnostics/burnin_sfs_by_s.png

Reading the output: within each |s| bin, TVD-vs-final should fall and then
plateau at a noise floor, and TVD-vs-previous should settle to that same floor.
A bin whose TVD is still trending down at the last checkpoint has not converged.
"""
import argparse
import glob
import os
import re

import numpy as np

# |s| bin edges in REAL (Q-unscaled) units. The first bin is everything weaker
# than 1e-4, which is where most of the m1 exponential mass sits and where
# equilibration is slowest.
S_BINS = [(0.0, 1e-4), (1e-4, 1e-3), (1e-3, 1e-2), (1e-2, np.inf)]


def bin_label(lo, hi):
    hi_s = "inf" if not np.isfinite(hi) else f"{hi:g}"
    return f"|s| {lo:g}-{hi_s}"


def cycle_of(path):
    """Tick number from a `....cycle_<N>.ts` filename."""
    m = re.search(r"\.cycle_(\d+)\.ts$", os.path.basename(path))
    if not m:
        raise ValueError(f"cannot parse a cycle number from {path!r}")
    return int(m.group(1))


def site_abs_s(ts, q_scaling):
    """|selection coefficient| in real units, per SITE (indexed by site id).

    Sites with no mutation record get NaN and are dropped. Where a site carries
    several stacked mutations the first is used, matching get_vars_df.
    """
    out = np.full(ts.num_sites, np.nan)
    for mut in ts.mutations():
        if np.isnan(out[mut.site]):
            ml = mut.metadata["mutation_list"]
            if ml:
                out[mut.site] = abs(float(ml[0]["selection_coeff"])) / q_scaling
    return out


def folded_log2_sfs(ts, n_bins=16):
    """Normalized folded SFS binned by log2(minor allele count).

    allele_frequency_spectrum(polarised=True, span_normalise=False, mode="site")
    returns counts indexed by derived allele count; folding and log2-binning
    makes the vector comparable across checkpoints without depending on sample
    size. Returns a zero vector when the subset has no segregating sites.
    """
    if ts.num_sites == 0:
        return np.zeros(n_bins)
    afs = ts.allele_frequency_spectrum(polarised=True, span_normalise=False, mode="site")
    n = len(afs) - 1                       # number of sample nodes
    counts = np.arange(len(afs))
    minor = np.minimum(counts, n - counts)  # fold
    keep = (minor > 0) & (afs > 0)
    if not keep.any():
        return np.zeros(n_bins)
    idx = np.minimum(np.floor(np.log2(minor[keep])).astype(int), n_bins - 1)
    binned = np.bincount(idx, weights=afs[keep], minlength=n_bins)[:n_bins]
    total = binned.sum()
    return binned / total if total > 0 else binned


def tvd(p, q):
    """Total variation distance between two normalized SFS vectors."""
    if p.sum() == 0 or q.sum() == 0:
        return np.nan
    return 0.5 * np.abs(p - q).sum()


def analyse(path, q_scaling):
    """Per-|s|-bin SFS and scalar summaries for one checkpoint."""
    import tskit
    ts = tskit.load(path)
    s = site_abs_s(ts, q_scaling)
    all_ids = np.arange(ts.num_sites)

    rows = []
    for lo, hi in S_BINS:
        in_bin = np.isfinite(s) & (s >= lo) & (s < hi)
        # delete_sites() keeps the genealogy and drops only the sites we are not
        # binning on, so every statistic below stays a C-level call.
        sub = ts.delete_sites(all_ids[~in_bin])
        rows.append({
            "s_bin": bin_label(lo, hi),
            "n_sites": int(in_bin.sum()),
            "sfs": folded_log2_sfs(sub),
            "pi": float(sub.diversity(mode="site")) if sub.num_sites else np.nan,
            "theta_w": float(sub.segregating_sites(mode="site")) if sub.num_sites else np.nan,
            "tajimas_d": float(sub.Tajimas_D(mode="site")) if sub.num_sites > 1 else np.nan,
        })
    return rows, ts.num_sites


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--checkpoints", required=True,
                   help="Glob for the burn-in .cycle_<N>.ts checkpoints (quote it)")
    p.add_argument("--Q_scaling", type=float, required=True,
                   help="Q_scaling the burn-in ran under; selection coefficients are "
                        "divided by it to get real units (cattle round 2: 0.01, "
                        "round-3 deep history: 1)")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-plot", default=None)
    args = p.parse_args()

    paths = sorted(glob.glob(args.checkpoints), key=cycle_of)
    if not paths:
        raise SystemExit(f"no checkpoints matched {args.checkpoints!r}")
    print(f"{len(paths)} checkpoints, cycles {cycle_of(paths[0])}..{cycle_of(paths[-1])}")

    per_cycle = {}
    for path in paths:
        rows, n_sites = analyse(path, args.Q_scaling)
        per_cycle[cycle_of(path)] = rows
        print(f"  cycle {cycle_of(path):>7}: {n_sites} sites "
              f"({', '.join(str(r['n_sites']) for r in rows)} by |s| bin)")

    cycles = sorted(per_cycle)
    final = per_cycle[cycles[-1]]

    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)), exist_ok=True)
    with open(args.out_tsv, "w") as fh:
        fh.write("cycle\ts_bin\tn_sites\tpi\ttheta_w\ttajimas_d\ttvd_vs_final\ttvd_vs_prev\n")
        for i, c in enumerate(cycles):
            for b, row in enumerate(per_cycle[c]):
                d_final = tvd(row["sfs"], final[b]["sfs"])
                if i == 0:
                    d_prev = ""                      # no earlier checkpoint to compare against
                else:
                    d_prev = f"{tvd(row['sfs'], per_cycle[cycles[i - 1]][b]['sfs']):.5f}"
                fh.write(
                    f"{c}\t{row['s_bin']}\t{row['n_sites']}\t{row['pi']:.6e}\t"
                    f"{row['theta_w']:.6e}\t{row['tajimas_d']:.4f}\t"
                    f"{d_final:.5f}\t{d_prev}\n"
                )
    print(f"Wrote {args.out_tsv}")

    if args.out_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for b, (lo, hi) in enumerate(S_BINS):
            y = [tvd(per_cycle[c][b]["sfs"], final[b]["sfs"]) for c in cycles]
            ax.plot(cycles, y, marker="o", ms=3, label=bin_label(lo, hi))
        ax.set_xlabel("burn-in tick")
        ax.set_ylabel("TVD of folded SFS vs final checkpoint")
        ax.set_title("Burn-in convergence by selection class")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(args.out_plot)), exist_ok=True)
        fig.savefig(args.out_plot, dpi=150)
        print(f"Wrote {args.out_plot}")


if __name__ == "__main__":
    main()
