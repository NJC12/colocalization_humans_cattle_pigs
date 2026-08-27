"""Rebuild a missing `*_trait_partners_*.tsv` from the archived variant tables.

WHY A RECONSTRUCTION IS POSSIBLE AT ALL

The partner table records, per GWAS causal locus, whether a GTEx trait was
defined at the SAME causal variant. Every input to that question survives in the
`*_vars_*.tsv` files: a trait is named `tr<position>`, and a variant is causal in
a panel exactly when its `beta` is non-zero there. So the table is a function of
two files that are already archived, and no simulation has to be re-run --
which matters, because re-running stage 2 would redraw the trait set and
invalidate the stage-3/4/5 output already computed on top of it.

WHY IT IS NEEDED

`create_gwas_files_and_phenotypes.py` used to gate the table on
`causal_sampling: power`, back when power sampling was the only scheme that
topped up the GTEx causal set instead of intersecting it. The drawn-DFE arms
(H, I) top up under uniform sampling too, so runs made before the gate was
widened have no table despite needing one.

SCOPE

Only for arms whose GTEx set is topped up -- H and I, and power-sampled runs.
Where the set is INTERSECTED (A/E under uniform), every GWAS locus has a partner
by construction, stage 2 deliberately writes nothing, and reconstructing a
file of all-True rows would invent a distinction that does not exist. This
refuses to write for those unless --force is given.

NOT handled: the redistribution arms B and D. There the trait is named for the
RECIPIENT of a moved effect, and while the vars file's beta does sit on the
recipient, verifying that equivalence is a separate question from the one this
script answers. It refuses on those categories.

    python helpers/reconstruct_trait_partners.py <replicate_dir> [...] [--apply]
"""

import argparse
import glob
import os
import re
import sys

TOPUP_CATEGORIES = set("HI")          # drawn-DFE arms; power runs are detected separately
REDISTRIBUTION_CATEGORIES = set("BD")  # refuse: trait id is the recipient, not the donor


def causal_positions(path):
    """{int position} where beta != 0, i.e. the loci causal in that panel."""
    out = set()
    with open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        pi, bi = cols.index("position"), cols.index("beta")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                if float(f[bi]) != 0.0:
                    out.add(int(float(f[pi])))
            except (ValueError, IndexError):
                continue
    return out


def rebuild(rep_dir, apply_=False, force=False):
    rep = os.path.basename(rep_dir.rstrip("/"))
    cat = rep[0] if rep else "?"
    if cat in REDISTRIBUTION_CATEGORIES:
        return f"{rep}: SKIP -- category {cat} names traits for the effect RECIPIENT"

    made = []
    for gwas_vars in sorted(glob.glob(os.path.join(rep_dir, "*gwas_vars_*.tsv"))):
        prefix = os.path.basename(gwas_vars).split("_vars_")[0]      # hgwas / cgwas
        suffix = os.path.basename(gwas_vars).split("_vars_")[1]      # gwas_5_gtex_20_maf_0.01.tsv
        out = os.path.join(rep_dir, f"{prefix}_trait_partners_{suffix}")
        if os.path.exists(out):
            return f"{rep}: already has {os.path.basename(out)}"

        # The 1000-individual GTEx panel is the one the causal set is defined
        # against; the 500/250 subsamples are drawn FROM it afterwards.
        species = prefix[0]
        gtex_vars = os.path.join(rep_dir, f"{species}gtex_vars_{suffix}")
        if not os.path.exists(gtex_vars):
            return f"{rep}: no {os.path.basename(gtex_vars)} to compare against"

        if cat not in TOPUP_CATEGORIES and not force:
            return (f"{rep}: SKIP -- category {cat} intersects its GTEx causal set, so "
                    f"every locus has a partner by construction and stage 2 writes "
                    f"no table. Pass --force to write one anyway.")

        gwas_pos = sorted(causal_positions(gwas_vars))
        gtex_pos = causal_positions(gtex_vars)
        rows = [(f"tr{p}", f"tr{p}" if p in gtex_pos else "", p, p in gtex_pos)
                for p in gwas_pos]
        n_shared = sum(1 for r in rows if r[3])
        if apply_:
            with open(out, "w") as fh:
                fh.write("gwas_trait\tgtex_trait\tposition\tshared\n")
                for t, g, p, s in rows:
                    fh.write(f"{t}\t{g}\t{p}\t{s}\n")
        made.append(f"{os.path.basename(out)} ({n_shared}/{len(rows)} shared)")
    if not made:
        return f"{rep}: no *gwas_vars_* table found"
    return f"{rep}: {'wrote' if apply_ else 'would write'} " + ", ".join(made)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+", help="replicate directories")
    ap.add_argument("--apply", action="store_true", help="write; otherwise report only")
    ap.add_argument("--force", action="store_true",
                    help="also write for categories whose GTEx set is intersected")
    args = ap.parse_args()
    for d in args.dirs:
        print(rebuild(d, args.apply, args.force))
    if not args.apply:
        print("\nDry run -- nothing written. Re-run with --apply.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
