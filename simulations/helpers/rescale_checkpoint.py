#!/usr/bin/env python
"""Re-express a SLiM checkpoint under a different Q rescaling factor.

Used at the cattle deep-history handoff. The burn-in and epochs 2-7 run at
``Q_scaling = 1`` (real sizes, 17,000 -> 1,500), which makes a 30,000-generation
burn-in affordable; epochs 8-12 then have to run at ``Q_scaling = 0.01``, where
the population is inflated 100x so that 9,000 individuals can be sampled from a
terminal real Ne of 90. Both are discretizations of the same diffusion process --
with N' = N/Q, mu' = mu*Q, r' = r*Q and s' = s*Q, the products 4*N*mu, 4*N*r and
2*N*s are all preserved -- so the switch is a change of representation, not of
model. Three pieces of state have to convert:

1. Population size. Handled by SLiM: farm_selection_from_ep8.slim already calls
   setSubpopulationSize(ep_sizes[8]/Q) on the first tick after the read, and WF
   reproduction copes with a parental generation smaller than the offspring
   generation. That single step adds one generation of drift at the pre-handoff
   size -- about 0.6% of all drift accumulated between the handoff and sampling,
   roughly a ninth of what the (now fixed) epoch-8 double-run was contributing.

2. Selection coefficients of mutations already segregating. THIS SCRIPT. The
   Q=1 phase recorded s_real, but the Q=0.01 script's mutation types generate
   new mutations at 0.01*s_real, so the carried-over ones must be scaled to
   match or they would be 100x too strongly selected.

   Done here rather than with SLiM's setSelectionCoeff() because SLiM writes
   mutation metadata at treeSeqOutput() time and it is not obvious whether
   already-recorded entries would reflect the old or the new value. Assigning
   coefficients in Python and letting SLiM read them is the path this pipeline
   already relies on: farm_create_orig_pop_e2.py does exactly that to seed the
   burn-in, and create_gwas_files_and_phenotypes.py:relabel_ag_variants uses the
   same packset_metadata round-trip.

3. Tree-sequence time units. NOT handled here, deliberately. Node times are in
   ticks, and one tick is one real generation before the handoff but 0.01 after
   it, so the output has a piecewise time scale. Rewriting node times would mean
   editing SLiM's tick bookkeeping in the metadata as well, which is more places
   to get wrong. Instead create_gwas_files_and_phenotypes.py:add_neutral() lays
   the neutral mutations down in two passes with complementary time windows --
   see its --handoff_ticks argument. If you change that, check Watterson's theta
   from the final tree sequence still lands at 8.4e-9 per bp per real generation.

Example:
    python helpers/rescale_checkpoint.py \
        --input  farm_selection_Q_1.L_2000000.seed_20250303.ep7.ts \
        --output farm_selection_Q_0.01.L_2000000.seed_20250303.ep7.ts \
        --from-Q 1 --to-Q 0.01
"""
import argparse

import numpy as np
import tskit


def rescale_selection_coeffs(ts, factor):
    """Multiply every mutation's recorded selection_coeff by `factor`."""
    tables = ts.dump_tables()
    schema = tables.mutations.metadata_schema

    new_metadata = []
    n_muts = 0
    for row in tables.mutations:
        md = row.metadata
        for ml in md["mutation_list"]:
            ml["selection_coeff"] = float(ml["selection_coeff"]) * factor
            n_muts += 1
        new_metadata.append(md)

    tables.mutations.packset_metadata(
        [schema.validate_and_encode_row(md) for md in new_metadata]
    )
    return tables.tree_sequence(), n_muts


def _selco_summary(ts):
    vals = [ml["selection_coeff"]
            for m in ts.mutations()
            for ml in m.metadata["mutation_list"]]
    if not vals:
        return "no mutations"
    a = np.abs(np.asarray(vals, dtype=float))
    nz = a[a > 0]
    if nz.size == 0:
        return f"{len(vals)} mutations, all neutral"
    return (f"{len(vals)} mutations, |s| median {np.median(nz):.4e}, "
            f"max {nz.max():.4e}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Checkpoint written by the --from-Q run")
    p.add_argument("--output", required=True, help="Checkpoint to feed the --to-Q run")
    p.add_argument("--from-Q", dest="from_q", type=float, required=True,
                   help="Q_scaling the input was produced under (deep history: 1)")
    p.add_argument("--to-Q", dest="to_q", type=float, required=True,
                   help="Q_scaling the output will be consumed under (epochs 8-12: 0.01)")
    args = p.parse_args()

    if args.from_q <= 0 or args.to_q <= 0:
        raise SystemExit("--from-Q and --to-Q must both be positive")

    ts = tskit.load(args.input)
    print(f"Loaded {args.input}")
    print(f"  individuals={ts.num_individuals} sites={ts.num_sites} "
          f"mutations={ts.num_mutations} sequence_length={ts.sequence_length:.0f}")
    print(f"  before: {_selco_summary(ts)}")

    if args.from_q == args.to_q:
        print("--from-Q == --to-Q; nothing to rescale, copying through")
        ts.dump(args.output)
        return

    factor = args.to_q / args.from_q
    out, n_muts = rescale_selection_coeffs(ts, factor)
    print(f"Scaled selection_coeff on {n_muts} mutation records by {factor:g} "
          f"(Q {args.from_q:g} -> {args.to_q:g})")
    print(f"  after:  {_selco_summary(out)}")

    out.dump(args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
