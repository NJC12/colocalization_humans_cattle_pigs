#!/usr/bin/env python
"""Extract the sub-interval [0, end) from a SLiM tree sequence.

Used to derive a smaller-L burn-in from a larger one without re-running the
burn-in -- e.g. the first 2 Mb of the equilibrated 10 Mb cattle burn-in
(`cycle_25000.ts`). The genealogy/sites outside [0, end) are dropped and the
coordinate system is truncated to `end`, so the output is a valid L=end tree.

ALL individuals are retained (only the genome is truncated), so the population
size is unchanged -- `farm_selection.slim`'s size check (== 17000/Q) still
passes on the result.

Notes
-----
* `simplify=False`: we keep every node/edge overlapping [0, end) without
  renumbering, to perturb the SLiM-specific table structure as little as
  possible. SLiM 5 re-reads the result via `readFromPopulationFile`; the
  caller must `initialize` the genome to the same length (L = end).
* Caveat (scientific): the extracted region carries the linked/background-
  selection signature of having been embedded in the larger region, which
  differs from a standalone burn-in of length `end`.
"""
import argparse
import tskit


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Input SLiM .ts (e.g. the 10 Mb cycle_25000.ts)")
    p.add_argument("--output", required=True, help="Output .ts truncated to [0, end)")
    p.add_argument("--end", type=float, required=True, help="Sub-region end in bp (e.g. 2000000)")
    args = p.parse_args()

    ts = tskit.load(args.input)
    print(f"Loaded {args.input}: sequence_length={ts.sequence_length:.0f} "
          f"individuals={ts.num_individuals} nodes={ts.num_nodes} "
          f"sites={ts.num_sites} mutations={ts.num_mutations}")

    end = args.end
    if end > ts.sequence_length:
        raise SystemExit(f"--end {end:.0f} exceeds sequence_length {ts.sequence_length:.0f}")

    # Drop genealogy/sites outside [0, end); keep all individuals/nodes.
    ts = ts.keep_intervals([(0.0, end)], simplify=False)

    # Truncate the coordinate system to [0, end). Valid because keep_intervals
    # left no edges/sites beyond `end`.
    tables = ts.dump_tables()
    tables.sequence_length = end
    ts_out = tables.tree_sequence()

    ts_out.dump(args.output)
    print(f"Wrote {args.output}: sequence_length={ts_out.sequence_length:.0f} "
          f"individuals={ts_out.num_individuals} nodes={ts_out.num_nodes} "
          f"sites={ts_out.num_sites} mutations={ts_out.num_mutations}")


if __name__ == "__main__":
    main()
