"""Compare two tree sequences by GENETIC CONTENT, ignoring provenance.

WHY NOT `cmp`

A .trees file embeds a provenance record: an ISO timestamp, the software
versions, and the command line. Two runs of identical code from identical seeds
therefore differ in bytes while being the same simulation. Observed while
verifying the publication run's stage 1: a freshly re-run E1 came out 24 bytes
larger than the copy on disk, entirely in that record.

A byte comparison would have said "the code changed, re-run all 20" -- costly,
and wrong. Worse, in the other direction it teaches you to shrug at differences,
which is exactly the habit that lets a real one through.

WHAT COUNTS AS THE SAME

Everything a downstream stage can observe: the node, edge, site, mutation,
individual, population and migration tables, and the sequence length. tskit's
`TableCollection.equals(..., ignore_provenance=True)` is precisely that test.
`ignore_timestamps` is NOT enough on its own -- the record also carries versions
and the resolved command line.

When they differ, the tables are compared one at a time so the report says WHICH
table moved. A difference in `sites`/`mutations` alone means the overlay changed;
a difference in `edges`/`nodes` means the genealogy did, which is the serious one.

    python helpers/compare_tree_sequences.py A.ts B.ts [C.ts D.ts ...]

Exit 0 only if every pair matches ignoring provenance.
"""

import sys

import numpy as np
import tskit


TABLES = ("individuals", "nodes", "edges", "sites", "mutations",
          "populations", "migrations")


def compare(path_a, path_b):
    """-> (same_ignoring_provenance, detail)"""
    a = tskit.load(path_a)
    b = tskit.load(path_b)

    if a.sequence_length != b.sequence_length:
        return False, (f"sequence_length {a.sequence_length} != {b.sequence_length}")

    ta, tb = a.tables, b.tables
    if ta.equals(tb, ignore_provenance=True):
        byte_same = open(path_a, "rb").read() == open(path_b, "rb").read()
        note = "byte-identical" if byte_same else "differ only in provenance"
        return True, note

    moved = [name for name in TABLES
             if getattr(ta, name) != getattr(tb, name)]
    if not moved:
        return False, "no table differs -- top-level metadata or time_units"

    # WHICH COLUMN, not just which table. A table can differ in its row count
    # (a genuinely different simulation) or in one column (often metadata that
    # nothing downstream reads), and those need opposite responses. Reporting
    # only the table name sent me to diagnose this by hand the first time an
    # `individuals` table differed at an identical row count.
    details = []
    for name in moved:
        A, B = getattr(ta, name), getattr(tb, name)
        if len(A) != len(B):
            details.append(f"{name}: {len(A)} vs {len(B)} rows")
            continue
        cols = []
        for col in A.column_names:
            try:
                same = np.array_equal(np.asarray(getattr(A, col)),
                                      np.asarray(getattr(B, col)))
            except Exception:
                same = getattr(A, col) == getattr(B, col)
            if not same:
                cols.append(col)
        meta_only = bool(cols) and all("metadata" in c for c in cols)
        suffix = " (METADATA ONLY)" if meta_only else ""
        details.append(f"{name}: {len(A)} rows, columns differ: "
                       f"{', '.join(cols) or 'unknown'}{suffix}")
    return False, "; ".join(details)


def main(argv):
    if len(argv) < 2 or len(argv) % 2:
        raise SystemExit(__doc__)
    ok = True
    for i in range(0, len(argv), 2):
        a, b = argv[i], argv[i + 1]
        try:
            same, detail = compare(a, b)
        except Exception as e:              # a truncated or absent file
            print(f"ERROR\t{a}\t{b}\t{type(e).__name__}: {e}")
            ok = False
            continue
        print(f"{'SAME' if same else 'DIFFERENT'}\t{a}\t{b}\t{detail}")
        ok = ok and same
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
