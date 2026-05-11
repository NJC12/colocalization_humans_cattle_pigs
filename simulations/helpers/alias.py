#!/usr/bin/env python3
"""Symlink legacy filenames into the canonical layout the Snakefile expects.

Existing pre-Snakemake outputs use names like
``revision_farm_selection_mult_100_gen_23_muts_26_bottlenecked_sd24.full.ts``
while the new pipeline expects e.g.
``cattle_sel_bot_mult100_gen23_muts26_sd24.full.ts``. Run this once before the
first Snakemake invocation to add symlinks so the multi-directory search
returns matches without copying data.

Idempotent. Use ``--dry-run`` to preview.
"""

import argparse
import os
import re
import sys
from pathlib import Path


# (legacy_pattern, new_basename_template, optional flavor map) per known
# legacy form. Templates are .format()-style with named groups from the regex.
#
# Note: human and cattle baseline filenames already match the canonical naming
# in helpers/paths.py, so they don't need aliases. Only the cattle+selection
# revision files have a divergent legacy form.
LEGACY_RULES = [
    # cattle + selection (bottlenecked + not-bottlenecked)
    (
        re.compile(
            r"^revision_farm_selection_mult_(?P<mult>\d+)_gen_(?P<gen>\d+)_muts_(?P<muts>\d+)_(?P<flavor>bottlenecked|not_bottlenecked)_sd(?P<sd>\d+)\.(?P<rest>.+)$"
        ),
        "{base}_mult{mult}_gen{gen}_muts{muts}_sd{sd}.{rest}",
        {"bottlenecked": "cattle_sel_bot", "not_bottlenecked": "cattle_sel_nobot"},
    ),
]


def plan_aliases(search_dir):
    out = []
    for entry in Path(search_dir).iterdir():
        if not entry.is_file() and not entry.is_symlink():
            continue
        for regex, template, flavor_map in LEGACY_RULES:
            m = regex.match(entry.name)
            if not m:
                continue
            d = m.groupdict()
            if flavor_map is not None:
                d["base"] = flavor_map[d.pop("flavor")]
            new_name = template.format(**d)
            if new_name == entry.name:
                continue
            link = entry.parent / new_name
            out.append((entry, link))
            break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--search-dir", required=True,
                    help="Directory to scan for legacy filenames.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print the symlink plan without creating links.")
    args = ap.parse_args()

    plan = plan_aliases(args.search_dir)
    if not plan:
        print(f"No legacy filenames found in {args.search_dir}.")
        return 0

    for src, link in plan:
        if link.exists() or link.is_symlink():
            if link.is_symlink() and os.readlink(link) == src.name:
                continue
            print(f"SKIP (target exists, differs): {link}")
            continue
        if args.dry_run:
            print(f"DRY:   ln -s {src.name} {link}")
        else:
            link.symlink_to(src.name)
            print(f"linked {link} -> {src.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
