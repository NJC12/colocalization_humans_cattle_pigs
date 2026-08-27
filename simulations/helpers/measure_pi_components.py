#!/usr/bin/env python
"""Measure pi's neutral and selected halves on a stage-1 tree sequence.

This is the one measurement categories O and P need before they can run: the
keep fraction that halves pi is SOLVED from it (helpers/neutral_thinning.py:
keep_fraction_for_pi_target), not searched for, so one pass per replicate fixes
the parameter for the whole arm.

It replays the stage-2 preamble rather than reading the raw file, because the
tree sequence stage 2 actually thins is not the one on disk:

    cattle  add_neutral() overlays the 8.4e-9 neutral class HERE, in stage 2, in
            two passes across the Q=1 -> Q=0.01 handoff. Measuring the raw
            .full.ts would see only the 5.6e-9 DFE mutations and report a
            selected share of ~1.0.
    human   the overlay already happened in stage 1, so this is a no-op --
            already_includes_neutral is True in every human config.

Run on O2, where the tree sequences live. One replicate:

    python helpers/measure_pi_components.py \
        --ts_file /n/scratch/.../A1/stage1/hts_11.ts \
        --species human --Q_scaling 10 --already_includes_neutral True \
        --seed 11 --label A1

    python helpers/measure_pi_components.py \
        --ts_file /n/scratch/.../E1/stage1/farm_selection_from_ep8...full.ts \
        --species cattle --Q_scaling 0.01 --handoff_ticks 2400 \
        --deep_Q_scaling 1 --seed 51 --label E1

Writes one TSV row per invocation (header only when the output file is new), so
several replicates can append to one table:

    label species pi_total pi_neutral pi_selected selected_share \
        n_sites n_neutral n_selected keep_fraction_for_half_pi

Take the MEAN keep fraction across a species' replicates, round it, and pin it in
that category's config as `neutral_keep_fraction`. Pinning one number per species
rather than one per replicate is what makes O and P categories rather than ten
separate conditions.
"""

import argparse
import os
import sys

# sys.path[0] is helpers/, so add the repo root the way the stage-2 script's own
# `from helpers import ...` line assumes.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import tskit  # noqa: E402

# Importing the stage-2 script is safe: everything it does at run time is behind
# `if __name__ == '__main__'`. Reusing its functions rather than reimplementing
# them is the point -- a private copy of add_neutral could drift from the one the
# pipeline runs, and then the calibration would describe a tree sequence nobody
# builds.
import create_gwas_files_and_phenotypes as stage2  # noqa: E402
from helpers import neutral_thinning  # noqa: E402


COLUMNS = [
    "label", "species", "ts_file",
    "pi_total", "pi_neutral", "pi_selected", "selected_share",
    "n_sites", "n_neutral", "n_selected",
    "keep_fraction_for_half_pi",
]


def prepared_ts(args):
    """The tree sequence as stage 2 sees it, immediately before thinning."""
    ts = tskit.load(args.ts_file)
    if args.species == "cattle":
        # Unconditional on the cattle branch, exactly as stage 2 does it --
        # already_includes_neutral is not consulted there, and every cattle
        # config sets it True, so honouring it would switch the overlay off.
        ts = stage2.remove_fixed(stage2.add_neutral(
            ts, Q=1 / args.Q_scaling, seed=args.seed,
            handoff_ticks=args.handoff_ticks,
            deep_Q_scaling=args.deep_Q_scaling))
    elif not args.already_includes_neutral:
        ts = stage2.remove_fixed(stage2.add_neutral(ts, seed=args.seed))
    # generate_nucleotides/convert_alleles are deliberately skipped: they set
    # REF/ALT letters and change no genotype, position or frequency, so they
    # cannot move pi, and they are the slow part of the preamble.
    return stage2.remove_position_zero(ts)


def get_arguments():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--ts_file', required=True, help='Stage-1 tree sequence.')
    p.add_argument('--species', required=True, choices=('human', 'cattle'))
    p.add_argument('--Q_scaling', type=float, required=True,
                   help='The config key of the same name (human 10, cattle 0.01).')
    p.add_argument('--handoff_ticks', type=float, default=None,
                   help='Cattle only: ticks since the Q=1 -> Q=0.01 handoff (2400).')
    p.add_argument('--deep_Q_scaling', type=float, default=None,
                   help='Cattle only: the Q the deep phase ran under (1).')
    p.add_argument('--already_includes_neutral', type=stage2.str2bool, default=True,
                   help='Human only; True for every config in config/.')
    p.add_argument('--seed', type=int, required=True,
                   help="This replicate's stage2_seed, so the cattle overlay is "
                        "the same one stage 2 will lay down.")
    p.add_argument('--target', type=float, default=0.5,
                   help='Fraction of pi the thinned arm should carry (default 0.5).')
    p.add_argument('--label', default='', help='Replicate id, e.g. A1. For the TSV.')
    p.add_argument('--out_tsv', default=None,
                   help='Append a row here as well as printing it.')
    return p.parse_args()


if __name__ == '__main__':
    args = get_arguments()
    ts = prepared_ts(args)

    pi_neutral, pi_selected = neutral_thinning.pi_components(ts)
    pi_total = pi_neutral + pi_selected
    n_neutral = int(neutral_thinning.neutral_site_ids(ts).size)

    try:
        k = neutral_thinning.keep_fraction_for_pi_target(
            pi_neutral, pi_selected, target=args.target)
        k_str = f'{k:.6g}'
    except ValueError as exc:
        # Do not clamp. A negative k means the selected class alone already
        # carries the target, and the arm as specified is not buildable -- that
        # has to be visible, not rounded away.
        print(f'INFEASIBLE: {exc}', file=sys.stderr)
        k_str = 'NA'

    row = {
        'label': args.label, 'species': args.species, 'ts_file': args.ts_file,
        'pi_total': f'{pi_total:.6g}',
        'pi_neutral': f'{pi_neutral:.6g}',
        'pi_selected': f'{pi_selected:.6g}',
        'selected_share': f'{pi_selected / pi_total:.6g}' if pi_total else 'NA',
        'n_sites': ts.num_sites,
        'n_neutral': n_neutral,
        'n_selected': ts.num_sites - n_neutral,
        'keep_fraction_for_half_pi': k_str,
    }

    print('\t'.join(COLUMNS))
    print('\t'.join(str(row[c]) for c in COLUMNS))

    if args.out_tsv:
        new = not os.path.exists(args.out_tsv) or os.path.getsize(args.out_tsv) == 0
        with open(args.out_tsv, 'a') as fh:
            if new:
                fh.write('\t'.join(COLUMNS) + '\n')
            fh.write('\t'.join(str(row[c]) for c in COLUMNS) + '\n')
