"""Categories K and L: A's and E's genomes with H's and I's effect model.

The one thing that can go wrong silently.

`--synthetic_dfe_effects` used to mean "do not require selco != 0 on the causal
pool". On a coalescent that is the same as "require selco == 0", because a
coalescent has no selected variants -- so H and I never exercised the
difference. K and L run on FORWARD genomes, where the difference is the whole
category: dropping the filter would let them draw variants that already carry a
selection coefficient and then assign a second, unrelated one on top, producing
an arm that is neither A nor H and looks like both.

These tests pin the predicate, and pin that it stayed inert for H and I.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from helpers import paths


SIM_DIR = Path(paths.__file__).parent.parent
SRC = (SIM_DIR / "create_gwas_files_and_phenotypes.py").read_text()


def _load(fn):
    """Exec one top-level function out of the stage-2 script.

    Importing the module would parse argv and pull in tstrait/tskit, which the
    stubbed test environment does not carry.
    """
    m = re.search(rf'^def {fn}\(.*?\n(?=^def |\Z)', SRC, re.S | re.M)
    assert m, f"{fn} not found"
    ns = {'np': np, 'pd': pd}
    exec(compile(m.group(0), fn, 'exec'), ns)
    return ns[fn]


causal_eligible = _load('causal_eligible')
flank_eligible = _load('flank_eligible')


def pool():
    """A forward-genome pool: some selected variants, some strictly neutral."""
    return pd.DataFrame({
        'maf':      [0.20, 0.20, 0.20, 0.20, 0.20],
        'position': [7e5,  8e5,  9e5,  1.0e6, 1.1e6],
        'selco':    [0.0, -0.01, 0.0,  0.02,  0.0],
    })


# ------------------------------------------------------- the causal predicate

def test_the_selected_arms_still_draw_only_selected_variants():
    """A and E: unchanged behaviour."""
    got = causal_eligible(pool(), 5e5, 1.5e6, 0.0, False)
    assert sorted(got.selco.tolist()) == [-0.01, 0.02]


def test_the_synthetic_arms_draw_only_NEUTRAL_variants():
    """K and L: the category definition. Not 'all variants'."""
    got = causal_eligible(pool(), 5e5, 1.5e6, 0.0, True)
    assert got.selco.tolist() == [0.0, 0.0, 0.0]
    assert len(got) == 3


def test_the_two_pools_are_disjoint():
    """Whatever else changes, an arm cannot be both. If these ever overlap, a
    variant could be drawn as causal under one rule and re-coefficiented under
    the other."""
    p = pool()
    sel = set(causal_eligible(p, 5e5, 1.5e6, 0.0, False).index)
    neu = set(causal_eligible(p, 5e5, 1.5e6, 0.0, True).index)
    assert not (sel & neu)
    assert sel | neu == set(p.index)


def test_the_flank_pool_follows_the_same_rule():
    """The flank loci are causal too -- GTEx-only, but causal -- so they must be
    neutral in K/L for the same reason the central ones are."""
    fl = pd.DataFrame({'maf': [0.2, 0.2, 0.2, 0.2],
                       'position': [1e5, 1e5, 1.9e6, 1.9e6],
                       'selco': [0.0, -0.01, 0.0, 0.03]})
    assert flank_eligible(fl, 0.0, 5e5, 1.5e6, True).selco.tolist() == [0.0, 0.0]
    assert sorted(flank_eligible(fl, 0.0, 5e5, 1.5e6, False).selco.tolist()) == [-0.01, 0.03]


def test_a_coalescent_genome_is_unaffected_by_the_change():
    """H and I: every variant has selco == 0, so 'all' and 'the neutral ones'
    are the same set and the tightening cannot move their results.

    This is the regression that would otherwise be invisible -- it only shows up
    as different H/I numbers months later.
    """
    coal = pd.DataFrame({'maf': [0.2] * 4,
                         'position': [7e5, 8e5, 9e5, 1e6],
                         'selco': [0.0] * 4})
    got = causal_eligible(coal, 5e5, 1.5e6, 0.0, True)
    assert len(got) == len(coal)


def test_the_maf_floor_and_window_still_compose():
    """The selco rule must not have swallowed the other two predicates."""
    p = pool()
    p.loc[0, 'maf'] = 0.0005
    p.loc[2, 'position'] = 2e5          # outside the central window
    got = causal_eligible(p, 5e5, 1.5e6, 0.001, True)
    assert got.index.tolist() == [4]


# ----------------------------------------------------------- the wiring

def test_the_predicate_runs_before_the_drawn_coefficients_land():
    """synthetic_dfe.apply_to OVERWRITES selco with the drawn values. If it ran
    first, `selco == 0` would filter on the draw rather than on the genome and
    would select nothing at all."""
    i_pool = SRC.index("pool = causal_eligible(gwas_vars, trait_pos_lo")
    i_apply = SRC.index("pool = synthetic_dfe.apply_to(pool, *synth)")
    assert i_pool < i_apply


@pytest.mark.parametrize("letter,species", [("K", "human"), ("L", "cattle")])
def test_the_categories_are_known_to_the_summariser(letter, species):
    from helpers import summarize_coloc
    assert summarize_coloc.DEMOGRAPHY[letter] == species


@pytest.mark.parametrize("letter,seed", [("K", 111), ("L", 121)])
def test_the_seed_bands_are_disjoint_from_every_earlier_letter(letter, seed):
    """seed = 10*letter_index + replicate, K the 11th letter and L the 12th."""
    idx = ord(letter) - ord('A') + 1
    assert 10 * idx + 1 == seed
    assert seed > 109          # clear of J's band, and of the 11-99 above it


def test_the_submit_script_gives_them_the_override_and_a_fresh_build(guard_letters):
    sh = (SIM_DIR / "submit_2Mb_r3_cmaf_fm001.sh").read_text()
    assert 'K|L) echo "synthetic_dfe_effects=True"' in sh
    assert 'K) echo "11${n}" ;; L) echo "12${n}"' in sh
    # They must NOT adopt stage 2 -- adopting A's or E's would hand them
    # selected causal loci, which is the one thing they are defined by lacking.
    assert 'STAGE2_SRC["$ID"]=""' in sh
    assert {"K", "L"} <= guard_letters(sh, 'STAGE2_SRC["$ID"]=""')
    # and they must build their own stage 1, since the seeds are new
    assert {"H", "I", "J", "K", "L"} <= guard_letters(sh, 'STAGE1_SRC["$ID"]=""')
