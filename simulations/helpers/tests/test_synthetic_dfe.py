"""Tests for the drawn-effect-size path used by category H.

Two things are being checked, and they are different in kind:

  * helpers/synthetic_dfe.py draws from the right distribution, and draws the
    SAME coefficient for a position in every panel that carries it; and
  * create_gwas_files_and_phenotypes.py's pool predicates behave under
    --synthetic_dfe_effects, where no variant has a selection coefficient and a
    `selco != 0` filter would empty every pool.

Like test_causal_selection.py, this stubs the heavy simulation imports so the
pandas filters can be checked without installing SLiM's Python stack.
"""

import sys
import types

import numpy as np
import pandas as pd
import pytest

_HEAVY = ("tskit", "pyslim", "msprime", "tstrait", "seaborn",
          "matplotlib", "matplotlib.pyplot")

for _name in _HEAVY:
    if _name not in sys.modules:
        _mod = types.ModuleType(_name)
        if _name == "matplotlib":
            _mod.pyplot = types.ModuleType("matplotlib.pyplot")
            sys.modules["matplotlib.pyplot"] = _mod.pyplot
        sys.modules[_name] = _mod

import create_gwas_files_and_phenotypes as stage2  # noqa: E402
from helpers import synthetic_dfe  # noqa: E402


LO, HI = 5e5, 1.5e6  # the L = 2 Mb central window


def make_vars(positions, selco, maf):
    return pd.DataFrame({
        "id": np.arange(len(positions)),
        "type": 0,
        "selco": np.asarray(selco, dtype=float),
        "daf": np.asarray(maf, dtype=float),
        "time": 0.0,
        "position": np.asarray(positions, dtype=float),
        "maf": np.asarray(maf, dtype=float),
        "Vs": 0.0,
    })


def neutral_vars(positions, maf=0.2):
    """A panel from the H arm: every variant has selco == 0."""
    positions = np.asarray(positions, dtype=float)
    return make_vars(positions, np.zeros(len(positions)),
                     np.full(len(positions), maf, dtype=float))


# ------------------------------------------------------------------- the draw

def test_draw_returns_one_coefficient_and_type_per_element():
    s, t = synthetic_dfe.draw(1000, np.random.default_rng(0))
    assert s.shape == (1000,) and t.shape == (1000,)
    assert set(np.unique(t)) <= set(synthetic_dfe.MUTATION_TYPES)


def test_draw_of_zero_is_empty_not_an_error():
    """The flank pool can legitimately be empty at a high causal_min_maf."""
    s, t = synthetic_dfe.draw(0, np.random.default_rng(0))
    assert len(s) == 0 and len(t) == 0


def test_draw_never_returns_a_zero_coefficient():
    """beta = sqrt(|s|) * scaling, so a zero would be a trait with no genetic
    component -- and the whole point of dropping the neutral class is that it
    cannot happen."""
    s, _ = synthetic_dfe.draw(20000, np.random.default_rng(1))
    assert (s != 0).all()


def test_component_proportions_match_the_dfe():
    """The neutral class is removed and the remainder renormalized, which for
    this DFE is the g1 vector verbatim -- neutral was never a component, it is a
    separate mutation rate."""
    _, t = synthetic_dfe.draw(200000, np.random.default_rng(2))
    for mtype, expected in zip(synthetic_dfe.MUTATION_TYPES, synthetic_dfe.PROPORTIONS):
        seen = float((t == mtype).mean())
        assert seen == pytest.approx(expected, abs=4 * np.sqrt(expected / 200000) + 1e-4)


def test_component_signs_and_means_match_slim():
    """m1/m2 deleterious (negative), m3 positive, with the means the SLiM
    initializeMutationType calls and farm_create_orig_pop_e2.py both use."""
    s, t = synthetic_dfe.draw(400000, np.random.default_rng(3))
    m1, m2, m3 = (s[t == k] for k in synthetic_dfe.MUTATION_TYPES)
    assert (m1 < 0).all() and (m2 < 0).all() and (m3 > 0).all()
    assert np.abs(m1).mean() == pytest.approx(synthetic_dfe.M1_MEAN, rel=0.05)
    assert np.abs(m2).mean() == pytest.approx(synthetic_dfe.M2_MEAN, rel=0.25)
    assert m3.mean() == pytest.approx(synthetic_dfe.M3_MEAN, rel=0.30)


def test_the_draw_is_reproducible_from_the_seed():
    a, _ = synthetic_dfe.draw(500, np.random.default_rng(7))
    b, _ = synthetic_dfe.draw(500, np.random.default_rng(7))
    assert np.array_equal(a, b)


# ------------------------------------------------- one draw per position, shared

def test_a_position_gets_the_same_coefficient_in_every_panel():
    """In the selected categories selco is mutation metadata: identical in the
    GWAS and GTEx panels, with only the multiplier differing. A per-panel draw
    would give a locus causal in both two unrelated effect sizes."""
    gwas = neutral_vars([6e5, 7e5, 8e5])
    gtex = neutral_vars([7e5, 8e5, 9e5])
    smap, tmap = synthetic_dfe.assign_by_position([gwas, gtex], seed=11)
    g = synthetic_dfe.apply_to(gwas, smap, tmap)
    t = synthetic_dfe.apply_to(gtex, smap, tmap)
    for pos in (7e5, 8e5):
        assert (float(g.loc[g.position == pos, "selco"].iloc[0])
                == float(t.loc[t.position == pos, "selco"].iloc[0]))


def test_the_map_does_not_depend_on_the_order_the_frames_are_passed():
    a, b = neutral_vars([6e5, 8e5]), neutral_vars([7e5, 9e5])
    assert (synthetic_dfe.assign_by_position([a, b], seed=5)[0]
            == synthetic_dfe.assign_by_position([b, a], seed=5)[0])


def test_a_position_absent_from_the_map_raises_rather_than_defaulting():
    """Silently keeping the tree sequence's 0 would make a trait with no genetic
    component, which downstream reads as a colocalization failure rather than as
    a bug."""
    covered = neutral_vars([6e5])
    smap, tmap = synthetic_dfe.assign_by_position([covered], seed=1)
    with pytest.raises(KeyError):
        synthetic_dfe.apply_to(neutral_vars([6e5, 7e5]), smap, tmap)


def test_apply_to_does_not_mutate_its_input():
    frame = neutral_vars([6e5, 7e5])
    smap, tmap = synthetic_dfe.assign_by_position([frame], seed=2)
    synthetic_dfe.apply_to(frame, smap, tmap)
    assert (frame["selco"] == 0).all()


def test_apply_to_also_records_the_component_it_came_from():
    frame = neutral_vars([6e5, 7e5])
    smap, tmap = synthetic_dfe.assign_by_position([frame], seed=2)
    out = synthetic_dfe.apply_to(frame, smap, tmap)
    assert set(out["type"]) <= set(synthetic_dfe.MUTATION_TYPES)


# --------------------------------------------------------- the pool predicates

def test_every_pool_is_empty_without_the_synthetic_flag():
    """The failure this guards against: running an H tree sequence through the
    ordinary path yields no causal variants at all, and 50 traits silently
    become 0."""
    v = neutral_vars([6e5, 7e5, 1.4e6])
    assert stage2.causal_eligible(v, LO, HI, 0, synthetic=False).empty
    flanks = neutral_vars([1e5, 1.9e6])
    assert stage2.flank_eligible(flanks, 0, LO, HI, synthetic=False).empty


def test_the_flank_pool_takes_neutral_variants_under_the_flag():
    v = neutral_vars([1e5, 1.9e6, 8e5])          # two flanks and one central
    flank = stage2.flank_eligible(v, 0, LO, HI, synthetic=True)
    assert sorted(flank["position"]) == [1e5, 1.9e6]


def test_the_flank_pool_applies_the_causal_floor_and_drops_monomorphic_sites():
    v = make_vars([1e5, 1.2e5, 1.9e6], [0.0] * 3, [0.0, 0.005, 0.3])
    assert list(stage2.flank_eligible(v, 0, LO, HI, synthetic=True)["position"]) \
        == [1.2e5, 1.9e6]                         # maf == 0 gone even at floor 0
    assert list(stage2.flank_eligible(v, 0.01, LO, HI, synthetic=True)["position"]) \
        == [1.9e6]


def test_the_central_floor_is_measured_in_the_panel_it_is_given():
    """causal_eligible only ever sees one panel; the GTEx cross-check is the
    caller's job and only happens on the non-top-up path."""
    v = make_vars([6e5, 7e5], [0.0, 0.0], [0.02, 0.005])
    assert list(stage2.causal_eligible(v, LO, HI, 0.01, synthetic=True)["position"]) == [6e5]


def test_gtex_topup_candidates_respect_the_floor_and_the_flag():
    """A top-up locus is as causative as a drawn one, so a floor that excluded the
    latter but not the former would put variants in the GTEx causal set that the
    config forbade."""
    chosen = neutral_vars([6e5])
    gtex = make_vars([6e5, 7e5, 8e5], [0.0] * 3, [0.3, 0.005, 0.3])
    shared, topup = stage2.select_gtex_topup(
        gtex, chosen, LO, HI, 3, seed=1, causal_min_maf=0.01, synthetic=True)
    assert list(shared["position"]) == [6e5]
    assert list(topup["position"]) == [8e5]       # 7e5 is below the floor


def test_the_power_weight_uses_the_drawn_coefficient():
    """The 'draw s for the whole pool, then compute power' ordering. If the draw
    did not land on the frame before select_central_power, every candidate would
    weigh sqrt(0) * scaling == 0 and the inclusion probabilities would be
    undefined."""
    rng = np.random.default_rng(4)
    pos = np.sort(rng.choice(np.arange(int(LO) + 1, int(HI)), size=600, replace=False))
    pool = neutral_vars(pos, maf=0.3)
    smap, tmap = synthetic_dfe.assign_by_position([pool], seed=9)
    with_s = synthetic_dfe.apply_to(pool, smap, tmap)

    chosen, diag = stage2.select_central_power(
        with_s, 50, scaling=35, sampling_n=8000, sig_p=5e-8,
        min_power=0.05, min_pool_multiple=2, seed=3, label="synthetic")
    assert chosen.shape[0] == 50
    assert (diag["beta"] > 0).all()
    assert diag["power"].max() > 0
    # And the drawn coefficient travelled with the row, so beta is derived from it.
    row = chosen.iloc[0]
    assert float(np.sqrt(abs(row["selco"])) * 35) == pytest.approx(
        float(diag.loc[diag["position"] == row["position"], "beta"].iloc[0]))


# ------------------------------------------- the causal-MAF path token

def test_the_maf_token_agrees_between_the_config_and_the_argparse_float():
    """paths.stage2_inner and create_gwas_files_and_phenotypes.py's output-directory
    f-string must render causal_min_maf identically, or Snakemake declares one
    directory as the rule's output and the script writes another -- failing only
    after stage 2 has done all its work.

    They agreed by accident while every value was 0.01 or 0.001. Zero is where it
    broke: the config carries YAML int 0, the script parses --causal_min_maf as a
    float, and str() renders those "0" and "0.0"."""
    from helpers import paths
    for value in (0, 0.0, 0.01, 0.001, 0.05, 0.5):
        assert paths.causal_maf_token(value) == paths.causal_maf_token(float(value))
    assert paths.causal_maf_token(0) == "0"
    assert paths.causal_maf_token(0.0) == "0"
    assert paths.causal_maf_token(0.01) == "0.01"
    assert paths.causal_maf_token(0.001) == "0.001"


def test_the_legacy_floor_still_emits_no_path_segment():
    """0.01 is the path-suppression sentinel, NOT the default. Every output that
    predates the knob was produced at 0.01 and must stay reachable."""
    from helpers import paths
    assert paths.causal_maf_segment({"causal_min_maf": 0.01}) == ""
    assert paths.causal_maf_segment({"causal_min_maf": 0}) == ".cmaf_0"
    assert paths.causal_maf_segment({"causal_min_maf": 0.0}) == ".cmaf_0"
    assert paths.causal_maf_segment({}) == ".cmaf_0"          # the new default
    assert paths.causal_maf_segment({"causal_min_maf": 0.001}) == ".cmaf_0.001"
