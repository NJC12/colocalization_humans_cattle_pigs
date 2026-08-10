"""Tests for the trait-selection logic in create_gwas_files_and_phenotypes.py.

That module imports the whole simulation stack (tskit, pyslim, msprime, tstrait)
at the top even though the selection functions touch none of it, so the fixture
below stubs those imports out.  The alternative -- installing SLiM's Python stack
just to check a pandas filter -- would put these tests out of reach on any machine
that only ever reads the results.

Covers the new power-sampling path and the pre-existing uniform one, which had no
tests at all.
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


LO, HI = 5e5, 1.5e6  # the L = 2 Mb central window


def make_vars(positions, selco, maf):
    """A minimal *_vars_* frame: the columns the selection code reads."""
    return pd.DataFrame({
        "id": np.arange(len(positions)),
        "type": 2,
        "selco": np.asarray(selco, dtype=float),
        "daf": np.asarray(maf, dtype=float),
        "time": 0.0,
        "position": np.asarray(positions, dtype=float),
        "maf": np.asarray(maf, dtype=float),
        "Vs": 0.0,
    })


def detectable_pool(n=600, seed=7):
    """A pool shaped like the real one: a detectable third, an invisible rest.

    200 variants clear power 0.05 at sampling_n = 8000, which is exactly what the
    default guard (2 x 50 central traits = 100) has slack against.
    """
    rng = np.random.default_rng(seed)
    pos = np.sort(rng.choice(np.arange(int(LO) + 1, int(HI)), size=n, replace=False))
    strong = np.arange(n) % 3 == 0
    maf = np.where(strong, rng.uniform(0.1, 0.5, n), rng.uniform(1e-4, 1e-3, n))
    selco = np.where(strong, rng.uniform(1e-3, 1e-2, n), rng.uniform(1e-8, 1e-7, n))
    return make_vars(pos, selco, maf)


def graded_pool(n=800, seed=13):
    """A pool with a continuous frequency/effect spread, for the sampling_n sweep."""
    rng = np.random.default_rng(seed)
    pos = np.sort(rng.choice(np.arange(int(LO) + 1, int(HI)), size=n, replace=False))
    maf = 10.0 ** rng.uniform(-3.5, np.log10(0.5), n)
    selco = 10.0 ** rng.uniform(-5.0, -2.0, n)
    return make_vars(pos, selco, maf)


# ------------------------------------------------------------------ the pool

def test_central_pool_drops_the_flanks_and_the_neutral_variants():
    v = make_vars([4e5, 6e5, 7e5, 1.4e6, 1.6e6, 8e5],
                  [1e-3, 1e-3, 0.0, 1e-3, 1e-3, 1e-3],
                  [0.2] * 6)
    pool = stage2.causal_eligible(v, LO, HI, 0)
    assert sorted(pool["position"]) == [6e5, 8e5, 1.4e6]


def test_central_pool_keeps_every_polymorphic_variant_at_a_zero_floor():
    """causal_min_maf = 0 is the default, and means no frequency floor at all."""
    v = make_vars([6e5, 7e5], [1e-3, 1e-3], [1e-5, 0.4])
    assert stage2.causal_eligible(v, LO, HI, 0).shape[0] == 2


def test_central_pool_applies_the_floor_when_one_is_given():
    """The floor gates the central pool in BOTH sampling schemes -- it used to be
    dropped under power on the grounds that the detection-power weight subsumed it,
    which made 'the pool at causal_min_maf=0.01' mean two different things."""
    v = make_vars([6e5, 7e5], [1e-3, 1e-3], [1e-5, 0.4])
    assert list(stage2.causal_eligible(v, LO, HI, 0.01)["position"]) == [7e5]


def test_central_pool_drops_variants_monomorphic_in_this_panel():
    """remove_fixed runs before the GWAS/GTEx split, so a site can be fixed within
    one panel. Its power is undefined, not small, and an effect on it would make a
    phenotype of pure noise -- detection_power rejects maf == 0 outright.

    Load-bearing at causal_min_maf = 0, where `maf >= floor` alone would keep it."""
    v = make_vars([6e5, 7e5, 8e5], [1e-3] * 3, [0.0, 0.3, 0.0])
    assert list(stage2.causal_eligible(v, LO, HI, 0)["position"]) == [7e5]


def test_central_pool_keeps_neutral_variants_under_synthetic_effects():
    """Category H: nothing in the tree sequence carries a selection coefficient, so
    a `selco != 0` pool would be empty and every trait would be lost."""
    v = make_vars([6e5, 7e5, 8e5], [0.0, 0.0, 0.0], [0.2] * 3)
    assert stage2.causal_eligible(v, LO, HI, 0, synthetic=True).shape[0] == 3
    assert stage2.causal_eligible(v, LO, HI, 0, synthetic=False).empty


def test_the_draw_never_sees_a_monomorphic_variant():
    pool = detectable_pool()
    dirty = pd.concat([pool, make_vars([1.45e6], [1e-2], [0.0])])
    cleaned = stage2.causal_eligible(dirty, LO, HI, 0)
    chosen, _ = _draw(cleaned)          # would raise inside detection_power otherwise
    assert 1.45e6 not in set(chosen["position"])


def test_a_locus_monomorphic_in_the_gtex_panel_is_not_a_partner():
    """'Present in the GTEx samples' has to mean segregating there -- an effect on
    a monomorphic column is not an eQTL."""
    chosen = make_vars([6e5, 7e5], [1e-3] * 2, [0.3] * 2)
    gtex = make_vars([6e5, 7e5, 9e5], [1e-3] * 3, [0.0, 0.3, 0.3])
    shared, topup = stage2.select_gtex_topup(gtex, chosen, LO, HI, 2, seed=1)
    assert list(shared["position"]) == [7e5]
    assert list(topup["position"]) == [9e5]   # and never the monomorphic 6e5


def test_central_pool_boundary_positions_go_to_the_flanks():
    # select_flank_gtex uses <= / >=, so the pool must use strict > / < or the
    # boundary variants would land in both sets.
    v = make_vars([LO, HI], [1e-3, 1e-3], [0.3, 0.3])
    assert stage2.causal_eligible(v, LO, HI, 0).empty


# --------------------------------------------------------- the power-weighted draw

def _draw(pool, n=50, seed=11, **kw):
    kw.setdefault("scaling", 10)
    kw.setdefault("sampling_n", 8000)
    kw.setdefault("sig_p", 5e-8)
    kw.setdefault("min_power", 0.05)
    kw.setdefault("min_pool_multiple", 2)
    return stage2.select_central_power(pool, n, seed=seed, **kw)


def test_the_draw_returns_the_requested_number_of_pool_rows():
    pool = detectable_pool()
    chosen, diag = _draw(pool)
    assert chosen.shape[0] == 50
    assert set(chosen["position"]).issubset(set(pool["position"]))
    assert diag.shape[0] == pool.shape[0]  # diagnostics cover the whole pool


def test_the_draw_keeps_the_rows_in_position_order():
    chosen, _ = _draw(detectable_pool())
    assert list(chosen["position"]) == sorted(chosen["position"])


def test_the_draw_is_reproducible_from_the_seed():
    pool = detectable_pool()
    a, _ = _draw(pool, seed=99)
    b, _ = _draw(pool, seed=99)
    c, _ = _draw(pool, seed=100)
    assert list(a["position"]) == list(b["position"])
    assert list(a["position"]) != list(c["position"])


def test_the_draw_is_enriched_for_detectable_variants():
    pool = detectable_pool()
    _, diag = _draw(pool)
    sel = diag.loc[diag["selected"]]
    assert sel["power"].median() > 0.9
    assert diag["power"].median() < 0.01
    # Half the pool is invisible; essentially none of it should be drawn.
    assert (sel["power"] < 0.05).sum() <= 1


def test_the_diagnostics_flag_exactly_the_drawn_rows():
    pool = detectable_pool()
    chosen, diag = _draw(pool)
    assert diag["selected"].sum() == chosen.shape[0]
    assert sorted(diag.loc[diag["selected"], "position"]) == sorted(chosen["position"])


def test_a_larger_sampling_n_reaches_further_down_the_frequency_spectrum():
    """The knob's whole purpose: ask what a bigger study would find and the drawn
    loci get rarer, without the simulated GWAS itself changing."""
    pool = graded_pool()
    small, _ = _draw(pool, sampling_n=8000)
    large, _ = _draw(pool, sampling_n=2_000_000)
    assert large["maf"].median() < small["maf"].median()


def test_the_draw_refuses_a_pool_with_too_little_power():
    pool = detectable_pool()
    with pytest.raises(SystemExit, match="Refusing to draw"):
        _draw(pool, sampling_n=50)


def test_the_refusal_names_the_counts_and_the_knobs_to_turn():
    pool = detectable_pool()
    with pytest.raises(SystemExit) as exc:
        _draw(pool, sampling_n=50)
    msg = str(exc.value)
    assert "sampling_gwas_n=50" in msg and "--sampling_min_power" in msg


def test_the_guard_scales_with_the_requested_pool_multiple():
    pool = detectable_pool()
    _draw(pool, n=50, min_pool_multiple=2)          # 200 detectable, 100 required
    with pytest.raises(SystemExit):
        _draw(pool, n=50, min_pool_multiple=5)      # 250 required


def test_the_recipient_maf_drives_the_weight_when_effects_are_redistributed():
    """Categories B/D: the effect lands on a neutral recipient, so the recipient's
    frequency is the only one any downstream test sees."""
    pool = make_vars([6e5, 7e5, 8e5, 9e5], [1e-2] * 4, [0.4] * 4)
    # Donors all look identical and common; two of their recipients are singletons.
    recipient_maf = {6e5: 0.4, 7e5: 1e-4, 8e5: 0.4, 9e5: 1e-4}
    _, diag = _draw(pool, n=2, min_pool_multiple=1, maf_by_position=recipient_maf)
    by_pos = dict(zip(diag["position"], diag["power"]))
    assert by_pos[6e5] > 0.9 and by_pos[8e5] > 0.9
    assert by_pos[7e5] < 0.05 and by_pos[9e5] < 0.05
    assert list(diag.loc[diag["selected"], "position"]) == [6e5, 8e5]


# ---------------------------------------------------------- the GTEx central set

def test_the_gtex_set_shares_what_it_can_and_tops_up_the_rest():
    chosen = make_vars([6e5, 7e5, 8e5, 9e5], [1e-3] * 4, [0.3] * 4)
    # The GTEx panel is missing two of the drawn loci and carries four spares.
    gtex = make_vars([6e5, 8e5, 1.0e6, 1.1e6, 1.2e6, 1.3e6],
                     [1e-3] * 6, [0.3] * 6)
    shared, topup = stage2.select_gtex_topup(gtex, chosen, LO, HI, 4, seed=1)
    assert sorted(shared["position"]) == [6e5, 8e5]
    assert topup.shape[0] == 2
    assert set(topup["position"]).isdisjoint({6e5, 7e5, 8e5, 9e5})
    assert shared.shape[0] + topup.shape[0] == 4


def test_the_top_up_never_reuses_a_locus_the_gwas_already_drew():
    chosen = make_vars([6e5, 7e5], [1e-3] * 2, [0.3] * 2)
    gtex = make_vars([6e5, 7e5, 8e5], [1e-3] * 3, [0.3] * 3)
    shared, topup = stage2.select_gtex_topup(gtex, chosen, LO, HI, 3, seed=1)
    assert sorted(shared["position"]) == [6e5, 7e5]
    assert list(topup["position"]) == [8e5]


def test_the_top_up_pool_excludes_flanks_and_neutral_variants():
    chosen = make_vars([6e5], [1e-3], [0.3])
    gtex = make_vars([6e5, 3e5, 1.8e6, 9e5, 1.0e6],
                     [1e-3, 1e-3, 1e-3, 0.0, 1e-3], [0.3] * 5)
    _, topup = stage2.select_gtex_topup(gtex, chosen, LO, HI, 5, seed=1)
    assert list(topup["position"]) == [1.0e6]


def test_nothing_is_topped_up_when_every_locus_is_shared():
    chosen = make_vars([6e5, 7e5], [1e-3] * 2, [0.3] * 2)
    gtex = make_vars([6e5, 7e5, 8e5], [1e-3] * 3, [0.3] * 3)
    shared, topup = stage2.select_gtex_topup(gtex, chosen, LO, HI, 2, seed=1)
    assert shared.shape[0] == 2 and topup.empty


def test_sharing_is_decided_on_the_recipient_position_when_redistributing():
    """combine_phenos_to_df resolves the RECIPIENT's site id in the GTEx panel, so
    that is the position whose presence decides whether a locus can be shared."""
    chosen = make_vars([6e5, 7e5], [1e-3] * 2, [0.3] * 2)
    redist = {600000: 9e5, 700000: 1.0e6}       # donor -> neutral recipient
    gtex = make_vars([6e5, 9e5, 1.1e6], [1e-3, 0.0, 1e-3], [0.3] * 3)
    shared, topup = stage2.select_gtex_topup(gtex, chosen, LO, HI, 2, seed=1,
                                             redistribution=redist)
    # 6e5's recipient (9e5) is in the panel; 7e5's (1.0e6) is not. The donor frame
    # is what is returned, because the same map remaps it on the GTEx side too.
    assert list(shared["position"]) == [6e5]
    assert list(topup["position"]) == [1.1e6]


def test_no_two_gtex_traits_can_end_up_with_the_same_name():
    """A drawn locus's DONOR position is still an ordinary GTEx variant, so without
    the trait-position exclusion the top-up can pick it and remap it straight onto
    the recipient the drawn locus already occupies."""
    chosen = make_vars([6e5, 7e5], [1e-3] * 2, [0.3] * 2)
    redist = {600000: 9e5, 700000: 1.0e6}
    gtex = make_vars([6e5, 7e5, 9e5, 1.0e6, 1.1e6],
                     [1e-3, 1e-3, 0.0, 0.0, 1e-3], [0.3] * 5)
    shared, topup = stage2.select_gtex_topup(gtex, chosen, LO, HI, 3, seed=1,
                                             redistribution=redist)
    names = [redist.get(int(p), p)
             for p in list(shared["position"]) + list(topup["position"])]
    assert len(names) == len(set(names)) == 3


# ------------------------------------------------------------- the partner table

def test_the_partner_table_marks_every_locus_shared_under_uniform_sampling():
    chosen = make_vars([6e5, 7e5], [1e-3] * 2, [0.3] * 2)
    tbl = stage2.trait_partner_table(chosen, chosen)
    assert list(tbl["gwas_trait"]) == ["tr600000", "tr700000"]
    assert list(tbl["gtex_trait"]) == ["tr600000", "tr700000"]
    assert tbl["shared"].all()


def test_the_partner_table_leaves_the_gtex_id_empty_when_there_is_no_partner():
    chosen = make_vars([6e5, 7e5], [1e-3] * 2, [0.3] * 2)
    shared = chosen.iloc[:1]
    tbl = stage2.trait_partner_table(chosen, shared)
    assert list(tbl["shared"]) == [True, False]
    assert list(tbl["gtex_trait"]) == ["tr600000", ""]


def test_the_partner_table_names_traits_for_the_recipient_when_redistributing():
    chosen = make_vars([6e5], [1e-3], [0.3])
    tbl = stage2.trait_partner_table(chosen, chosen, redistribution={600000: 9e5})
    assert list(tbl["gwas_trait"]) == ["tr900000"]
    assert list(tbl["gtex_trait"]) == ["tr900000"]


# ---------------------------------------------- the pre-existing uniform helpers

def test_subsample_traits_returns_the_pool_unchanged_when_it_is_too_small():
    v = make_vars([6e5, 7e5], [1e-3] * 2, [0.3] * 2)
    assert stage2.subsample_traits(v, 5, seed=1).equals(v)
    assert stage2.subsample_traits(v, None, seed=1).equals(v)


def test_subsample_traits_keeps_the_rows_in_position_order():
    v = make_vars(np.arange(6e5, 6e5 + 100), [1e-3] * 100, [0.3] * 100)
    kept = stage2.subsample_traits(v, 10, seed=3)
    assert kept.shape[0] == 10
    assert list(kept["position"]) == sorted(kept["position"])


# ------------------------------------------------------------ the composed flow

def _panels(seed=5, n_central=900, n_flank=300, drop_frac=0.15):
    """A GWAS and a GTEx vars frame, with the GTEx panel missing some GWAS variants.

    Stands in for the real split, where the 1000-individual GTEx panel does not
    segregate every variant the 8000-individual GWAS panel does.
    """
    rng = np.random.default_rng(seed)
    central = graded_pool(n_central, seed=seed)
    flank_pos = np.concatenate([rng.uniform(0, LO, n_flank // 2),
                                rng.uniform(HI, 2e6, n_flank // 2)])
    flank = make_vars(np.sort(flank_pos), 10.0 ** rng.uniform(-5, -2, n_flank),
                      rng.uniform(0.02, 0.5, n_flank))
    gwas = pd.concat([central, flank]).sort_values('position').reset_index(drop=True)
    keep = rng.random(gwas.shape[0]) > drop_frac
    gtex = gwas.loc[keep].reset_index(drop=True)
    return gwas, gtex


def test_the_full_power_path_produces_a_well_formed_trait_set():
    """Mirrors the __main__ wiring: draw -> share -> top up -> add flanks."""
    gwas, gtex = _panels()
    n = 50

    pool = stage2.causal_eligible(gwas, LO, HI, 0)
    chosen, diag = _draw(pool, n=n, seed=42)
    shared, topup = stage2.select_gtex_topup(gtex, chosen, LO, HI, n, seed=43)
    flank = stage2.select_flank_gtex(gtex, 0.01, LO, HI, n, seed=44)
    gtex_traits = pd.concat([shared, topup, flank])
    partners = stage2.trait_partner_table(chosen, shared)

    assert chosen.shape[0] == n
    assert shared.shape[0] + topup.shape[0] == n         # central GTEx always fills
    assert flank.shape[0] == n
    assert gtex_traits.shape[0] == 2 * n

    # Every GTEx trait must be resolvable in the GTEx panel and uniquely named.
    assert set(gtex_traits['position']).issubset(set(gtex['position']))
    assert gtex_traits['position'].duplicated().sum() == 0

    # The partner table accounts for every GWAS locus and agrees with `shared`.
    assert partners.shape[0] == n
    assert int(partners['shared'].sum()) == shared.shape[0]
    assert set(partners.loc[partners['shared'], 'position']) == {
        int(p) for p in shared['position']}

    # The draw did its job.
    assert diag.loc[diag['selected'], 'power'].median() > 10 * diag['power'].median()


def test_the_uniform_path_still_pairs_every_locus():
    """The invariant the scorers rely on, asserted rather than assumed."""
    gwas, gtex = _panels()
    gtex_maf = dict(zip(gtex['position'], gtex['maf']))
    eligible = gwas.loc[(gwas['maf'] >= 0.01) & (gwas['selco'] != 0)
                        & (gwas['position'] > LO) & (gwas['position'] < HI)
                        & (gwas['position'].isin(gtex['position']))]
    eligible = eligible[[gtex_maf[p] >= 0.01 for p in eligible['position']]]
    chosen = stage2.subsample_traits(eligible, 50, seed=1)
    partners = stage2.trait_partner_table(chosen, chosen)
    assert partners['shared'].all()
    assert set(chosen['position']).issubset(set(gtex['position']))


def test_select_flank_gtex_takes_the_edges_and_applies_the_causal_floor():
    v = make_vars([3e5, 4e5, 9e5, 1.7e6, 1.8e6],
                  [1e-3] * 5, [0.3, 0.001, 0.3, 0.3, 0.001])
    flank = stage2.select_flank_gtex(v, 0.01, LO, HI, 10, seed=1)
    assert sorted(flank["position"]) == [3e5, 1.7e6]


# ------------------------------------------- the phenotype loop, and its two frames
#
# combine_phenos_to_df is handed `positions` and `vars`, and under power sampling
# with neutral_trait_vars those are frames from DIFFERENT panels: select_gtex_topup
# returns the shared central loci as donor rows out of the GWAS frame, while `vars`
# is the GTEx one. The donor need not segregate in the GTEx panel -- measured on the
# 30x psamp arm, 5316 of 7012 central donors do not -- so anything resolved against
# `vars` by donor position is a KeyError waiting to happen. Only the recipient may
# be looked up there; the donor's selco rides in on its own row.

class _FakeTS:
    """Just enough tree sequence for combine_phenos_to_df's IID column."""

    def __init__(self, n=4):
        self._inds = [types.SimpleNamespace(id=i) for i in range(n)]

    def individuals(self):
        return self._inds


def _record_calls(monkeypatch, n_ind=4):
    """Replace generate_phenos_from_pos with a recorder.

    The real one calls tstrait, which is a stub module in these tests. What is
    under test is the plumbing above it: which donor, which recipient and which
    selection coefficient each trait is built from. The stub still resolves the
    site id through `site_by_pos`, so a recipient missing from the panel fails
    here exactly as it would in the real function.
    """
    calls = []

    def fake(position, ts, vars, scaling=1, recipient_position=None,
             flip_seed=False, seed=0, site_by_pos=None, selco=None):
        if recipient_position is None:
            recipient_position = position
        site_id = int(site_by_pos[recipient_position])
        calls.append({"position": position, "recipient": recipient_position,
                      "selco": selco, "site_id": site_id, "scaling": scaling})
        tr_id = "tr" + str(int(recipient_position))
        return types.SimpleNamespace(
            trait=pd.DataFrame({"position": [float(recipient_position)],
                                "site_id": [site_id],
                                "effect_size": [float(selco)],
                                "trait_id": [tr_id]}),
            phenotype=pd.DataFrame({"phenotype": np.zeros(n_ind)}))

    monkeypatch.setattr(stage2, "generate_phenos_from_pos", fake)
    return calls


def _redistributed_panels():
    """Donors that segregate only in the GWAS panel; recipients in both.

    The shape power sampling produces under neutral_trait_vars: the donor pool is
    not intersected with the GTEx panel, so the donors here are rare and absent
    from it, while select_gtex_topup has already guaranteed the recipients are
    present. 2e5 is a flank locus -- a GTEx-only trait, never a redistribution
    donor, which must keep its own effect.
    """
    gwas = make_vars([6e5, 7e5, 9e5, 1.0e6],
                     [1e-3, 4e-3, 0.0, 0.0],
                     [0.002, 0.003, 0.3, 0.4])
    gtex = make_vars([2e5, 9e5, 1.0e6], [5e-3, 0.0, 0.0], [0.25, 0.3, 0.4])
    donors = gwas.loc[gwas["selco"] != 0]
    return gwas, gtex, donors, {600000: 9e5, 700000: 1.0e6}


def test_the_gtex_loop_reads_the_donor_selco_off_the_row_not_the_panel(monkeypatch):
    """The B x power regression: the donor is not in the GTEx frame at all.

    This used to raise KeyError: <donor position> from selco_by_pos[position].
    """
    gwas, gtex, donors, redist = _redistributed_panels()
    assert not set(donors["position"]) & set(gtex["position"])   # the crash shape

    calls = _record_calls(monkeypatch)
    stage2.combine_phenos_to_df(donors, _FakeTS(), gtex, scaling=30,
                                redistribution=redist, flip_seed=True, seed=11)
    assert [c["selco"] for c in calls] == [1e-3, 4e-3]


def test_both_panels_give_a_shared_locus_the_same_effect_size(monkeypatch):
    """selco is mutation metadata, so the GWAS and GTEx traits must agree on beta."""
    gwas, gtex, donors, redist = _redistributed_panels()

    calls = _record_calls(monkeypatch)
    stage2.combine_phenos_to_df(donors, _FakeTS(), gwas, scaling=30,
                                redistribution=redist, seed=11)
    from_gwas = [c["selco"] for c in calls]

    calls = _record_calls(monkeypatch)
    stage2.combine_phenos_to_df(donors, _FakeTS(), gtex, scaling=30,
                                redistribution=redist, flip_seed=True, seed=11)
    assert [c["selco"] for c in calls] == from_gwas


def test_the_effect_lands_on_the_recipient_when_redistributing(monkeypatch):
    gwas, gtex, donors, redist = _redistributed_panels()
    site_of = dict(zip(gtex["position"], gtex["id"]))

    calls = _record_calls(monkeypatch)
    key, phenos = stage2.combine_phenos_to_df(donors, _FakeTS(), gtex, scaling=30,
                                              redistribution=redist, flip_seed=True,
                                              seed=11)
    assert [c["position"] for c in calls] == [6e5, 7e5]          # donors
    assert [c["recipient"] for c in calls] == [9e5, 1.0e6]       # recipients
    assert [c["site_id"] for c in calls] == [site_of[9e5], site_of[1.0e6]]
    # Traits are named for the recipient, and that is what reaches the sbams header.
    assert list(key["trait_id"]) == ["tr900000", "tr1000000"]
    assert [c for c in phenos.columns if c not in ("FID", "IID")] == \
        ["tr900000", "tr1000000"]


def test_a_position_that_is_not_a_donor_keeps_its_own_effect(monkeypatch):
    """Flank and top-up rows travel through the same loop under a live map."""
    gwas, gtex, donors, redist = _redistributed_panels()
    flank = gtex.loc[gtex["position"] == 2e5]
    positions = pd.concat([donors, flank])

    calls = _record_calls(monkeypatch)
    stage2.combine_phenos_to_df(positions, _FakeTS(), gtex, scaling=30,
                                redistribution=redist, flip_seed=True, seed=11)
    assert calls[-1]["position"] == 2e5 and calls[-1]["recipient"] == 2e5
    assert calls[-1]["selco"] == 5e-3


def test_a_trait_position_missing_from_the_panel_fails_loudly(monkeypatch):
    """The recipient IS resolved against `vars`, so its absence is a real error.

    select_gtex_topup guarantees it never happens; if a caller ever breaks that,
    say so rather than raising KeyError from three frames down.
    """
    gwas, gtex, donors, redist = _redistributed_panels()
    short = gtex.loc[gtex["position"] != 1.0e6]                  # drop one recipient

    _record_calls(monkeypatch)
    with pytest.raises(SystemExit, match="not in this panel"):
        stage2.combine_phenos_to_df(donors, _FakeTS(), short, scaling=30,
                                    redistribution=redist, flip_seed=True, seed=11)


def test_the_loop_is_unchanged_without_redistribution(monkeypatch):
    """Categories A/E/F/G: the trait is its own donor and lives in both frames."""
    gtex = make_vars([9e5, 1.0e6], [1e-3, 4e-3], [0.3, 0.4])

    calls = _record_calls(monkeypatch)
    stage2.combine_phenos_to_df(gtex, _FakeTS(), gtex, scaling=30, seed=11)
    assert [c["position"] for c in calls] == [c["recipient"] for c in calls]
    assert [c["selco"] for c in calls] == [1e-3, 4e-3]


# ------------------------------------------------------- the redistribution map

def test_build_redistribution_map_pairs_one_to_one():
    donors = [6e5, 7e5, 8e5]
    recipients = [9e5, 1.0e6, 1.1e6, 1.2e6]
    kept, redist = stage2.build_redistribution_map(donors, recipients, seed=1)
    assert kept == [600000, 700000, 800000]                     # int keys
    assert sorted(redist) == kept
    assert len(set(redist.values())) == len(redist)             # distinct recipients
    assert set(redist.values()) <= set(recipients)


def test_build_redistribution_map_truncates_donors_when_recipients_run_out():
    donors = [6e5, 7e5, 8e5, 9e5]
    kept, redist = stage2.build_redistribution_map(donors, [1.0e6, 1.1e6], seed=1)
    assert len(kept) == 2 and len(redist) == 2
    assert kept == sorted(kept)                                  # position order held
    assert set(kept) <= {600000, 700000, 800000, 900000}


def test_build_redistribution_map_is_reproducible_from_the_seed():
    donors, recipients = [6e5, 7e5, 8e5], [9e5, 1.0e6, 1.1e6, 1.2e6]
    a = stage2.build_redistribution_map(donors, recipients, seed=4)
    b = stage2.build_redistribution_map(donors, recipients, seed=4)
    c = stage2.build_redistribution_map(donors, recipients, seed=5)
    assert a == b
    assert a[1] != c[1]
