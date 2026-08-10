"""Tests for the cattle epoch schedule category I's stage 1 runs backwards.

The schedule is a pair of literal vectors that exist in three places now: two
SLiM scripts (farm_selection.slim, farm_selection_from_ep8.slim) and
helpers/cattle_demography.py. SLiM cannot import Python, so the copies can drift,
and a drifted copy fails silently -- I would simply be a different demography
from E/F/G while claiming to be the same one.

What makes that catchable is that the round-3 cattle configs independently state
one derived quantity: `handoff_ticks: 2400`, the ticks spanned by epochs 8-12 at
Q_scaling 0.01. Any transcription slip in the last five epochs moves it. The
tests below pin that, the total, and the backwards-time boundaries the coalescent
actually consumes.

No heavy imports here on purpose -- helpers/cattle_demography.py is deliberately
msprime-free so this can run in the same stubbed environment as the other tests.
"""

import glob
import os

import pytest
import yaml

from helpers import cattle_demography as cd


CONFIG_DIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))), "config")

Q = 0.01  # the value every round-3 cattle config runs at


# --------------------------------------------------------------- the table

def test_the_two_vectors_line_up_and_are_padded_the_way_slim_indexes_them():
    """Index 0 is dummy padding so epoch N sits at index N, matching SLiM."""
    assert len(cd.EP_SIZES) == len(cd.EP_LENGTHS) == cd.LAST_EPOCH + 1
    assert cd.EP_SIZES[0] is None and cd.EP_LENGTHS[0] is None
    assert all(cd.EP_SIZES[k] > 0 for k in range(1, cd.LAST_EPOCH + 1))


def test_the_population_only_ever_shrinks():
    """Holstein-Friesian: a monotone decline from 17,000 to 90. Not a property
    the code relies on, but an inverted pair is exactly what a transcription
    slip looks like."""
    sizes = [cd.EP_SIZES[k] for k in range(1, cd.LAST_EPOCH + 1)]
    assert sizes == sorted(sizes, reverse=True)
    assert sizes[0] == 17000 and sizes[-1] == 90


# ------------------------------------------------- the cross-check with SLiM

def test_epochs_8_to_12_come_to_the_handoff_ticks_the_configs_declare():
    """2400 = 600+600+600+300+300. E/F/G's tree sequences span the Q=1 ->
    Q=0.01 boundary at exactly this depth, and their configs say so in a key
    this module has no other connection to."""
    assert cd.epochs_8_to_12_ticks(Q) == 2400
    assert cd.epochs_8_to_12_ticks(1.0) == 24        # 24 real generations


@pytest.mark.parametrize("path", sorted(glob.glob(os.path.join(
    CONFIG_DIR, "cattle_*_from_midpoint_2Mb*_r3.yaml"))
    + glob.glob(os.path.join(CONFIG_DIR, "cattle_sel_*_2Mb*_r3.yaml"))))
def test_every_split_q_cattle_config_agrees_with_the_table(path):
    """The live assertion: read handoff_ticks straight out of the configs that
    declare it and require the table to reproduce it at their own Q_scaling."""
    cfg = yaml.safe_load(open(path))
    if "handoff_ticks" not in cfg:
        pytest.skip(f"{os.path.basename(path)} declares no handoff_ticks")
    assert cd.epochs_8_to_12_ticks(float(cfg["Q_scaling"])) == cfg["handoff_ticks"]


def test_the_total_is_the_published_post_burn_in_length():
    """3,354 real generations from the start of epoch 2 to sampling. Epoch 1 is
    excluded: backwards in time it is the ancestral state and has no end."""
    assert cd.total_ticks(1.0) == 3354
    assert cd.total_ticks(Q) == 335400


# ------------------------------------------- what the coalescent consumes

def test_size_changes_walk_backwards_from_epoch_11_to_epoch_1():
    changes = cd.size_changes(Q)
    assert [epoch for _, _, epoch in changes] == list(range(11, 0, -1))
    times = [t for t, _, _ in changes]
    assert times == sorted(times), "times must increase going back"
    assert times[0] > 0, "the terminal size is the initial state, not a change"


def test_the_boundaries_are_the_ones_the_stage_1_log_prints():
    """Pinned literally: these are the numbers a reader checks the log against,
    and the handoff at 2400 is visible among them."""
    assert [t for t, _, _ in cd.size_changes(Q)] == [
        300, 600, 1200, 1800, 2400, 15400, 45400, 65400, 175400, 235400, 335400]


def test_the_boundaries_are_integers_not_floating_point_debris():
    """`3 / 0.01` is one of the divisions that can land at 299.99999999999994.
    A boundary off by an ULP is harmless to msprime and poisonous to a reader
    diffing the log against the SLiM script."""
    for t, size, _ in cd.size_changes(Q):
        assert t == int(t), t
        assert size == int(size), size


def test_the_terminal_population_is_exactly_the_sample_the_cattle_arms_take():
    """9,000 = gwas_size 8,000 + the largest GTEx panel 1,000, which is also
    E/F's whole SLiM population. rules/stage1_cattle_neutral.smk refuses to
    over-sample past this."""
    assert cd.terminal_size(Q) == 9000
    assert cd.ancestral_size(Q) == 1700000


def test_the_schedule_is_a_pure_rescaling_in_Q():
    """Every size scales by 1/Q and every duration by 1/Q, which is what makes
    the Hudson coalescent give the identical tree sequence at any Q. If this
    ever stops holding, the Q-invariance check in the stage-1 log stops meaning
    anything."""
    fine, coarse = cd.size_changes(0.01), cd.size_changes(1.0)
    for (t_f, n_f, e_f), (t_c, n_c, e_c) in zip(fine, coarse):
        assert e_f == e_c
        assert t_f == pytest.approx(t_c * 100)
        assert n_f == pytest.approx(n_c * 100)


def test_describe_names_every_epoch_once():
    text = cd.describe(Q)
    for k in range(cd.FIRST_EPOCH, cd.LAST_EPOCH + 1):
        assert f"epoch  {k:>2}" in text
    assert "9,000" in text and "1,700,000" in text
