"""`cattle_baseline_segment`: the cattle deep history, made visible in the filename.

WHY THIS EXISTS. Every cattle replicate through the 5-replicate publication set
resumed from ONE deep history (`cattle_baseline_seed: 20250303`) -- 33,130 of
33,154 generations shared, with only epochs 8-12 simulated per replicate. Giving
each block of five replicates its own history is what makes cattle replicates
independent populations. But the `cattle_sel` pipeline (categories F and G) named
its stage-1 outputs `{base}_mult{m}_gen{g}_muts{n}_sd{stage1_seed}.full.ts`, with
no reference to the deep history at all: two F runs at one `stage1_seed` off
different ancestries would write the same `.full.ts`, the same `.m4_marks.tsv`,
the same `stage2_run_tag`, and after publishing the same files in the same
directory. `cattle_baseline_seed` is not in `params_record.STAGE2_KEYS` either,
so nothing downstream could separate them.

The segment is suppressed at the legacy seed, so it closes that gap while moving
ZERO existing paths. `test_the_legacy_deep_history_emits_no_filename_segment` is
the assertion that keeps that true: 40 finished F/G runs are named for these
files, and a run whose stage-1 filename moves is a run that has been orphaned.
"""

import pytest

from helpers import paths, search_dirs


F_BASE = dict(basename="cattle_sel_bot", selection_multiplier=100,
              selected_generations=23, num_muts_selected=5,
              Q_scaling=0.01, L=2000000, pipeline="cattle_sel")
G_BASE = dict(F_BASE, basename="cattle_sel_nobot")

LEGACY = paths.LEGACY_CATTLE_BASELINE_SEED   # 20250303
NEW = 20260901


def cfg(base, seed, cb):
    return {**base, "stage1_seed": seed, "cattle_baseline_seed": cb}


# --------------------------------------------------------------- the guarantee

@pytest.mark.parametrize("base,prefix,seed", (
    [(F_BASE, "cattle_sel_bot", s) for s in (61, 62, 63, 64, 65)] +
    [(G_BASE, "cattle_sel_nobot", s) for s in (71, 72, 73, 74, 75)]
))
def test_the_legacy_deep_history_emits_no_filename_segment(base, prefix, seed):
    """Literal names for the F/G replicates that have already run.

    Not a round-trip through the same function that builds them -- the literal
    string is the point. 40 published runs are named for these, and the whole
    argument for adding the segment is that it costs nothing already on disk.
    """
    stem = f"{prefix}_mult100_gen23_muts5_sd{seed}"
    assert paths.stage1_cattle_sel_full(cfg(base, seed, LEGACY)) == f"{stem}.full.ts"
    assert paths.stage1_cattle_sel_marks(cfg(base, seed, LEGACY)) == f"{stem}.m4_marks.tsv"


def test_the_segment_is_empty_at_the_legacy_seed_however_it_is_spelled():
    """The config carries an int; a --config override arrives as a string."""
    for spelling in (20250303, "20250303"):
        assert paths.cattle_baseline_segment({"cattle_baseline_seed": spelling}) == ""


def test_an_absent_key_emits_nothing():
    """Human and the neutral pipelines never set it, and must not gain a segment."""
    assert paths.cattle_baseline_segment({}) == ""
    assert paths.cattle_baseline_segment({"cattle_baseline_seed": None}) == ""


# ------------------------------------------------------------- the new capability

def test_a_new_deep_history_changes_the_cattle_sel_filename():
    a = paths.stage1_cattle_sel_full(cfg(F_BASE, 61, LEGACY))
    b = paths.stage1_cattle_sel_full(cfg(F_BASE, 61, NEW))
    assert a != b
    assert b == f"cattle_sel_bot_mult100_gen23_muts5_sd61_cb{NEW}.full.ts"


def test_the_marks_file_moves_with_the_full_file():
    """`rules/stage1_cattle_sel.smk` requires BOTH; a pair that disagreed would
    make the rule adopt one and simulate the other."""
    full = paths.stage1_cattle_sel_full(cfg(F_BASE, 610, NEW))
    marks = paths.stage1_cattle_sel_marks(cfg(F_BASE, 610, NEW))
    assert full[: -len(".full.ts")] == marks[: -len(".m4_marks.tsv")]


@pytest.mark.parametrize("cb", [LEGACY, NEW])
def test_the_from_midpoint_name_already_carried_the_deep_history(cb):
    """E and L never had this gap. Regression only."""
    name = paths.stage1_cattle_baseline_from_midpoint_full(
        {"Q_scaling": 0.01, "L": 2000000, "cattle_baseline_seed": cb,
         "stage1_seed": 51})
    assert f"cb_{cb}" in name


# ------------------------------------------------- it does not break seed lookup

@pytest.mark.parametrize("seed", [61, 610, 620])
def test_the_seed_extractor_still_reads_the_stage1_seed(seed):
    """The segment sits AFTER `_sd{seed}`, so `_extract_seed`'s first match is
    still the stage-1 seed and not the 8-digit deep-history seed."""
    name = paths.stage1_cattle_sel_full(cfg(F_BASE, seed, NEW))
    assert search_dirs._extract_seed(name) == str(seed)


def test_the_collision_glob_is_scoped_to_one_deep_history():
    """`_seed_collision_check` rewrites `sd{N}` to `sd*`. With the segment after
    it, the glob still ends in the deep history, so a sibling history's file is
    not reported as a seed conflict for this one."""
    name = paths.stage1_cattle_sel_full(cfg(F_BASE, 610, NEW))
    glob_pat = search_dirs._SEED_GLOB_REPLACE.sub("sd*", name)
    assert glob_pat == f"cattle_sel_bot_mult100_gen23_muts5_sd*_cb{NEW}.full.ts"
    # and it does NOT match a file from a different history at another seed
    import fnmatch
    other = paths.stage1_cattle_sel_full(cfg(F_BASE, 615, 20260902))
    assert not fnmatch.fnmatch(other, glob_pat)


def test_two_deep_histories_at_one_stage1_seed_are_now_distinguishable():
    """The failure this segment exists to prevent, stated directly."""
    names = {paths.stage1_cattle_sel_full(cfg(F_BASE, 61, cb))
             for cb in (LEGACY, 20260901, 20260902, 20260903)}
    assert len(names) == 4
