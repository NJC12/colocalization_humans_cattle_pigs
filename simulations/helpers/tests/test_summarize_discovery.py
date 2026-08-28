"""`summarize_coloc.runs_from_manifest` against the arm-major publication tree.

The publication layout is `<root>/<arm>/<run_dir>`, RUNS.tsv sits at `<root>`, and
`--roots` is given the ARM directories. The run_dir names are arm-relative and the
same 30 repeat under every arm, so a manifest row for one arm resolves to a real
directory under any other arm. Discovery has to filter on the arm or every run is
counted once per arm and labelled with all of them.
"""

import os

from helpers import summarize_coloc as sc


ARMS = ["causal_maf001_paired", "causal_maf001_unpaired",
        "causal_power_n200000", "causal_power_n30000"]
RUNS = ["human_baseline_negative_selection_rep1",
        "cattle_baseline_negative_selection_rep1"]
LETTER = {"human_baseline_negative_selection_rep1": "A",
          "cattle_baseline_negative_selection_rep1": "E"}


def _tree(tmp_path, arms=ARMS, runs=RUNS):
    root = tmp_path / "pub"
    root.mkdir()
    lines = ["arm\trun_dir\tletter"]
    for a in arms:
        for r in runs:
            (root / a / r).mkdir(parents=True)
            lines.append(f"{a}\t{r}\t{LETTER[r]}")
    (root / "RUNS.tsv").write_text("\n".join(lines) + "\n")
    return root


def test_an_arm_root_sees_only_its_own_runs(tmp_path):
    root = _tree(tmp_path)
    for arm in ARMS:
        found = sc.runs_from_manifest(str(root / arm))
        assert found is not None, arm
        assert len(found) == len(RUNS), (arm, found)
        assert {a for _, _, a in found} == {arm}, (arm, found)
        for rep, _, _ in found:
            assert os.path.basename(os.path.dirname(rep)) == arm


def test_every_run_is_visited_exactly_once_across_all_arm_roots(tmp_path):
    """The bug this pins: 4 arms x 2 runs must be 8 visits, not 32."""
    root = _tree(tmp_path)
    seen = [rep for arm in ARMS for rep, _, _ in sc.runs_from_manifest(str(root / arm))]
    assert len(seen) == len(ARMS) * len(RUNS)
    assert len(set(seen)) == len(seen), "a directory was visited more than once"


def test_letters_come_from_the_manifest_not_the_path(tmp_path):
    """The interpretable names carry no letter, which is why RUNS.tsv exists."""
    root = _tree(tmp_path)
    found = sc.runs_from_manifest(str(root / ARMS[0]))
    assert {l for _, l, _ in found} == {"A", "E"}


def test_a_root_that_is_not_an_arm_is_left_unfiltered(tmp_path):
    """A flat single-arm tree keeps working: RUNS.tsv beside the root, run dirs
    directly beneath it, and no arm segment in the path."""
    root = tmp_path / "flat"
    root.mkdir()
    (root / "RUNS.tsv").write_text(
        "arm\trun_dir\tletter\nsome_arm\thuman_baseline_negative_selection_rep1\tA\n")
    (root / "human_baseline_negative_selection_rep1").mkdir()
    found = sc.runs_from_manifest(str(root))
    assert len(found) == 1
    assert found[0][2] == "some_arm"


def test_absent_manifest_returns_none_so_the_glob_fallback_runs(tmp_path):
    """None, not [] -- the caller distinguishes 'no manifest' from 'empty'."""
    d = tmp_path / "nothing"
    d.mkdir()
    assert sc.runs_from_manifest(str(d)) is None


def test_rows_for_runs_that_have_not_landed_are_skipped(tmp_path):
    root = _tree(tmp_path, runs=RUNS)
    import shutil
    shutil.rmtree(root / ARMS[0] / RUNS[1])
    found = sc.runs_from_manifest(str(root / ARMS[0]))
    assert [os.path.basename(r) for r, _, _ in found] == [RUNS[0]]
