"""`helpers/publication_deep_histories.tsv` and the two readers over it.

WHY. Through the 5-replicate publication set every cattle replicate resumed from
ONE deep history (`cattle_baseline_seed: 20250303`): 33,130 of 33,154 generations
shared, only epochs 8-12 simulated per replicate. Cattle "replicates" were
therefore not independent populations, and their spread did not measure what the
human arm's spread measures. This table assigns each replicate a deep history and
is the record of which seeds were used.

Two readers again -- bash for the launcher, Python for the manifest -- so the
drift risk that `test_publication_tables.py` exists to close applies here too.
"""

import os
import subprocess

import pytest

from helpers import paths
from helpers import write_run_manifests as wrm


SIM_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LIB = os.path.join(SIM_DIR, "lib", "publication_common.sh")

DH_HEADER, DH_ROWS = wrm.read_tsv(wrm.DEEP_HISTORIES_TSV)
DH = {int(r["replicate"]): r for r in DH_ROWS}

_, CAT_ROWS = wrm.read_tsv(wrm.CATEGORIES_TSV)
CATTLE = [r["letter"] for r in CAT_ROWS if r["species"] == "cattle"]


def sh(snippet):
    out = subprocess.run(
        ["bash", "-c", f'REPO="{SIM_DIR}"; CELL="g5t20"; source "{LIB}"; {snippet}'],
        capture_output=True, text=True)
    assert out.returncode == 0, f"{snippet}\n{out.stdout}\n{out.stderr}"
    return out.stdout.strip()


# ------------------------------------------------------------------- the table

def test_it_covers_every_replicate_exactly_once():
    assert sorted(DH) == list(range(1, 21))
    assert len(DH_ROWS) == 20


def test_the_first_five_replicates_keep_the_round_three_history():
    """80 finished cattle runs are named for the stage-1 files this history
    produced. paths.LEGACY_CATTLE_BASELINE_SEED is suppressed at this value for
    the same reason, and the two must not drift apart."""
    for rep in (1, 2, 3, 4, 5):
        assert int(DH[rep]["cattle_baseline_seed"]) == paths.LEGACY_CATTLE_BASELINE_SEED


def test_the_blocks_are_contiguous_and_equal_sized():
    """Five replicates per history. Not load-bearing for correctness, but an
    unbalanced design is a statistics problem, so it should be deliberate."""
    from collections import Counter
    counts = Counter(r["cattle_baseline_seed"] for r in DH_ROWS)
    assert sorted(counts.values()) == [5, 5, 5, 5]
    for cb, n in counts.items():
        reps = sorted(int(r["replicate"]) for r in DH_ROWS
                      if r["cattle_baseline_seed"] == cb)
        assert reps == list(range(reps[0], reps[0] + 5)), (cb, reps)


def test_the_deep_history_seeds_are_distinct_and_could_not_be_a_stage1_seed():
    """8-digit YYYYMMDD, so a deep-history seed can never be mistaken for a
    stage-1 seed (2-3 digits), and < 2**32 so msprime accepts it -- and so does
    `seed + 1`, which farm_create_orig_pop_e2.py uses for the mutation overlay."""
    seeds = {int(r["cattle_baseline_seed"]) for r in DH_ROWS}
    assert len(seeds) == 4
    for s in seeds:
        assert 10_000_000 <= s <= 99_999_999
        assert s + 1 < 2 ** 32


def test_the_handoff_filename_is_the_one_paths_py_will_ask_for():
    """The table spells the name; the rule builds it from cattle_baseline_seed.
    If these disagree the lookup misses and Snakemake starts a two-hour SLiM run
    instead of resuming."""
    for rep, row in DH.items():
        want = paths.stage1_cattle_baseline_handoff(
            {"Q_scaling": 0.01, "L": 2000000,
             "cattle_baseline_seed": row["cattle_baseline_seed"]})
        assert row["handoff_file"] == want, rep


# ------------------------------------------------- bash and Python must agree

@pytest.mark.parametrize("rep", range(1, 21))
def test_shell_and_python_agree_on_the_deep_history(rep):
    assert sh(f"cattle_baseline_seed_of {rep}") == DH[rep]["cattle_baseline_seed"]
    assert sh(f"deep_history_handoff_of {rep}") == DH[rep]["handoff_file"]


def test_the_shell_replicate_list_is_the_whole_table():
    assert sh("replicate_list").split() == [str(r) for r in sorted(DH)]


def test_an_unknown_replicate_is_refused_not_invented():
    """Silence here would mean falling back to the config YAML's 20250303, which
    collapses the independent histories back to one without any error."""
    out = subprocess.run(
        ["bash", "-c",
         f'REPO="{SIM_DIR}"; source "{LIB}"; cattle_baseline_seed_of 21'],
        capture_output=True, text=True)
    assert out.returncode != 0
    assert out.stdout.strip() == ""


# --------------------------------------------------------- reaching the manifest

def _manifest(tmp_path, reps=None):
    args = ["--out-dir", str(tmp_path),
            "--scratch-root", "/scratch", "--publish-root", "/publish"]
    if reps is not None:
        args += ["--reps"] + [str(r) for r in reps]
    import sys
    old = sys.argv
    sys.argv = ["write_run_manifests.py"] + args
    try:
        wrm.main()
    finally:
        sys.argv = old
    rows = []
    with open(os.path.join(tmp_path, "RUNS.tsv")) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            rows.append(dict(zip(head, line.rstrip("\n").split("\t"))))
    return rows


def test_the_default_replicate_list_is_the_whole_table(tmp_path):
    """THE TRUNCATION BUG: --reps defaulted to a literal [1,2,3,4,5], and
    publish_to_data2.sh regenerates the manifests with no --reps -- so publishing
    an arm rewrote the published RUNS.tsv back down to five replicates."""
    rows = _manifest(tmp_path)
    assert len(rows) == 4 * 6 * 20
    assert sorted({int(r["replicate"]) for r in rows}) == list(range(1, 21))


def test_the_manifest_carries_the_deep_history(tmp_path):
    rows = _manifest(tmp_path, reps=[1, 6, 20])
    for r in rows:
        if r["species"] == "human":
            assert r["cattle_baseline_seed"] == "-"
            assert r["deep_history_handoff"] == "-"
        else:
            rep = int(r["replicate"])
            assert r["cattle_baseline_seed"] == DH[rep]["cattle_baseline_seed"]
            assert r["deep_history_handoff"] == DH[rep]["handoff_file"]


def test_the_deep_history_is_in_the_cattle_stage1_filenames(tmp_path):
    """The point of the whole change: a cattle stage-1 file must name the
    ancestry it came from, or two histories are indistinguishable on disk."""
    rows = _manifest(tmp_path, reps=[6])
    for r in rows:
        if r["species"] != "cattle":
            continue
        assert f"cb{r['cattle_baseline_seed']}" in r["stage1_file"].replace("cb_", "cb")


def test_cattle_sel_runs_name_their_marks_file(tmp_path):
    """rules/stage1_cattle_sel.smk requires the marks file AND the .full.ts; the
    launcher cannot pre-flight what the manifest does not name."""
    rows = _manifest(tmp_path, reps=[6])
    for r in rows:
        if r["pipeline"] == "cattle_sel":
            assert r["stage1_marks_file"].endswith(".m4_marks.tsv")
            assert (r["stage1_marks_file"][: -len(".m4_marks.tsv")]
                    == r["stage1_file"][: -len(".full.ts")])
        else:
            assert r["stage1_marks_file"] == "-"


def test_every_replicate_of_a_cattle_category_gets_a_distinct_stage1_file(tmp_path):
    rows = _manifest(tmp_path)
    for letter in CATTLE:
        files = [r["stage1_file"] for r in rows
                 if r["letter"] == letter and r["arm"] == "causal_maf001_paired"]
        assert len(set(files)) == len(files) == 20, letter
