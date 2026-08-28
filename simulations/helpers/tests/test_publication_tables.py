"""The publication category/arm tables, and the two readers that must agree.

`helpers/publication_categories.tsv` and `helpers/publication_arms.tsv` are read
from BOTH bash (`lib/publication_common.sh`, which builds the --config tokens and
the seeds) and Python (`helpers/write_run_manifests.py`, which builds RUNS.tsv and
whose `stage1_file` column is the launcher's only guard on stage-1 adoption for
the human categories). Two readers over one table is the drift risk these tests
close.

The pairing invariant is the one worth stating plainly. K adopts A's tree
sequence and L adopts E's, at the DONOR's seed, so that A-K and E-L hold the
genealogy fixed and isolate the effect assignment. Stage-1 lookup is seed-strict,
so if a paired category's seed ever left its donor's band the adoption would not
fail -- `search_dirs` would simply not find the file, Snakemake would decide it
has to build one, and a multi-hour genetic simulation would start. The pairing
would be silently gone and the contrast would silently become noisy.
"""

import os
import subprocess

import pytest
import yaml

from helpers import write_run_manifests as wrm


SIM_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LIB = os.path.join(SIM_DIR, "lib", "publication_common.sh")

CAT_HEADER, CAT_ROWS = wrm.read_tsv(wrm.CATEGORIES_TSV)
ARM_HEADER, ARM_ROWS = wrm.read_tsv(wrm.ARMS_TSV)
CATS = {r["letter"]: r for r in CAT_ROWS}
ARMS = {r["arm"]: r for r in ARM_ROWS}

# Keys the launcher asserts against the YAML rather than passing as --config, so
# there is one source of truth for them.
ASSERT_KEYS = ["gwas_scaling", "gtex_scaling", "gtex_size", "gwas_size", "L"]


def sh(snippet, cell="g5t20"):
    """Run a snippet against the sourced shell library.

    The assignments are separate statements, not a `VAR=x source file` prefix:
    for a builtin that form does not reliably outlive the command, so CELL would
    read empty inside the functions and ${CELL} would interpolate to nothing.
    """
    out = subprocess.run(
        ["bash", "-c", f'REPO="{SIM_DIR}"; CELL="{cell}"; source "{LIB}"; {snippet}'],
        capture_output=True, text=True)
    assert out.returncode == 0, f"{snippet}\nstdout={out.stdout}\nstderr={out.stderr}"
    return out.stdout.strip()


# ------------------------------------------------------------- the table itself

def test_the_six_publication_categories_are_present():
    assert set(CATS) == {"A", "E", "K", "L", "F", "G"}


def test_the_four_publication_arms_are_present():
    assert set(ARMS) == {"causal_maf001_paired", "causal_maf001_unpaired",
                         "causal_power_n200000", "causal_power_n30000"}


def test_directory_names_are_interpretable_and_unique():
    """The whole point of the rename: no reader should need the letter table."""
    names = [r["dir_name"] for r in CAT_ROWS]
    assert len(set(names)) == len(names)
    for n in names:
        assert n.startswith(("human_", "cattle_")), n
        assert n == n.lower() and " " not in n, n


@pytest.mark.parametrize("letter", sorted(CATS))
def test_every_config_exists(letter):
    cfg = CATS[letter]["config"].replace("${CELL}", "g5t20")
    assert os.path.exists(os.path.join(SIM_DIR, cfg)), cfg


@pytest.mark.parametrize("letter", sorted(CATS))
def test_the_pipeline_column_matches_the_config(letter):
    cfg = CATS[letter]["config"].replace("${CELL}", "g5t20")
    with open(os.path.join(SIM_DIR, cfg)) as fh:
        loaded = yaml.safe_load(fh)
    assert loaded["pipeline"] == CATS[letter]["pipeline"]
    assert loaded["species"] == CATS[letter]["species"]


# ---------------------------------------------------------- the pairing invariant

@pytest.mark.parametrize("letter", sorted(CATS))
def test_a_paired_category_sits_in_its_donors_seed_band(letter):
    donor = CATS[letter]["stage1_donor"]
    if donor == "self":
        return
    assert donor in CATS, donor
    assert CATS[letter]["seed_prefix"] == CATS[donor]["seed_prefix"], (
        f"{letter} adopts {donor}'s stage 1 but sits in a different seed band; "
        f"the seed-strict lookup would not adopt, it would re-simulate")


def test_k_and_l_are_paired_to_a_and_e():
    """Stated explicitly rather than only as an invariant, because this is the
    design decision the background-selection contrast rests on."""
    assert CATS["K"]["stage1_donor"] == "A"
    assert CATS["L"]["stage1_donor"] == "E"
    assert CATS["K"]["seed_prefix"] == CATS["A"]["seed_prefix"] == "1"
    assert CATS["L"]["seed_prefix"] == CATS["E"]["seed_prefix"] == "5"


def test_a_shared_config_category_carries_its_distinguishing_override():
    """K and L use A's and E's config files. Without the override they would be A
    and E, written under K's and L's names, with nothing downstream looking wrong."""
    by_config = {}
    for letter, r in CATS.items():
        by_config.setdefault(r["config"], []).append(letter)
    for cfg, letters in by_config.items():
        if len(letters) == 1:
            continue
        for letter in letters:
            extra = CATS[letter]["extra_config"]
            donor = CATS[letter]["stage1_donor"]
            if donor == "self":
                continue  # the config's own category
            assert extra and extra != "-", (
                f"{letter} shares {cfg} with {letters} but sets no override")
            assert "synthetic_dfe_effects=True" in extra, (letter, extra)


def test_no_config_token_contains_a_space():
    """EXTRA_CONFIG is expanded UNQUOTED on the snakemake command line so each
    key=value becomes its own --config arg."""
    for r in CAT_ROWS:
        e = r["extra_config"]
        if e and e != "-":
            for tok in e.split():
                assert "=" in tok and not tok.endswith("="), (r["letter"], tok)


# ------------------------------------------------- bash and Python must agree

@pytest.mark.parametrize("letter", sorted(CATS))
@pytest.mark.parametrize("rep", [1, 5])
def test_shell_and_python_derive_the_same_seed(letter, rep):
    assert sh(f"seed_of {letter}{rep}") == f"{CATS[letter]['seed_prefix']}{rep}"


@pytest.mark.parametrize("letter", sorted(CATS))
def test_shell_and_python_resolve_the_same_config(letter):
    assert sh(f"config_of {letter}") == CATS[letter]["config"].replace("${CELL}", "g5t20")


@pytest.mark.parametrize("letter", sorted(CATS))
def test_shell_and_python_agree_on_dir_name_and_donor(letter):
    assert sh(f"category_dir_of {letter}") == CATS[letter]["dir_name"]
    donor = CATS[letter]["stage1_donor"]
    assert sh(f"stage1_donor_of {letter}") == (letter if donor == "self" else donor)


@pytest.mark.parametrize("letter", sorted(CATS))
def test_shell_and_python_agree_on_the_category_override(letter):
    want = CATS[letter]["extra_config"]
    assert sh(f"category_extra {letter}") == ("" if want == "-" else want)


def test_the_cell_placeholder_is_honoured_by_the_shell():
    """A future cell must not silently resolve to g5t20."""
    assert "x20" in sh("config_of A", cell="x20")


# --------------------------------------------------- the arms, and the YAML split

@pytest.mark.parametrize("arm", sorted(ARMS))
def test_asserted_keys_match_every_config(arm):
    """These are recorded in ARMS.tsv but NOT passed as --config, so the YAML is
    the single source of truth and a drifted YAML must be caught before launch."""
    for letter, c in CATS.items():
        cfg = c["config"].replace("${CELL}", "g5t20")
        with open(os.path.join(SIM_DIR, cfg)) as fh:
            loaded = yaml.safe_load(fh)
        for k in ASSERT_KEYS:
            assert float(ARMS[arm][k]) == float(loaded[k]), (arm, letter, k)


def test_the_arms_differ_only_in_how_causal_loci_are_chosen():
    """If they differ in anything else the design is not a one-factor comparison."""
    varying = {k for k in ARM_HEADER
               if k not in ("arm", "description")
               and len({ARMS[a][k] for a in ARMS}) > 1}
    assert varying == {"causal_min_maf", "causal_sampling", "sampling_gwas_n",
                       "require_gtex_partner"}, varying


def test_the_paired_arm_is_the_only_one_that_requires_a_partner():
    for arm, r in ARMS.items():
        want = arm == "causal_maf001_paired"
        assert (r["require_gtex_partner"] == "True") is want, arm


def test_the_two_uniform_arms_differ_in_exactly_one_key():
    a, b = ARMS["causal_maf001_paired"], ARMS["causal_maf001_unpaired"]
    diff = {k for k in ARM_HEADER if k not in ("arm", "description") and a[k] != b[k]}
    assert diff == {"require_gtex_partner"}, diff


def test_the_two_power_arms_differ_in_exactly_one_key():
    a, b = ARMS["causal_power_n200000"], ARMS["causal_power_n30000"]
    diff = {k for k in ARM_HEADER if k not in ("arm", "description") and a[k] != b[k]}
    assert diff == {"sampling_gwas_n"}, diff


def test_power_arms_require_an_explicit_trait_count():
    """causal_sampling: power refuses to run without one -- at causal_min_maf 0
    the pool is every polymorphic central variant."""
    for arm, r in ARMS.items():
        if r["causal_sampling"] == "power":
            assert int(r["n_central_traits"]) > 0, arm


# The power guard is a launch-time cliff, not a warning. Encoding the measurement
# here is the only place it is cheap to check: the guard itself lives inside
# stage 2, so a too-small sampling_gwas_n is otherwise discovered by thirty
# controllers that each run for half an hour and then refuse.
POWER_N_FLOOR = 20000

def test_power_arms_clear_the_eligible_pool_guard():
    """No power arm may sit below the measured feasibility floor.

    `select_central_power` refuses the draw unless the pool holds
    `sampling_min_pool_multiple * n_central_traits` candidates at
    power >= `sampling_min_power`. The arm was first written at
    sampling_gwas_n=8000 and stage 2 refused three of the five human baseline
    replicates -- 39, 49 and 39 eligible against a required 50.

    Swept over the per-candidate (maf_for_power, beta) diagnostics that the
    n=200000 arm recorded for all 30 runs (the pool does not depend on
    sampling_gwas_n, only the weight does), the number of runs the guard would
    refuse is:

        n= 8,000  3 of 30      n= 40,000  0 of 30
        n=20,000  0 of 30      n= 50,000  0 of 30
        n=30,000  0 of 30      n=200,000  0 of 30

    so 20,000 is the smallest round value at which every run clears it, and the
    published arm uses 30,000, where the tightest run
    (human_baseline_negative_selection_rep2) has 99 eligible against 50.
    """
    for arm, r in ARMS.items():
        if r["causal_sampling"] != "power":
            continue
        n = int(float(r["sampling_gwas_n"]))
        assert n >= POWER_N_FLOOR, (
            f"arm {arm} draws causal loci by power at sampling_gwas_n={n}, below the "
            f"measured floor of {POWER_N_FLOOR}. At n=8000 the human pool supplied "
            f"39-49 candidates at power>=0.05 against the {r['sampling_min_pool_multiple']}"
            f"x{r['n_central_traits']} the guard requires, and stage 2 exits. Raise n, or "
            f"lower sampling_min_power / sampling_min_pool_multiple deliberately and "
            f"move this floor with a fresh sweep."
        )


def test_the_guard_requirement_is_what_the_floor_was_measured_against():
    """The floor above is only valid for the guard settings it was swept at."""
    for arm, r in ARMS.items():
        if r["causal_sampling"] != "power":
            continue
        assert float(r["sampling_min_power"]) == 0.05, arm
        assert float(r["sampling_min_pool_multiple"]) == 2, arm
        assert float(r["sampling_sig_p"]) == 5e-8, arm
        assert int(r["n_central_traits"]) == 25, arm
