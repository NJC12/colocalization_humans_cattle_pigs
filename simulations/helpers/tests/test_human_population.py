"""Invariants of the sampled-population key that category J introduced.

Category J is category A with one parameter changed: the AFR branch of
OutOfAfrica_2T12 is sampled instead of the EUR one. That is a small enough
difference to be dangerous. The traps below all fail late and quietly:

  * a EUR tree sequence satisfying an AFR run's stage-1 input. Adoption via
    ``stage1_search_dirs`` matches on filename plus embedded seed, and an A run
    is otherwise indistinguishable from a J one -- same script, same
    ``pipeline``, same ``species``, same downstream stages. If both wrote
    ``hts_{seed}.ts``, the only thing standing between them would be a seed
    band, and the tree sequence would be adopted rather than simulated;
  * the reverse: a ``population``-aware filename that does not fall back to the
    bare ``hts_{seed}.ts``, which would orphan every human tree sequence written
    before the key existed;
  * ``paths.py`` and ``human_simulation_o2.py`` disagreeing about the name, so
    the rule copies a file that is not there;
  * a dot in the stage-1 filename, which ``stage2_run_tag`` would carry into a
    path that ``summarize_coloc.read_params`` then prefix-globs.
"""

import glob
import os
import re

import pytest
import yaml

from helpers import human_demography as hd
from helpers import paths


SIM_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CONFIG_DIR = os.path.join(SIM_DIR, "config")

ALL_CONFIGS = sorted(glob.glob(os.path.join(CONFIG_DIR, "*.yaml")))
HUMAN_AFR = sorted(glob.glob(os.path.join(CONFIG_DIR, "human_afr_*.yaml")))


def load(path):
    return yaml.safe_load(open(path))


def human_configs():
    return [p for p in ALL_CONFIGS if load(p).get("pipeline") == "human"]


# ------------------------------------------------------------ the naming rule

def test_eur_keeps_the_historical_filename():
    """The fallback that keeps every pre-category-J tree sequence adoptable."""
    assert hd.stage1_ts_name("EUR", 11) == "hts_11.ts"


def test_afr_gets_its_own_filename():
    assert hd.stage1_ts_name("AFR", 101) == "hts_AFR_101.ts"


def test_the_two_populations_never_share_a_filename():
    """The whole point. Same seed, different population -> different file."""
    for seed in (11, 101, 55):
        assert hd.stage1_ts_name("AFR", seed) != hd.stage1_ts_name("EUR", seed)


def test_an_unknown_population_is_refused_rather_than_passed_through():
    """Otherwise it surfaces as a KeyError inside stdpopsim, in a submitted job."""
    for bad in ("YRI", "eur", "", None):
        with pytest.raises(ValueError):
            hd.stage1_ts_name(bad, 11)
    with pytest.raises(ValueError):
        hd.validate_population("CEU")


def test_the_sample_dict_names_every_population():
    """stdpopsim is handed both populations whichever one is sampled, matching
    the ``{"AFR": 0, "EUR": n}`` shape the script used before this was a knob."""
    d = hd.sample_dict("AFR", 9000)
    assert d == {"AFR": 9000, "EUR": 0}
    assert set(d) == set(hd.POPULATIONS)
    assert hd.sample_dict("EUR", 9000) == {"AFR": 0, "EUR": 9000}


# --------------------------------------------------------- paths.py agreement

@pytest.mark.parametrize("path", human_configs(), ids=os.path.basename)
def test_every_human_config_resolves_through_the_naming_rule(path):
    cfg = load(path)
    # The model has to be resolved before the population: since categories M/N
    # the default population is the MODEL's default, not a global constant.
    model = cfg.get("demographic_model") or hd.DEFAULT_MODEL
    expected = hd.stage1_ts_name(
        cfg.get("population") or hd.default_population(model),
        cfg["stage1_seed"], model)
    assert paths.stage1_human_ts(cfg) == expected
    assert paths.stage1_full_filename(cfg) == expected


def test_a_config_without_the_key_still_means_eur():
    """No config written before category J carries ``population``, and
    helpers/tests feeds raw YAML that has not been through the Snakefile's
    setdefault -- so paths.py must not require the key."""
    assert paths.stage1_human_ts({"stage1_seed": 11}) == "hts_11.ts"


def test_the_stage1_filename_has_no_extra_dots():
    """stage2_run_tag strips two extensions and then appends dot-joined
    segments; summarize_coloc.read_params prefix-globs the result."""
    name = hd.stage1_ts_name("AFR", 101)
    assert name.count(".") == 1
    tag = paths.stage2_run_tag({"pipeline": "human", "population": "AFR",
                                "stage1_seed": 101, "causal_min_maf": 0.01})
    assert tag == "hts_AFR_101"


def test_the_script_and_paths_agree_on_the_filename():
    """human_simulation_o2.py must not spell the name out for itself. Read the
    source rather than importing it -- it parses argv and imports stdpopsim at
    module level."""
    src = open(os.path.join(SIM_DIR, "human_simulation_o2.py")).read()
    assert "human_demography.stage1_ts_name(" in src
    assert not re.search(r"hts_\{seed\}\.ts", src), \
        "the filename is built in helpers/human_demography.py, not here"


# ------------------------------------------------------- the category J config

def test_the_category_j_config_exists():
    """A missing cell is a config_of() miss in four submit scripts, which fails
    loudly -- but only once someone runs that cell."""
    assert HUMAN_AFR, "expected at least config/human_afr_2Mb_g5t20_r3.yaml"


@pytest.mark.parametrize("path", HUMAN_AFR, ids=os.path.basename)
def test_category_j_is_category_a_with_one_parameter_changed(path):
    """The identity IS the experiment. If these drift apart, A - J stops being
    the ancestry contrast and becomes an uninterpretable mixture."""
    j = load(path)
    a = load(path.replace("human_afr_", "human_"))

    assert j["population"] == "AFR"
    assert a.get("population", hd.DEFAULT_POPULATION) == "EUR"
    assert j["pipeline"] == a["pipeline"] == "human"
    assert j["species"] == a["species"] == "human"

    # Everything except the population, the name and the seeds must match.
    differ = {"basename", "population", "stage1_seed", "stage2_seed",
              "workdir", "publishdir"}
    for k in sorted(set(a) | set(j)):
        if k in differ:
            continue
        assert a.get(k) == j.get(k), f"{os.path.basename(path)} diverges from A on {k!r}"


@pytest.mark.parametrize("path", HUMAN_AFR, ids=os.path.basename)
def test_the_seed_is_on_the_j_band(path):
    """seed = 10*letter_index + replicate, and J is the tenth letter, so 101-109.
    Disjoint from the 11-99 that A through I occupy."""
    assert 101 <= int(load(path)["stage1_seed"]) <= 109


@pytest.mark.parametrize("path", HUMAN_AFR, ids=os.path.basename)
def test_category_j_builds_its_own_stage1(path):
    """There is no AFR predecessor to adopt, and a EUR stage1 dir cannot supply
    one. A populated stage1_search_dirs here would look like a dependency this
    arm does not have."""
    assert load(path).get("stage1_search_dirs") in ([], None)
