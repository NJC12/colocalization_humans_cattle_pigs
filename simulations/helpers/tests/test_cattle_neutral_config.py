"""Invariants of the category-I configs that nothing else would catch in time.

Every trap below fails late and quietly rather than at parse time:

  * a `mutation_rate` copied from category H (1.4e-8) would run this arm at 2.5x
    E/F/G's variant density, because stage 2's cattle branch overlays 8.4e-9 on
    top of whatever stage 1 produced;
  * a stray `handoff_ticks` would make get_vars_df un-rescale the times
    piecewise across a boundary this arm's single-Q tree sequence does not have;
  * a `cattle_baseline_search_dirs` would point at a checkpoint this pipeline
    never reads, and read as a dependency it does not have;
  * `gwas_size + gtex_size` above the terminal population would have Hudson
    over-sample a population of 9,000 -- legal, and no longer an approximation
    of the forward arms.
"""

import glob
import os

import pytest
import yaml

from helpers import cattle_demography as cd
from helpers import paths


SIM_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CONFIG_DIR = os.path.join(SIM_DIR, "config")

NEUTRAL_CATTLE = sorted(glob.glob(os.path.join(CONFIG_DIR, "cattle_neutral_*.yaml")))
ALL_CONFIGS = sorted(glob.glob(os.path.join(CONFIG_DIR, "*.yaml")))


def load(path):
    return yaml.safe_load(open(path))


def test_the_configs_exist_at_all():
    """A missing cell is a config_of() miss in four submit scripts, which fails
    loudly -- but only once someone runs that cell."""
    assert len(NEUTRAL_CATTLE) >= 4, [os.path.basename(p) for p in NEUTRAL_CATTLE]


@pytest.mark.parametrize("path", ALL_CONFIGS, ids=os.path.basename)
def test_every_config_resolves_to_a_stage1_filename(path):
    """stage1_full_filename raises on an unknown pipeline. Running it over every
    config is what catches a new pipeline wired into the Snakefile but not into
    paths.py -- the failure mode is a ValueError deep in a submitted job."""
    cfg = load(path)
    assert paths.stage1_full_filename(cfg)


@pytest.mark.parametrize("path", NEUTRAL_CATTLE, ids=os.path.basename)
def test_the_stage1_filename_is_distinct_from_every_other_arm(path):
    """cnts_ / nts_ / hts_ must not collide: stage-1 adoption matches on
    filename plus embedded seed, so a shared prefix lets one arm's tree sequence
    silently satisfy another's input."""
    cfg = load(path)
    name = paths.stage1_full_filename(cfg)
    assert name == f"cnts_{cfg['stage1_seed']}.ts"
    assert not name.startswith(("hts_", "nts_"))
    assert paths.stage1_marks_filename(cfg) is None   # m4 marks are cattle_sel only


@pytest.mark.parametrize("path", NEUTRAL_CATTLE, ids=os.path.basename)
def test_the_mutation_rate_is_the_stage1_component_not_the_total(path):
    """5.6e-9 here + 8.4e-9 from stage 2's unconditional cattle overlay = 1.4e-8,
    E/F/G's total. Category H's 1.4e-8 in this slot would land at 2.24e-8."""
    assert load(path)["mutation_rate"] == pytest.approx(5.6e-9)


@pytest.mark.parametrize("path", NEUTRAL_CATTLE, ids=os.path.basename)
def test_there_is_no_split_q_handoff(path):
    """Those keys describe E/F/G's piecewise time scale. This arm is one
    coalescent at one Q, and stage 2 does the right thing only when both are
    absent."""
    cfg = load(path)
    assert "handoff_ticks" not in cfg
    assert "deep_Q_scaling" not in cfg


@pytest.mark.parametrize("path", NEUTRAL_CATTLE, ids=os.path.basename)
def test_there_is_no_baseline_checkpoint_dependency(path):
    """I simulates all twelve epochs itself; E/F/G's ep7 handoff is not part of
    its DAG. A leftover key here would also make the ^(E|F|G)$ pre-flight guard
    look wrongly scoped."""
    cfg = load(path)
    assert "cattle_baseline_seed" not in cfg
    assert "cattle_baseline_search_dirs" not in cfg
    assert cfg["stage1_search_dirs"] == []


@pytest.mark.parametrize("path", NEUTRAL_CATTLE, ids=os.path.basename)
def test_the_sample_fits_inside_the_terminal_population(path):
    """9,000 = 8,000 GWAS + the largest GTEx panel, and also epoch 12 at
    Q_scaling 0.01. rules/stage1_cattle_neutral.smk refuses to exceed it."""
    cfg = load(path)
    largest_gtex = 1000 if int(cfg["gtex_size"]) == -1 else int(cfg["gtex_size"])
    n = int(cfg["gwas_size"]) + largest_gtex
    assert n <= cd.terminal_size(float(cfg["Q_scaling"]))


@pytest.mark.parametrize("path", NEUTRAL_CATTLE, ids=os.path.basename)
def test_it_uses_the_drawn_dfe_and_not_the_redistribution_model(path):
    """The two neutral models are mutually exclusive; the Snakefile rejects the
    pair, but only after stage 1 has already run."""
    cfg = load(path)
    assert cfg["synthetic_dfe_effects"] is True
    assert cfg["neutral_trait_vars"] is False


@pytest.mark.parametrize("path", NEUTRAL_CATTLE, ids=os.path.basename)
def test_it_is_routed_as_cattle_all_the_way_down(path):
    """`species` -- not `pipeline` -- picks the --cattle_ts_file branch, the c
    panel prefix and plain .vcf. Setting it to human would route an I run
    through the human branch and mislabel every downstream file."""
    cfg = load(path)
    assert cfg["species"] == "cattle"
    assert paths.species_prefix(cfg) == "c"
    assert paths.stage2_gwas_category(cfg) == "cgwas"
    assert paths.stage2_vcf(cfg, "cgwas").endswith("cgwas.vcf")


@pytest.mark.parametrize("path", NEUTRAL_CATTLE, ids=os.path.basename)
def test_the_seed_is_on_the_I_band(path):
    """I is the 9th letter, so I{N} -> 9{N}. seed_of() in the submit scripts
    derives the same thing from the letter; a config off the band would make the
    two disagree and land the run in a directory named for another replicate."""
    cfg = load(path)
    assert 91 <= int(cfg["stage1_seed"]) <= 99
    assert cfg["stage1_seed"] == cfg["stage2_seed"]


def test_the_cells_differ_only_in_the_multipliers_and_the_workdir():
    """The four configs are one arm at four phenotype-scaling cells. Anything
    else differing between them is a copy that drifted."""
    IGNORE = {"gwas_scaling", "gtex_scaling", "workdir", "publishdir"}
    base = load(os.path.join(CONFIG_DIR, "cattle_neutral_2Mb_r3.yaml"))
    for path in NEUTRAL_CATTLE:
        cfg = load(path)
        assert cfg.keys() == base.keys(), os.path.basename(path)
        differing = {k for k in base if cfg[k] != base[k]}
        assert differing <= IGNORE, (os.path.basename(path), differing - IGNORE)
