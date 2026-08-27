"""Categories M and N: the Wang 2014 Finnish founder pair.

WHAT CAN GO WRONG SILENTLY HERE

1. THE SAMPLING CAP. stdpopsim's SLiM engine refuses to sample more individuals
   than the Q-RESCALED deme holds at the sampling tick, and it finds out at the
   very end of a multi-hour forward simulation. FIN's present-day 34,458 is
   2,266 individuals at Q=10 -- and 7,491 at Q=4 -- against the 9,000 stage 1
   asks for; that is why these configs run Q_scaling 3, where it holds 9,988.
   `test_the_sampled_deme_can_supply_the_sample` is the cheap version of that
   check, and it runs against every config. Note it asserts against MEASURED
   numbers: present_size/Q is 13-34% optimistic, because stdpopsim rounds epoch
   boundaries to whole ticks and truncates the growth phase's tail.

2. THE FILENAME. Human stage-1 adoption matches on filename plus embedded seed,
   and search_dirs' seed pattern (`sd\\d+|seed_\\d+`) matches no `hts_*` name at
   all -- so for this pipeline the filename is the ONLY guard. M and N share a
   demographic model and differ only by deme, so the name must carry both. A
   name that carried only the model would let M adopt N's tree sequence.

3. THE M/N IDENTITY. M minus N is the founder event and nothing else. If the
   two configs drift apart on anything but the deme, the name and the seeds,
   that contrast becomes an uninterpretable mixture -- the same failure the
   A/J test guards against.
"""

import glob
import os

import pytest
import yaml

from helpers import human_demography as hd
from helpers import paths

SIM_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CONFIG_DIR = os.path.join(SIM_DIR, "config")

ALL_CONFIGS = sorted(glob.glob(os.path.join(CONFIG_DIR, "*.yaml")))
FIN_CONFIGS = sorted(glob.glob(os.path.join(CONFIG_DIR, "human_fin_*.yaml")))
NFE_CONFIGS = sorted(glob.glob(os.path.join(CONFIG_DIR, "human_nfe_*.yaml")))

MODEL = "FinnishWang2014"


def load(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


# ------------------------------------------------------------ the registry

def test_the_default_model_is_unchanged():
    """Every arm through category L reads these three aliases. Moving them
    would rename tree sequences that already exist."""
    assert hd.DEFAULT_MODEL == "OutOfAfrica_2T12"
    assert hd.DEMOGRAPHIC_MODEL == "OutOfAfrica_2T12"
    assert hd.POPULATIONS == ("AFR", "EUR")
    assert hd.DEFAULT_POPULATION == "EUR"
    assert hd.SPECIES == "HomSap"


@pytest.mark.parametrize("model", sorted(hd.MODELS))
def test_every_registry_entry_is_well_formed(model):
    spec = hd.MODELS[model]
    assert spec["kind"] in {"catalog", "demes"}
    assert spec["populations"], "a model with no populations cannot be sampled"
    assert spec["default_population"] in spec["populations"]
    assert isinstance(spec["generation_time"], (int, float))
    for key in ("description", "long_description"):
        assert spec[key], f"{model} needs a {key} for stdpopsim.DemographicModel"
    if spec["kind"] == "catalog":
        assert spec["catalog_id"]
    else:
        assert os.path.exists(hd.demes_path(model)), hd.demes_path(model)


def test_model_tags_are_unique_and_path_safe():
    """The tag goes into a filename that stage2_run_tag splits on dots and that
    summarize_coloc prefix-globs. Only the default model may have an empty tag;
    two models sharing one would collide."""
    tags = [hd.MODELS[m]["tag"] for m in hd.MODELS]
    assert len(tags) == len(set(tags))
    assert hd.MODELS[hd.DEFAULT_MODEL]["tag"] == ""
    for model, tag in zip(hd.MODELS, tags):
        if model == hd.DEFAULT_MODEL:
            continue
        assert tag and tag == tag.lower()
        assert not (set(tag) & set("._ /")), f"{tag!r} is not path-safe"


def test_a_catalog_model_has_no_demes_file():
    with pytest.raises(ValueError, match="catalog"):
        hd.demes_path(hd.DEFAULT_MODEL)


def test_an_unknown_model_is_refused_by_name():
    with pytest.raises(ValueError, match="demographic_model must be one of"):
        hd.validate_model("Finnish")


# ------------------------------------------------------------- the naming rule

def test_the_pre_existing_names_are_unchanged():
    """Back-compat pins. A/B/C/D/K and J already have tree sequences on disk."""
    assert hd.stage1_ts_name("EUR", 11) == "hts_11.ts"
    assert hd.stage1_ts_name("AFR", 101) == "hts_AFR_101.ts"
    assert paths.stage1_human_ts({"stage1_seed": 11}) == "hts_11.ts"


def test_the_new_names_carry_the_model_and_the_deme():
    assert hd.stage1_ts_name("FIN", 131, MODEL) == "hts_finwang_FIN_131.ts"
    assert hd.stage1_ts_name("NFE", 141, MODEL) == "hts_finwang_NFE_141.ts"


def test_m_and_n_cannot_adopt_each_other_or_anyone_else():
    """The adoption invariant. Every stage-1 name at the same seed must differ,
    because the seed half of the match never fires for human filenames."""
    seed = 131
    names = {
        hd.stage1_ts_name("EUR", seed),
        hd.stage1_ts_name("AFR", seed),
        hd.stage1_ts_name("FIN", seed, MODEL),
        hd.stage1_ts_name("NFE", seed, MODEL),
        f"nts_{seed}.ts",     # category H
        f"cnts_{seed}.ts",    # category I
    }
    assert len(names) == 6


@pytest.mark.parametrize("population,seed", [("FIN", 131), ("NFE", 141)])
def test_the_stage1_filename_has_no_extra_dots(population, seed):
    """stage2_run_tag strips two extensions and then appends dot-joined
    segments; summarize_coloc.read_params prefix-globs the result."""
    name = hd.stage1_ts_name(population, seed, MODEL)
    assert name.count(".") == 1
    tag = paths.stage2_run_tag({"pipeline": "human", "demographic_model": MODEL,
                                "population": population, "stage1_seed": seed,
                                "causal_min_maf": 0.01})
    assert tag == f"hts_finwang_{population}_{seed}"


def test_a_population_of_the_wrong_model_is_refused():
    """EUR is not a population of the Wang model, and FIN is not one of 2T12.
    Without this the typo becomes a KeyError inside stdpopsim's sampling dict,
    minutes into a submitted job."""
    with pytest.raises(ValueError, match="population must be one of"):
        hd.stage1_ts_name("EUR", 131, MODEL)
    with pytest.raises(ValueError, match="population must be one of"):
        hd.stage1_ts_name("FIN", 131)


def test_sample_dict_covers_every_deme_of_the_chosen_model():
    assert hd.sample_dict("FIN", 9000, MODEL) == {"NFE": 0, "FIN": 9000}
    assert hd.sample_dict("NFE", 9000, MODEL) == {"NFE": 9000, "FIN": 0}
    assert hd.sample_dict("EUR", 9000) == {"AFR": 0, "EUR": 9000}


# --------------------------------------------------------------- the demes file

def demes_graph():
    with open(hd.demes_path(MODEL)) as fh:
        return yaml.safe_load(fh)


def test_the_demes_file_is_the_model_the_configs_claim():
    g = demes_graph()
    assert g["time_units"] == "generations", \
        "the registry's generation_time assumes the graph is already in generations"
    assert [d["name"] for d in g["demes"]] == list(hd.MODELS[MODEL]["populations"])


def test_the_present_day_sizes_are_the_ones_q_was_chosen_against():
    """Pinning these is the point. Q_scaling 3 is not a preference -- it is the
    smallest integer at which FIN can supply 9,000 individuals. Editing the
    growth rates in the demes file moves that number, and nothing else in the
    pipeline would notice."""
    sizes = {d["name"]: _present_size(d) for d in demes_graph()["demes"]}
    assert sizes == {"NFE": 111394, "FIN": 34458}


def _present_size(deme):
    """The deme's size at time 0, which is what the sampling cap is measured on."""
    last = deme["epochs"][-1]
    assert last["end_time"] == 0, f"{deme['name']} does not reach the present"
    return last.get("end_size", last.get("start_size"))


# ------------------------------------------------- the cap, over every config

def _n_samples(cfg):
    """Mirror of rules/stage1_human.smk:_human_stage1_n_samples."""
    gtex = int(cfg.get("gtex_size", -1))
    return int(cfg["gwas_size"]) + (1000 if gtex == -1 else gtex)


@pytest.mark.parametrize("path", ALL_CONFIGS, ids=os.path.basename)
def test_the_sampled_deme_can_supply_the_sample(path):
    """THE TEST THAT WOULD HAVE CAUGHT THE BLOCKER.

    stdpopsim's SLiM engine errors with "Request to sample N individuals from
    pX ... but only M individuals will be alive". Asserted against
    hd.SLIM_CAPACITY, which holds numbers MEASURED from SLiM rather than
    computed -- `present_size / Q` overshoots by 13-34% because stdpopsim rounds
    epoch boundaries to whole ticks and truncates the tail of an exponential
    growth phase. See the table on SLIM_CAPACITY.

    A config whose (model, population, Q) has no measured entry FAILS rather
    than falling back to the arithmetic. Guessing is what the arithmetic is, and
    the whole point of this test is that the guess was 15% optimistic.
    """
    cfg = load(path)
    if cfg.get("pipeline") != "human":
        return
    model = cfg.get("demographic_model") or hd.DEFAULT_MODEL
    if hd.MODELS[model]["kind"] != "demes":
        return
    population = cfg.get("population") or hd.default_population(model)
    key = (model, population, int(cfg["Q_scaling"]))
    assert key in hd.SLIM_CAPACITY, (
        f"{os.path.basename(path)}: no measured capacity for {key}. Measure it "
        "with a SLiM dry run (recipe on helpers.human_demography.SLIM_CAPACITY) "
        "rather than dividing the present-day size by Q -- that overestimates."
    )
    available, wanted = hd.SLIM_CAPACITY[key], _n_samples(cfg)
    assert wanted <= available, (
        f"{os.path.basename(path)}: asks for {wanted:,} {population} individuals "
        f"but Q_scaling={cfg['Q_scaling']} leaves only {available:,} alive at "
        "the sampling tick"
    )


@pytest.mark.parametrize("key", sorted(hd.SLIM_CAPACITY))
def test_the_measured_capacity_is_below_the_naive_one(key):
    """Documents the DIRECTION of the error, and catches a transcription slip.

    Tick rounding can only lose individuals relative to present_size/Q, never
    add them. A measured number above the naive one means someone typed it
    wrong -- and typed it wrong in the dangerous direction."""
    model, population, Q = key
    with open(hd.demes_path(model)) as fh:
        graph = yaml.safe_load(fh)
    present = {d["name"]: _present_size(d) for d in graph["demes"]}[population]
    naive = round(present / Q)
    assert 0 < hd.SLIM_CAPACITY[key] <= naive, (
        f"{key}: measured {hd.SLIM_CAPACITY[key]:,} is above the naive "
        f"{naive:,}, which tick rounding cannot produce"
    )


# ----------------------------------------------------------------- the configs

def test_both_arms_of_the_pair_exist():
    """A missing cell is a config_of() miss in four submit scripts, which fails
    loudly -- but only once someone runs that cell."""
    assert FIN_CONFIGS, "expected at least config/human_fin_2Mb_g5t20_r3.yaml"
    assert NFE_CONFIGS, "expected at least config/human_nfe_2Mb_g5t20_r3.yaml"
    assert len(FIN_CONFIGS) == len(NFE_CONFIGS), "the pair must stay in step"


@pytest.mark.parametrize("path", FIN_CONFIGS, ids=os.path.basename)
def test_m_is_n_with_one_parameter_changed(path):
    """The identity IS the experiment. If these drift apart, M - N stops being
    the founder-event contrast."""
    m = load(path)
    n = load(path.replace("human_fin_", "human_nfe_"))

    assert m["population"] == "FIN"
    assert n["population"] == "NFE"
    assert m["demographic_model"] == n["demographic_model"] == MODEL
    assert m["pipeline"] == n["pipeline"] == "human"
    assert m["species"] == n["species"] == "human"

    differ = {"basename", "population", "stage1_seed", "stage2_seed",
              "workdir", "publishdir"}
    for k in sorted(set(m) | set(n)):
        if k in differ:
            continue
        assert m.get(k) == n.get(k), \
            f"{os.path.basename(path)} diverges from its NFE twin on {k!r}"


@pytest.mark.parametrize("path", FIN_CONFIGS + NFE_CONFIGS, ids=os.path.basename)
def test_the_pair_is_category_a_plus_a_demography_and_a_scaling(path):
    """Everything except the demography, the scaling forced by it, the name and
    the seeds must match A. An accidental extra divergence fails here rather
    than turning up as an unexplained difference in the results."""
    arm = load(path)
    a = load(path.replace("human_fin_", "human_").replace("human_nfe_", "human_"))

    differ = {"basename", "demographic_model", "population", "Q_scaling",
              "stage1_seed", "stage2_seed", "workdir", "publishdir"}
    for k in sorted(set(a) | set(arm)):
        if k in differ:
            continue
        assert a.get(k) == arm.get(k), \
            f"{os.path.basename(path)} diverges from A on {k!r}"


@pytest.mark.parametrize("path", FIN_CONFIGS + NFE_CONFIGS, ids=os.path.basename)
def test_the_scaling_is_the_forced_one(path):
    """Q_scaling 3 is load-bearing, not incidental. It is also in
    params_record.STAGE2_KEYS, so raising it silently would be refused by the
    provenance guard rather than corrupting a stage 2 -- but only after the
    stage-1 job has already failed on the sampling cap."""
    assert load(path)["Q_scaling"] == 3


@pytest.mark.parametrize("path", FIN_CONFIGS + NFE_CONFIGS, ids=os.path.basename)
def test_the_pair_builds_its_own_stage1(path):
    """There is no earlier Finnish tree sequence anywhere to adopt."""
    assert load(path).get("stage1_search_dirs") in ([], None)


@pytest.mark.parametrize("path,lo,hi", [(p, 131, 139) for p in FIN_CONFIGS]
                                       + [(p, 141, 149) for p in NFE_CONFIGS],
                         ids=lambda v: os.path.basename(v) if isinstance(v, str) else str(v))
def test_the_seeds_are_on_their_own_bands(path, lo, hi):
    """seed = 10*letter_index + replicate, M the 13th letter and N the 14th."""
    seed = int(load(path)["stage1_seed"])
    assert lo <= seed <= hi
    assert seed > 129, "must clear L's 121-129 band and everything above it"


# ------------------------------------------------------ downstream is wired up

@pytest.mark.parametrize("letter,seed", [("M", 131), ("N", 141)])
def test_the_seed_bands_follow_the_repo_rule(letter, seed):
    idx = ord(letter) - ord('A') + 1
    assert 10 * idx + 1 == seed


@pytest.mark.parametrize("letter", ["M", "N"])
def test_the_categories_are_known_to_the_summariser(letter):
    """A missing entry degrades to demography='?' rather than erroring."""
    from helpers import summarize_coloc
    assert summarize_coloc.DEMOGRAPHY[letter] == "human"


@pytest.mark.parametrize("script", [
    "submit_2Mb_r3_cmaf_replicates.sh",
    "submit_2Mb_r3_cmaf_fm001.sh",
    "submit_2Mb_r3_x20_psamp_fm001.sh",
    "submit_2Mb_r3_cmaf01_control.sh",
])
def test_every_submit_script_can_resolve_the_pair(script):
    sh = open(os.path.join(SIM_DIR, script)).read()
    assert 'M) echo "13${n}" ;; N) echo "14${n}"' in sh
    assert "config/human_fin_2Mb_${CELL}_r3.yaml" in sh
    assert "config/human_nfe_2Mb_${CELL}_r3.yaml" in sh


@pytest.mark.parametrize("script,marker", [
    ("submit_2Mb_r3_cmaf_fm001.sh", 'STAGE1_SRC["$ID"]=""'),
    ("submit_2Mb_r3_cmaf_fm001.sh", 'STAGE2_SRC["$ID"]=""'),
    ("submit_2Mb_r3_x20_psamp_fm001.sh", 'STAGE1_SRC["$ID"]=""'),
    ("submit_2Mb_r3_cmaf01_control.sh", 'STAGE1_SRC["$ID"]=""'),
])
def test_the_pair_is_never_handed_someone_elses_upstream(script, marker,
                                                        guard_letters):
    """M and N have new seeds AND a new demography, so nothing existing can
    satisfy them -- and a search dir that finds nothing fails silently."""
    sh = open(os.path.join(SIM_DIR, script)).read()
    assert {"M", "N"} <= guard_letters(sh, marker)


def test_the_script_builds_the_model_from_the_registry():
    """human_simulation_o2.py must not hardcode either the catalog id or the
    demes path -- both live in helpers/human_demography.MODELS."""
    src = open(os.path.join(SIM_DIR, "human_simulation_o2.py")).read()
    assert "human_demography.demes_path(" in src
    assert "msprime.Demography.from_demes(" in src
    assert "human_demography.stage1_ts_name(" in src
    assert 'species.get_demographic_model(model_spec["catalog_id"])' in src, \
        "the catalog id comes from the registry entry, not from a literal"
