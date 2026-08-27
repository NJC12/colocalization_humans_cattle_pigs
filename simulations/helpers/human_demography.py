"""The human demographic models and the population sampled from them.

WHY THIS FILE EXISTS

`HomSap` / `OutOfAfrica_2T12` (Tennessen et al. 2012) was a pair of string
literals in two places -- human_simulation_o2.py (categories A/B/C/D, forward in
SLiM) and human_neutral_simulation.py (category H, a pure coalescent) -- which
between them are supposed to be the same demography with and without selection.
Category J added a third consumer and a reason for the arms to differ on purpose,
so the constants moved here where the divergence is one named parameter rather
than a diff between two files.

THE POPULATION

OutOfAfrica_2T12 has two populations, AFR and EUR, and simulates both forward
regardless of which one is sampled. Every arm through category I sampled EUR:

    samples = {"AFR": 0, "EUR": n_samples}

Category J samples the other side of the same run. Nothing else about the
simulation changes -- same model, same DFE, same recombination and mutation
rates, same sample size, same trait architecture, same stages 2-5 -- so A minus J
is the ancestry contrast and nothing else. Note that the AFR sample is not an
isolated African population: 2T12 carries migration between the two, so this is
the African branch of a two-population model rather than Africa_1T12.

THE MODEL

Categories M and N needed something J did not: a demography that is not in
stdpopsim's catalog at all. The Wang et al. 2014 NFE + Finnish founder model is
a `demes` graph shipped with this repo, so `population` alone stopped being
enough to say what stage 1 ran -- FIN is not a population OF OutOfAfrica_2T12.
MODELS below is the registry, and `demographic_model` is the config key that
picks an entry. It is pure data: a catalog entry names a stdpopsim model, a
demes entry names a YAML file relative to the simulations/ directory, and the
consumer (human_simulation_o2.py) decides how to turn either into a
stdpopsim.DemographicModel. Adding a parameterization is adding a dict entry.

Kept free of stdpopsim, msprime, demes and yaml (everything but `os`) for the
same reason helpers/cattle_demography.py is: the naming rule below has to be
testable without the simulation stack, since helpers/tests stubs those modules
out, and helpers/paths.py must stay importable anywhere Snakemake parses the
workflow -- including on a login node where nothing is installed.

THE FILENAME RULE

stage1_ts_name() is the single definition of what stage 1 writes. It matters
more than it looks. Stage-1 adoption via `stage1_search_dirs` matches on
filename plus embedded seed (helpers/search_dirs.py) -- and for human names the
seed half is a no-op, because search_dirs' seed pattern is `sd\\d+|seed_\\d+` and
`hts_11.ts` matches neither. So for this pipeline the FILENAME IS THE ONLY
GUARD. That is also why params_record.STAGE2_KEYS omits every demography key: it
delegates to this name.

    default model, EUR    ->  hts_{seed}.ts              A/B/C/D/K, unchanged
    default model, other  ->  hts_{pop}_{seed}.ts        J: hts_AFR_101.ts
    any other model       ->  hts_{tag}_{pop}_{seed}.ts  M: hts_finwang_FIN_131.ts

The tag is needed because two categories can share a model (M and N both run
Wang 2014); the population is needed because that is all M and N differ by. EUR
under the default model keeps the bare `hts_{seed}.ts` so every tree sequence
produced before either key existed stays adoptable.

No dots beyond the extension, ever: helpers/paths.stage2_run_tag() strips two
extensions off this name and then appends dot-joined segments, and
summarize_coloc.read_params() prefix-globs on the result.
"""

import os

#: Directory the `demes_file` paths in MODELS are relative to (simulations/).
_SIM_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: stdpopsim species. Every human arm, whatever its demography, is HomSap: the
#: species supplies the contig machinery, ploidy and generation time, and only
#: the demographic model is swapped.
SPECIES = "HomSap"

#: The human demographic models this pipeline can run.
#:
#: kind="catalog"  -> stdpopsim's own catalog, via species.get_demographic_model
#: kind="demes"    -> a demes YAML in this repo, via msprime.Demography.from_demes
#:
#: `tag` is the token that goes into the stage-1 filename. The default model's
#: tag is deliberately empty so its names are unchanged; every other model needs
#: a short, lowercase, dot-free one.
MODELS = {
    "OutOfAfrica_2T12": {
        "kind": "catalog",
        "catalog_id": "OutOfAfrica_2T12",
        "tag": "",
        "populations": ("AFR", "EUR"),
        "default_population": "EUR",
        "generation_time": 25,
        "description": "Tennessen et al. 2012 two-population out-of-Africa",
        "long_description": (
            "stdpopsim's HomSap/OutOfAfrica_2T12. Categories A/B/C/D/K sample "
            "EUR, category J samples AFR."
        ),
    },
    "FinnishWang2014": {
        "kind": "demes",
        "demes_file": "demography/finnish_demo_wang_2014.yaml",
        "tag": "finwang",
        "populations": ("NFE", "FIN"),
        "default_population": "FIN",
        "generation_time": 25,
        "description": "Wang et al. 2014 NFE + Finnish founder model",
        "long_description": (
            "Wang, Agarwala, Flannick, Chiang, Altshuler & Hirschhorn (2014), "
            "Am J Hum Genet 94:710-720, Table S2 class 3. A non-Finnish "
            "European population with a bottleneck to N=2,000 followed by 270 "
            "generations of 1.5%/generation growth, plus a Finnish founder "
            "population split off 100 generations ago at N=1,000 with one-way "
            "gene flow from NFE. Category M samples FIN, category N samples "
            "NFE; M minus N is the founder event with everything else held "
            "fixed. Three parameters are midpoints of ranges the paper "
            "reports -- see the notes in the demes file."
        ),
    },
}

#: The model every arm through category L ran, and what a config without a
#: `demographic_model` key still means.
DEFAULT_MODEL = "OutOfAfrica_2T12"

#: (model, population, Q_scaling) -> individuals alive at the sampling tick.
#:
#: MEASURED BY ASKING SLiM, not computed. stdpopsim's SLiM engine refuses to
#: sample more individuals than the rescaled deme holds at the sampling tick,
#: and `present_size / Q` is NOT that number: stdpopsim rounds epoch boundaries
#: to whole ticks, which truncates the tail of an exponential growth phase. For
#: the Wang model the naive figure overshoots by 13-34%:
#:
#:     Q     FIN alive   naive 34,458/Q      NFE alive   naive 111,394/Q
#:     10        2,266            3,446         10,653            11,139
#:      5        4,531            6,892         21,306            22,279
#:      4        7,491            8,614         27,437            27,848
#:      3        9,988           11,486         36,583            37,131
#:      1       34,458           34,458        111,394           111,394
#:
#: which is why M and N run Q=3 against the 9,000 stage 1 samples -- a margin of
#: 988, where the naive arithmetic would have promised 2,486. Q=4 looked like it
#: might fit and does not.
#:
#: TO ADD AN ENTRY, ask SLiM rather than dividing. Request an absurd sample and
#: read the number out of the error; a dry run stops right after setup, which is
#: where the check lives, so it costs seconds:
#:
#:     engine.simulate(model, contig, sample_dict(pop, 10**7, model_id),
#:                     slim_scaling_factor=Q, dry_run=True, seed=1)
#:     # ERROR: Request to sample 10000000 individuals from p1 (FIN) at tick
#:     #        27089, but only 9988 individuals will be alive.
SLIM_CAPACITY = {
    ("FinnishWang2014", "FIN", 3): 9988,
    ("FinnishWang2014", "NFE", 3): 36583,
}

#: Back-compatible aliases for the default model. human_neutral_simulation.py
#: keeps its own copies on purpose (category H must not follow A if A's
#: demography is ever changed); these are for callers that predate MODELS.
DEMOGRAPHIC_MODEL = DEFAULT_MODEL
POPULATIONS = MODELS[DEFAULT_MODEL]["populations"]
DEFAULT_POPULATION = MODELS[DEFAULT_MODEL]["default_population"]


def validate_model(model):
    """Return ``model``, or raise ValueError naming the registered models."""
    if model not in MODELS:
        raise ValueError(
            f"demographic_model must be one of {list(MODELS)}, got {model!r}"
        )
    return model


def spec(model=DEFAULT_MODEL):
    """The registry entry for ``model``."""
    return MODELS[validate_model(model)]


def populations(model=DEFAULT_MODEL):
    """Every population ``model`` defines, in the order it lists them."""
    return spec(model)["populations"]


def default_population(model=DEFAULT_MODEL):
    """The population a config gets if it names ``model`` but no population."""
    return spec(model)["default_population"]


def model_tag(model=DEFAULT_MODEL):
    """The filename token for ``model``; ``""`` for the default one."""
    return spec(model)["tag"]


def demes_path(model):
    """Absolute path to ``model``'s demes YAML.

    Resolved against simulations/ rather than stored absolute, so no config
    carries a machine-specific path and the same registry works locally and on
    O2. Raises for catalog models, which have no file.
    """
    s = spec(model)
    if s["kind"] != "demes":
        raise ValueError(
            f"demographic_model {model!r} is a stdpopsim catalog model "
            f"({s['catalog_id']}); it has no demes file"
        )
    return os.path.join(_SIM_DIR, s["demes_file"])


def validate_population(population, model=DEFAULT_MODEL):
    """Return ``population``, or raise ValueError naming the valid choices.

    A typo would otherwise surface as a KeyError inside stdpopsim's sampling
    dict, several minutes into a submitted job.
    """
    valid = populations(model)
    if population not in valid:
        raise ValueError(
            f"population must be one of {list(valid)} (the populations "
            f"{model} defines), got {population!r}"
        )
    return population


def sample_dict(population, n_samples, model=DEFAULT_MODEL):
    """stdpopsim's ``samples`` argument: ``n_samples`` from one population, 0 elsewhere.

    Every population of the model appears whichever is sampled, matching the
    shape human_simulation_o2.py used before the population was a parameter
    (``{"AFR": 0, "EUR": n_samples}``).

    NOTE THE CAPACITY LIMIT this feeds into: stdpopsim's SLiM engine refuses to
    sample more individuals than the Q-RESCALED deme holds at the sampling tick
    ("Request to sample N individuals from pX ... but only M will be alive").
    That is what pins categories M and N to Q_scaling 3 -- FIN holds 9,988
    individuals there against the 9,000 requested, and 7,491 at Q=4. See
    SLIM_CAPACITY above; helpers/tests/test_finnish_demography.py asserts it for
    every config that names a demes model.
    """
    validate_population(population, model)
    return {p: (n_samples if p == population else 0)
            for p in populations(model)}


def stage1_ts_name(population, seed, model=DEFAULT_MODEL):
    """The stage-1 tree sequence filename. See THE FILENAME RULE above.

    default model + EUR  -> ``hts_{seed}.ts``            (keeps A/B/C/D adoptable)
    default model + AFR  -> ``hts_AFR_{seed}.ts``        (category J)
    any other model      -> ``hts_{tag}_{pop}_{seed}.ts``  ``hts_finwang_FIN_131.ts``
    """
    validate_population(population, model)
    tag = model_tag(model)
    if tag:
        return f"hts_{tag}_{population}_{seed}.ts"
    if population == default_population(model):
        return f"hts_{seed}.ts"
    return f"hts_{population}_{seed}.ts"
