"""`require_gtex_partner`: the pairing of causal loci across the two panels.

Whether a central causal locus must segregate in the GTEx panel used to be
DERIVED rather than set:

    topup_gtex = power_sampling or synthetic_effects

which made two experimental conditions inexpressible -- a uniform draw WITHOUT
the intersection, and a synthetic-DFE or power draw WITH it. Both are arms of the
publication design, which varies the pairing on purpose, so the rule is now a
knob whose `auto` value reproduces the old behavior exactly.

Two properties carry the whole change and both are pinned here.

1. `auto` IS the old rule, in all four scheme combinations, and emits no path
   segment. Every stage-2..5 directory ever written is at `auto`, so if this
   slips the pipeline stops finding its own outputs.
2. A value that DIFFERS from the derived one gets its own path segment. Arms 1
   and 2 of the publication run differ in nothing else, so without the segment
   they would write the same filenames into a flat archive layout and land on
   top of each other.
"""

import itertools
import os
import re

import pytest

from helpers import params_record, paths


SIM_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
STAGE2_SCRIPT = os.path.join(SIM_DIR, "create_gwas_files_and_phenotypes.py")
SNAKEFILE = os.path.join(SIM_DIR, "Snakefile")
COMMON_SMK = os.path.join(SIM_DIR, "rules", "common.smk")

BASE = dict(pipeline="human", stage1_seed=11, causal_min_maf=0.001,
            n_central_traits=25)

SCHEMES = {
    "uniform_selected":  dict(),
    "uniform_synthetic": dict(synthetic_dfe_effects=True),
    "power_selected":    dict(causal_sampling="power", sampling_gwas_n=200000),
    "power_synthetic":   dict(causal_sampling="power", sampling_gwas_n=200000,
                              synthetic_dfe_effects=True),
}

# The derived rule, written out independently of the implementation so this is a
# real check and not a tautology: intersect iff uniform AND non-synthetic.
DERIVED = {
    "uniform_selected":  True,
    "uniform_synthetic": False,
    "power_selected":    False,
    "power_synthetic":   False,
}


# --------------------------------------------------------------- `auto` is inert

@pytest.mark.parametrize("scheme", sorted(SCHEMES))
def test_auto_reproduces_the_derived_rule(scheme):
    cfg = {**BASE, **SCHEMES[scheme]}
    assert paths.require_gtex_partner_auto(cfg) is DERIVED[scheme]


@pytest.mark.parametrize("scheme", sorted(SCHEMES))
def test_auto_emits_no_segment(scheme):
    """Not one of the ~55 existing arms may move."""
    cfg = {**BASE, **SCHEMES[scheme]}
    assert paths.require_gtex_partner_segment(cfg) == ""
    assert paths.require_gtex_partner_segment(
        {**cfg, "require_gtex_partner": "auto"}) == ""


@pytest.mark.parametrize("scheme", sorted(SCHEMES))
def test_a_missing_key_is_auto(scheme):
    """paths.py is handed raw config YAML in places that never ran the
    Snakefile's defaults, and a --config invocation may omit the key."""
    cfg = {**BASE, **SCHEMES[scheme]}
    assert "require_gtex_partner" not in cfg
    assert paths.require_gtex_partner(cfg) == paths.LEGACY_REQUIRE_GTEX_PARTNER
    assert paths.require_gtex_partner_segment(cfg) == ""


@pytest.mark.parametrize("scheme", sorted(SCHEMES))
def test_stating_the_derived_value_explicitly_is_also_inert(scheme):
    """The launcher passes this key for all 120 runs, including where it agrees
    with the derived rule. If agreement still emitted a segment, half those runs
    would land in directories no existing tool knows how to read."""
    cfg = {**BASE, **SCHEMES[scheme]}
    assert paths.require_gtex_partner_segment(
        {**cfg, "require_gtex_partner": DERIVED[scheme]}) == ""


# ------------------------------------------------- a real override is visible

@pytest.mark.parametrize("scheme", sorted(SCHEMES))
def test_overriding_the_derived_rule_emits_a_segment(scheme):
    cfg = {**BASE, **SCHEMES[scheme]}
    flipped = not DERIVED[scheme]
    seg = paths.require_gtex_partner_segment({**cfg, "require_gtex_partner": flipped})
    assert seg == (".gtexreq_1" if flipped else ".gtexreq_0")


def test_the_publication_arms_1_and_2_cannot_collide():
    """They differ in NOTHING else -- same floors, same cell, same trait counts,
    same seed, same draw. The segment is the only thing separating them."""
    for scheme in SCHEMES:
        cfg = {**BASE, **SCHEMES[scheme]}
        paired = paths.stage2_run_tag({**cfg, "require_gtex_partner": True})
        unpaired = paths.stage2_run_tag({**cfg, "require_gtex_partner": False})
        assert paired != unpaired, scheme


def test_the_segment_comes_last_so_old_tags_stay_prefixes():
    """summarize_coloc.read_params prefix-globs `{tag}.*.params.txt`, so a
    segment inserted before an existing one would orphan the record it is
    supposed to find."""
    cfg = {**BASE, "synthetic_dfe_effects": True, "neutral_keep_fraction": 0.306,
           "require_gtex_partner": True}
    tag = paths.stage2_run_tag(cfg)
    assert tag.endswith(".gtexreq_1")
    without = paths.stage2_run_tag({k: v for k, v in cfg.items()
                                    if k != "require_gtex_partner"})
    assert tag.startswith(without)


def test_the_full_publication_grid_is_distinct_within_an_arm():
    """A and K share a seed by design (K adopts A's tree sequence), so their run
    tags can and do coincide across arms -- that is safe, they live in different
    workdirs. What must NOT happen is two cells of the SAME arm colliding."""
    arms = {
        "causal_maf001_paired":   dict(causal_min_maf=0.001, require_gtex_partner=True),
        "causal_maf001_unpaired": dict(causal_min_maf=0.001, require_gtex_partner=False),
        "causal_power_n200000":   dict(causal_min_maf=0, causal_sampling="power",
                                       sampling_gwas_n=200000, require_gtex_partner=False),
        "causal_power_n30000":    dict(causal_min_maf=0, causal_sampling="power",
                                       sampling_gwas_n=30000, require_gtex_partner=False),
    }
    for arm, over in arms.items():
        a = paths.stage2_run_tag({**BASE, **over})
        k = paths.stage2_run_tag({**BASE, **over, "synthetic_dfe_effects": True})
        # A and K are different CATEGORIES in one arm; if their tags matched they
        # would still be in different workdirs, but stage-2 adoption keys on the
        # tag, so the pre-flight has to keep stage2_search_dirs empty. Assert the
        # tags differ wherever the pairing rule distinguishes them, and record
        # where it does not.
        if arm.startswith("causal_maf001"):
            assert a != k, arm
        else:
            # Arms 3 and 4: the pairing is False for both, so nothing separates
            # them. This is the documented reason previous_workdirs must not be
            # used between A and K.
            assert a == k, arm


# ------------------------------------------------------ wiring, not just paths

def test_the_key_is_a_stage2_key():
    """It decides which loci become causal, so a stage-2 directory built under a
    different value must be refused rather than silently reused."""
    assert "require_gtex_partner" in params_record.STAGE2_KEYS
    assert "require_gtex_partner" not in params_record.ANALYSIS_KEYS


def test_two_arms_get_different_stage2_uids():
    """This is what the Snakefile's parse-time guard compares."""
    cfg = {**BASE, "gwas_scaling": 5, "gtex_scaling": 20, "stage2_seed": 11}
    paired = params_record.uid({**cfg, "require_gtex_partner": True},
                               params_record.STAGE2_KEYS)
    unpaired = params_record.uid({**cfg, "require_gtex_partner": False},
                                 params_record.STAGE2_KEYS)
    assert paired != unpaired


def test_stage2_records_the_unresolved_value():
    """The values dict is config-keyed so its uid matches the one params_record
    computes from the config, and the config holds 'auto'. Recording the resolved
    bool there would make every 'auto' run's uid disagree with its own record."""
    src = open(STAGE2_SCRIPT).read()
    m = re.search(r"values=\{(.*?)\n        \}", src, re.S)
    assert m, "could not find the stage-2 values dict"
    assert "'require_gtex_partner': args.require_gtex_partner," in m.group(1)


def test_the_flag_is_suppressed_at_auto():
    """Emitting nothing at the legacy value is what keeps every pre-existing
    run's shell command byte-identical -- the same discipline as synthetic_flag,
    sampling_flag and thin_flag."""
    src = open(COMMON_SMK).read()
    assert 'partner_flag' in src
    assert '"" if config["require_gtex_partner"] == "auto" else' in src
    assert "{params.partner_flag}" in src


def test_the_snakefile_normalises_rather_than_trusting_the_parse():
    """The string "False" is truthy. If --config ever hands the unparsed string
    through, `if config[...]` would read arm 2 as arm 1 -- and since the two
    differ in nothing else, nothing downstream would look wrong."""
    src = open(SNAKEFILE).read()
    assert 'config.setdefault("require_gtex_partner"' in src
    assert '"false": False' in src


def test_the_partner_table_is_written_whenever_the_rule_is_stated():
    """Three of the four publication arms are unpaired, so no scorer may have to
    infer pairing from the file's ABSENCE -- which is ambiguous between
    "intersected, so all-True" and "written by a version that did not emit it"."""
    src = open(STAGE2_SCRIPT).read()
    assert "if topup_gtex or args.require_gtex_partner != 'auto':" in src


def test_the_derived_rule_survives_only_as_a_comment():
    """If the old assignment is still live somewhere, the knob does nothing."""
    src = open(STAGE2_SCRIPT).read()
    for line in src.splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        assert stripped != "topup_gtex = power_sampling or synthetic_effects", (
            "the derived rule is still executing; the knob is dead code")
