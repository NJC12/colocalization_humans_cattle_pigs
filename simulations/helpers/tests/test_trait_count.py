"""The central-trait count has to reach the path, because it reaches no filename.

`causal_min_maf` is in every stage-2 filename (`..._maf_0.001.tsv`) and in the
stage-2 directory, so two runs that differ in it are visibly different on disk.
`n_central_traits` is in NEITHER. Stage-2 outputs are named
`gwas_G_gtex_T_maf_M/` and `{panel}_vars_gwas_G_gtex_T_maf_M.tsv`, so a 25-locus
run and a 50-locus run at otherwise identical settings write byte-for-byte the
same names -- and the archive layout is FLAT (`<arm>/<replicate>/<file>`), so
they would land on top of each other with nothing to tell them apart.

That is what `trait_count_segment` exists to prevent, and what these tests pin.
"""

import glob
import os

import pytest
import yaml

from helpers import params_record, paths


SIM_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CONFIG_DIR = os.path.join(SIM_DIR, "config")
ALL_CONFIGS = sorted(glob.glob(os.path.join(CONFIG_DIR, "*.yaml")))

BASE = dict(pipeline="human", stage1_seed=11, causal_min_maf=0,
            causal_sampling="power", sampling_gwas_n=100000)


def test_the_historical_count_emits_no_segment():
    """Every config through round 3 uses 50, so 50 must stay invisible or every
    path already on disk moves."""
    assert paths.LEGACY_N_CENTRAL_TRAITS == 50
    assert paths.trait_count_segment({**BASE, "n_central_traits": 50}) == ""


def test_a_missing_key_is_treated_as_historical():
    """paths.py is fed raw config YAML in places that never ran the Snakefile's
    defaults, and a --config invocation may omit the key entirely."""
    assert paths.trait_count_segment(BASE) == ""
    assert paths.trait_count_segment({**BASE, "n_central_traits": None}) == ""


def test_a_halved_arm_gets_its_own_segment():
    assert paths.trait_count_segment({**BASE, "n_central_traits": 25}) == ".ntr_25"


def test_the_two_arms_cannot_share_a_stage2_directory():
    """The whole point: same cell, same floors, same sampling, different count."""
    fifty = paths.stage2_run_tag({**BASE, "n_central_traits": 50})
    twenty_five = paths.stage2_run_tag({**BASE, "n_central_traits": 25})
    assert fifty == "hts_11.cmaf_0.psamp_100000"
    assert twenty_five == "hts_11.cmaf_0.psamp_100000.ntr_25"
    assert fifty != twenty_five


def test_the_segment_composes_with_the_others_in_a_fixed_order():
    """causal MAF, then sampling scheme, then trait count. A reordering would
    silently orphan every directory written under the old order."""
    tag = paths.stage2_run_tag({**BASE, "causal_min_maf": 0.001,
                                "sampling_power_plateau": 0.5,
                                "n_central_traits": 25})
    assert tag == "hts_11.cmaf_0.001.psamp_100000_sat05.ntr_25"
    assert tag.index("cmaf_") < tag.index("psamp_") < tag.index("ntr_")


def test_the_count_is_a_stage2_determining_key():
    """Two runs differing only in it must not adopt each other's phenotypes."""
    assert "n_central_traits" in params_record.STAGE2_KEYS
    assert "n_flank_gtex_traits" in params_record.STAGE2_KEYS
    base = {k: 1 for k in params_record.STAGE2_KEYS}
    a = dict(base, n_central_traits=50)
    b = dict(base, n_central_traits=25)
    assert params_record.uid(a, params_record.STAGE2_KEYS) != \
           params_record.uid(b, params_record.STAGE2_KEYS)


def test_the_flank_count_is_guarded_even_though_it_is_not_in_the_path():
    """trait_count_segment carries the headline knob only, matching
    causal_sampling_segment's rule. The flank count is therefore protected by the
    provenance guard instead of by the path, so it must be in STAGE2_KEYS."""
    seg_a = paths.trait_count_segment({**BASE, "n_central_traits": 25,
                                       "n_flank_gtex_traits": 25})
    seg_b = paths.trait_count_segment({**BASE, "n_central_traits": 25,
                                       "n_flank_gtex_traits": 50})
    assert seg_a == seg_b == ".ntr_25"
    base = {k: 1 for k in params_record.STAGE2_KEYS}
    assert params_record.uid(dict(base, n_flank_gtex_traits=25), params_record.STAGE2_KEYS) != \
           params_record.uid(dict(base, n_flank_gtex_traits=50), params_record.STAGE2_KEYS)


@pytest.mark.parametrize("path", ALL_CONFIGS, ids=os.path.basename)
def test_every_existing_config_still_emits_no_trait_segment(path):
    """The regression that matters: if any shipped config stopped hashing to the
    empty segment, its stage-2..5 tree would move and every result under it
    would be orphaned."""
    assert paths.trait_count_segment(yaml.safe_load(open(path))) == ""


def test_the_submit_script_and_paths_agree_on_the_suffix():
    """The root name in submit_2Mb_r3_x20_psamp_fm001.sh and the path segment in
    paths.py are written independently and must say the same thing."""
    src = open(os.path.join(SIM_DIR, "submit_2Mb_r3_x20_psamp_fm001.sh")).read()
    assert '_ntr_${N_CENTRAL}' in src
    assert 'n_central_traits=${N_CENTRAL} n_flank_gtex_traits=${N_FLANK}' in src
    assert paths.trait_count_segment({**BASE, "n_central_traits": 25}).endswith("ntr_25")
