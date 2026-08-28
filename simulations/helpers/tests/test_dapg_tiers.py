"""DAP-G walltime tiers, and the branch the publication run actually uses.

`_dapg_runtime` picks a tier from whether the fine-mapping floor materially
thins the candidate set (`_fm_thins_candidates`, floor 0.005). The publication
run sets fm_min_maf = 0.001, which is BELOW that floor, so it takes the `else`
branch -- the one whose numbers were never measured.

68 stage3_dapg_locus jobs hit walltime at exactly 30:00 on that branch, all
gtex-side, every one succeeding on attempt 2. The 30 came from hgtex max 5:19
and cgtex max 6:29 measured at fm_min_maf = 0.01, where the filter thins the
human candidate set 8.7x. On this branch it does not thin at all, so the
candidate count is the unfiltered ~2847 and DAP-G's cost is at least linear in
it.

The retries hid the mis-sizing: the arms completed, nothing failed, and the only
trace was TIMEOUT rows in sacct. These tests pin the two properties that matter
-- that the branch the publication floors select is the unthinned one, and that
its tiers clear the observed 30-minute failures on the FIRST attempt.
"""

import os
import re

import pytest


SIM_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SMK = os.path.join(SIM_DIR, "rules", "common.smk")
SRC = open(SMK).read()

PUBLICATION_FM_MIN_MAF = 0.001
OBSERVED_TIMEOUT_MIN = 30


def _thinning_floor():
    m = re.search(r"_FM_THINNING_FLOOR\s*=\s*([0-9.]+)", SRC)
    assert m, "could not find _FM_THINNING_FLOOR"
    return float(m.group(1))


def _tiers():
    """{(thins, species, cat_kind): base_minutes} parsed from _dapg_runtime."""
    body = SRC[SRC.index("def _dapg_runtime"):SRC.index("def _dapg_mem_mb_base")]
    pairs = re.findall(r"base = (\d+) if wc\.cat\.endswith\(\"gwas\"\) else (\d+)", body)
    assert len(pairs) == 3, f"expected 3 tiers, found {len(pairs)}: {pairs}"
    return [(int(g), int(t)) for g, t in pairs]


def test_the_publication_floor_selects_the_unthinned_branch():
    """fm_min_maf 0.001 is below the thinning floor, so the publication run takes
    the branch whose numbers were never measured. If this ever stops being true
    the tier reasoning below no longer applies."""
    assert PUBLICATION_FM_MIN_MAF < _thinning_floor()


def test_three_tiers_are_present():
    assert len(_tiers()) == 3


def test_the_unthinned_tier_clears_the_observed_timeouts_on_attempt_one():
    """The `else` branch is the last of the three. Both of its bases must exceed
    the 30 minutes that 68 real jobs were killed at."""
    gwas_base, gtex_base = _tiers()[-1]
    assert gtex_base > OBSERVED_TIMEOUT_MIN, (
        f"unthinned gtex tier is {gtex_base} min; 68 jobs were killed at "
        f"{OBSERVED_TIMEOUT_MIN} min on this exact branch")
    assert gwas_base > OBSERVED_TIMEOUT_MIN


def test_the_unthinned_tier_is_never_smaller_than_a_thinned_one():
    """Thinning can only make DAP-G cheaper, so the unthinned branch must ask for
    at least as much as either thinned branch. Inverting this is how the gtex
    tier inherited a thinned-case number in the first place."""
    tiers = _tiers()
    unthinned_gwas, unthinned_gtex = tiers[-1]
    for gwas_base, gtex_base in tiers[:-1]:
        assert unthinned_gwas >= gwas_base
        assert unthinned_gtex >= gtex_base


def test_runtime_is_attempt_scaled():
    """A retry at an unchanged walltime is killed again; the 68 timeouts only
    recovered because of this."""
    body = SRC[SRC.index("def _dapg_runtime"):SRC.index("def _dapg_mem_mb_base")]
    assert "base * attempt" in body
