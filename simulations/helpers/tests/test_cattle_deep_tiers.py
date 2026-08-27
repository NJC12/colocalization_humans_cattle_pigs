"""Resource tiers for the cattle deep-history chain, which had never been driven
through Snakemake before.

Two bugs were found the first time it was, and these tests pin both fixes.

1. `stage1_cattle_create_pop` declared `slurm_partition = "priority,long"` with a
   4-hour walltime. O2 refuses a sub-12-hour job on `long`:

       sbatch: error: Job not submitted: please submit jobs that are less than
       12 hours long to the short (or priority) partition.
       sbatch: error: Batch job submission failed: Invalid partition name specified

   so the rule could not submit at all, at any Q. It went unnoticed because every
   cattle deep history to date was launched as a standalone sbatch rather than
   through the rule.

2. The burn-in and selection rules asked for 30 days and 150 GB unconditionally.
   Correct at Q=0.01 -- the population is 1.7M and the `cattle_burnin_hedge` run
   has spent ~32 days to reach 6.6% of its target. Wrong by three orders of
   magnitude at Q=1, the value the round-3 deep history actually uses: MEASURED
   38m48s for the burn-in (sacct 48879426) and ~40 min for epochs 1-11.

The partition rule is the load-bearing invariant: `short` accepts <= 12 h and
`long` refuses anything shorter, so a walltime and a partition can never be
chosen independently.
"""

import os
import re

import pytest


SIM_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SMK = os.path.join(SIM_DIR, "rules", "stage1_cattle_baseline.smk")
SRC = open(SMK).read()

# The helpers are module-level in a .smk file, which is not importable as Python
# (it has rule blocks). Exec just the helper block, which is plain Python.
_NS = {}
_block = SRC[SRC.index("_CATTLE_FULL_SCALE_Q"):SRC.index("def _stage1_create_pop_inputs")]
exec(_block, _NS)

partition = _NS["_cattle_deep_partition"]
runtime = _NS["_cattle_deep_runtime"]
mem_mb = _NS["_cattle_deep_mem_mb"]

SHORT_CAP_MIN = 12 * 60


@pytest.mark.parametrize("q", [1, 1.0, 3, 10])
def test_full_scale_and_up_goes_to_short(q):
    assert partition({"Q_scaling": q}) == "short"


@pytest.mark.parametrize("q", [0.01, 0.1, 0.5])
def test_rescaled_populations_keep_the_long_partition(q):
    assert partition({"Q_scaling": q}) == "priority,long"


@pytest.mark.parametrize("q", [0.01, 0.5, 1, 10])
@pytest.mark.parametrize("hours", [4, 8])
def test_walltime_and_partition_are_never_chosen_independently(q, hours):
    """The invariant O2 enforces: <= 12 h must NOT go to `long`, and a job longer
    than 12 h must not go to `short`."""
    cfg = {"Q_scaling": q}
    part = partition(cfg)
    mins = runtime(cfg, hours)
    if "long" in part:
        assert mins > SHORT_CAP_MIN, (
            f"Q={q}: asks {mins} min on {part}; `long` refuses anything under "
            f"{SHORT_CAP_MIN} min and rejects the whole partition list")
    else:
        assert mins <= SHORT_CAP_MIN, f"Q={q}: asks {mins} min on {part}"


def test_create_pop_would_now_submit():
    """The exact combination that failed: 4 hours, at the Q the round-3 deep
    history runs."""
    cfg = {"Q_scaling": 1}
    assert partition(cfg) == "short"
    assert 4 * 60 <= SHORT_CAP_MIN


def test_the_rescaled_tier_is_still_four_weeks():
    """The hedge case must not be shortened by this change -- ~32 days has already
    gone into one such run."""
    assert runtime({"Q_scaling": 0.01}, 8) == 30 * 24 * 60
    assert mem_mb({"Q_scaling": 0.01}, 32000) == 150000


def test_the_full_scale_tier_has_headroom_over_the_measurement():
    """Measured 38m48s for the burn-in. 8 h is ~12x that and still fits `short`."""
    mins = runtime({"Q_scaling": 1}, 8)
    assert mins >= 12 * 39, "less than 12x the measured burn-in"
    assert mins <= SHORT_CAP_MIN


def test_no_rule_hardcodes_the_partition_any_more():
    """A literal "priority,long" left in a resources block is how bug 1 happened."""
    for m in re.finditer(r"slurm_partition\s*=\s*(.+?),\n", SRC):
        val = m.group(1).strip()
        assert not val.startswith('"'), (
            f"slurm_partition is a literal ({val}); it must be chosen from Q so "
            f"it cannot contradict the walltime")
