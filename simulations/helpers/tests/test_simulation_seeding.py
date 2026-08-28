"""Every simulation call in the pipeline must be seeded.

Three separate unseeded calls were found while freezing the code for the
publication run, and each one meant a tree sequence could not be regenerated
from the seed it was named for:

  human_simulation_o2.py      engine.simulate(...) with no seed=, so stdpopsim
                              drew its own SLiM and recapitation seeds. Measured:
                              two runs at seed 11 gave 55,150 vs 54,900 nodes.
  farm_create_orig_pop_e2.py  msprime.sim_ancestry(...) with no random_seed, so
                              the coalescent was a fresh draw. Measured: two runs
                              at seed 20250303 gave 81,474 vs 81,230 nodes.
  farm_create_orig_pop_e2.py  msprime.sim_mutations(...) with no random_seed.

They were found one at a time, by re-running and comparing. That works but does
not scale, and the third was only found because the fix for the first two was
verified rather than assumed. This test is the systematic version: it walks the
AST of every pipeline script and fails on a simulation call that takes no seed.

A LOCAL wrapper named `simulate` that takes its seed positionally is fine and is
allowed for by name -- `human_neutral_simulation.simulate(seed, ...)` and its
cattle counterpart both do that, and their inner msprime calls are seeded.
"""

import ast
import glob
import os

import pytest


SIM_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# msprime / stdpopsim entry points that consume randomness.
SIM_CALLS = {"sim_ancestry", "sim_mutations", "simulate"}
SEED_KWARGS = {"random_seed", "seed"}

SCRIPTS = sorted(
    glob.glob(os.path.join(SIM_DIR, "*.py"))
    + glob.glob(os.path.join(SIM_DIR, "helpers", "*.py"))
)


def _local_function_names(tree):
    return {n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}


def _unseeded_calls(path):
    """[(lineno, name)] for simulation calls that pass no seed."""
    with open(path) as fh:
        tree = ast.parse(fh.read(), filename=path)
    local = _local_function_names(tree)
    out = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if isinstance(node.func, ast.Attribute):
            name, qualified = node.func.attr, True
        elif isinstance(node.func, ast.Name):
            name, qualified = node.func.id, False
        else:
            continue
        if name not in SIM_CALLS:
            continue
        # A bare call to a function defined in this same file is a local wrapper,
        # which is free to take its seed positionally.
        if not qualified and name in local:
            continue
        if {k.arg for k in node.keywords if k.arg} & SEED_KWARGS:
            continue
        out.append((node.lineno, name))
    return out


@pytest.mark.parametrize("path", SCRIPTS, ids=lambda p: os.path.basename(p))
def test_every_simulation_call_is_seeded(path):
    bad = _unseeded_calls(path)
    assert not bad, (
        f"{os.path.relpath(path, SIM_DIR)} has unseeded simulation calls at "
        + ", ".join(f"line {ln} ({name})" for ln, name in bad)
        + ". An unseeded call makes the output unregenerable from its seed, and "
          "nothing downstream notices -- the file is produced, named for a seed "
          "that did not determine it.")


def test_the_audit_would_catch_a_regression():
    """The test is only worth having if it actually fires."""
    src = "import msprime\nts = msprime.sim_ancestry(samples=10)\n"
    tree = ast.parse(src)
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and getattr(node.func, "attr", None) in SIM_CALLS:
            if not ({k.arg for k in node.keywords if k.arg} & SEED_KWARGS):
                found.append(node.lineno)
    assert found == [2]


def test_the_two_known_local_wrappers_are_not_flagged():
    """Regression guard on the exemption itself: if it ever stops matching, these
    two files start failing for the wrong reason and the exemption gets widened."""
    for name in ("human_neutral_simulation.py", "cattle_neutral_simulation.py"):
        path = os.path.join(SIM_DIR, name)
        if os.path.exists(path):
            assert _unseeded_calls(path) == [], name
