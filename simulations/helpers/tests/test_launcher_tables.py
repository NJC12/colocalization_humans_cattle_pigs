"""The four launchers each carry their own copy of the category dispatch tables.

`seed_of`, `config_of`, `stage1_donor_of` and `category_extra` are pure functions
of the category letter (and `$CELL`), and each of the four round-3 launchers has
its own hand-maintained copy. They have already drifted: only
`submit_2Mb_r3_cmaf_fm001.sh` knows K and L at all, and
`submit_2Mb_r3_cmaf01_control.sh` has no `category_extra` whatsoever.

Equality across the four is therefore NOT the property to assert -- it fails on
the first run and tells you nothing. What matters is:

1. Where two launchers both define a letter, they must define it the SAME way.
   Disagreement means the same category name produces two different simulations
   depending on which script you invoked.
2. A letter that shares another category's config MUST have an override in that
   launcher's `category_extra`. This is the dangerous one. K and L map to A's and
   E's config files; they are only different simulations because
   `category_extra` adds `synthetic_dfe_effects=True`. A launcher that knows K in
   `config_of` and not in `category_extra` runs A's simulation, writes it to K's
   directory, and nothing downstream looks wrong.

These are source-parsing tests for the same reason `guard_letters` in conftest.py
is: the tables live in bash, and the alternative is discovering the drift from a
figure.
"""

import os
import re

import pytest


SIM_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

LAUNCHERS = [
    "submit_2Mb_r3_cmaf_fm001.sh",
    "submit_2Mb_r3_x20_psamp_fm001.sh",
    "submit_2Mb_r3_cmaf_replicates.sh",
    "submit_2Mb_r3_cmaf01_control.sh",
]

# Categories whose config file is another category's. Each MUST be given its
# distinguishing override by any launcher that can dispatch it.
SHARED_CONFIG_CATEGORIES = {
    "B": "neutral_trait_vars=True",
    "K": "synthetic_dfe_effects=True",
    "L": "synthetic_dfe_effects=True",
}


def _read(name):
    with open(os.path.join(SIM_DIR, name)) as fh:
        return fh.read()


def _case_table(sh, func):
    """{letter: emitted expression} for a `case "${1:0:1}" in ... esac` function.

    Handles the two shapes these tables use: one letter per line, and several
    letters sharing a line (`H) echo "8${n}" ;; I) echo "9${n}" ;;`), plus
    alternations (`K|L) echo ...`).
    """
    m = re.search(rf"^{func}\(\)\s*\{{(.*?)^\}}", sh, re.S | re.M)
    if not m:
        return None
    body = m.group(1)
    # Strip comments so a letter mentioned in prose is not mistaken for a branch.
    body = "\n".join(l for l in body.splitlines() if not l.lstrip().startswith("#"))
    table = {}
    for letters, expr in re.findall(r"([A-Z](?:\|[A-Z])*)\)\s*(.*?);;", body, re.S):
        for letter in letters.split("|"):
            table[letter] = " ".join(expr.split())
    return table or None


ALL = {name: _read(name) for name in LAUNCHERS}


def test_every_launcher_has_a_seed_table():
    for name in LAUNCHERS:
        assert _case_table(ALL[name], "seed_of"), name


@pytest.mark.parametrize("func", ["seed_of", "config_of", "stage1_donor_of"])
def test_launchers_agree_wherever_they_overlap(func):
    """Not equality -- agreement on the intersection."""
    tables = {n: (_case_table(ALL[n], func) or {}) for n in LAUNCHERS}
    letters = set().union(*(set(t) for t in tables.values()))
    disagreements = []
    for letter in sorted(letters):
        seen = {n: t[letter] for n, t in tables.items() if letter in t}
        if len(set(seen.values())) > 1:
            disagreements.append((letter, seen))
    assert not disagreements, (
        f"{func} disagrees between launchers: " + "; ".join(
            f"{l}: " + ", ".join(f"{n}={e!r}" for n, e in s.items())
            for l, s in disagreements))


@pytest.mark.parametrize("name", LAUNCHERS)
def test_a_shared_config_category_always_gets_its_override(name):
    """The silent-wrong-answer case: K dispatched without
    synthetic_dfe_effects=True is category A wearing K's name."""
    sh = ALL[name]
    configs = _case_table(sh, "config_of") or {}
    extra = _case_table(sh, "category_extra") or {}
    for letter, override in SHARED_CONFIG_CATEGORIES.items():
        if letter not in configs:
            continue  # this launcher cannot dispatch it at all; fine
        assert letter in extra, (
            f"{name} can dispatch category {letter} (config_of has it) but its "
            f"category_extra does not, so it would run the config's OWN category "
            f"and write the result under {letter}'s name")
        key = override.split("=")[0]
        assert key in extra[letter], (
            f"{name}'s category_extra for {letter} does not set {key}: "
            f"{extra[letter]!r}")


def test_a_category_that_adopts_a_donors_stage1_never_builds_its_own():
    """O, P (and, in the publication launcher, K and L) reuse a parent's tree
    sequence at the parent's seed. If such a category ever runs its own stage 1
    the pairing silently becomes two independent simulations."""
    for name in LAUNCHERS:
        donors = _case_table(ALL[name], "stage1_donor_of")
        if not donors:
            continue
        seeds = _case_table(ALL[name], "seed_of") or {}
        for letter, donor_expr in donors.items():
            m = re.search(r"([A-Z])\$\{1:1\}", donor_expr)
            if not m:
                continue  # the identity branch
            parent = m.group(1)
            assert letter in seeds and parent in seeds, (name, letter, parent)
            # Same seed band as the parent, or the seed-strict stage-1 lookup
            # could not find the donor's file in the first place.
            assert seeds[letter] == seeds[parent], (
                f"{name}: {letter} adopts {parent}'s stage 1 but their seed "
                f"expressions differ ({seeds[letter]!r} vs {seeds[parent]!r})")
