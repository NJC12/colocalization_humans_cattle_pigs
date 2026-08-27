"""Make `helpers` importable from the tests without installing the package,
and share the one helper more than one test module needs.

`simulations/` is the package root (see the sys.path comment in
`simulations/create_gwas_files_and_phenotypes.py`), so put it on the path.
"""

import os
import re
import sys

import pytest

_SIMULATIONS = os.path.dirname(  # simulations/
    os.path.dirname(  # simulations/helpers/
        os.path.dirname(os.path.abspath(__file__))  # simulations/helpers/tests/
    )
)
if _SIMULATIONS not in sys.path:
    sys.path.insert(0, _SIMULATIONS)


@pytest.fixture
def guard_letters():
    """Read the letters out of a submit script's ``^(A|B|...)$`` category guard.

    Several scripts gate "this category builds its own stage 1 / stage 2" on an
    alternation of category letters, and more than one test module needs to
    assert that a letter is in one. Asserting on the SET rather than on the
    literal string matters: the assertions used to pin '^(H|I|J|K|L)$' verbatim,
    so every new category broke a test that had nothing to do with it -- a false
    alarm that teaches people to edit the assertion without reading it.

    ``marker`` is the line the guard protects, e.g. ``STAGE1_SRC["$ID"]=""``;
    the guard is the last ``if`` before it.
    """
    def _letters(sh, marker):
        before = sh[:sh.index(marker)].rsplit("if ", 1)[-1]
        m = re.search(r"\^\(([A-Z|]+)\)\$", before)
        assert m, f"no letter guard found before {marker!r}"
        return set(m.group(1).split("|"))
    return _letters
