"""Make `helpers` importable from the tests without installing the package.

`simulations/` is the package root (see the sys.path comment in
`simulations/create_gwas_files_and_phenotypes.py`), so put it on the path.
"""

import os
import sys

_SIMULATIONS = os.path.dirname(  # simulations/
    os.path.dirname(  # simulations/helpers/
        os.path.dirname(os.path.abspath(__file__))  # simulations/helpers/tests/
    )
)
if _SIMULATIONS not in sys.path:
    sys.path.insert(0, _SIMULATIONS)
