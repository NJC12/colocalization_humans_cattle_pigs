"""Multi-directory input search for the simulation Snakefile.

Each pipeline stage may have several locations on disk where a previously
generated artifact might already exist. ``find_or_canonical`` walks those
locations in order and returns the first hit. If nothing is found it returns
the canonical path under the workdir, which Snakemake then schedules the
producing rule for.

If a hit is found whose embedded seed does not match the expected seed, raise
``WorkflowError`` (or ``ValueError`` outside Snakemake) -- per design we never
silently substitute a mismatched file as input to a 30-day downstream job.
"""

import os
import re
from pathlib import Path

try:
    from snakemake.exceptions import WorkflowError
except ImportError:
    WorkflowError = ValueError

DEFAULT_SEED_PATTERN = r"sd(\d+)|seed_(\d+)"


def _extract_seed(filename, pattern=DEFAULT_SEED_PATTERN):
    m = re.search(pattern, filename)
    if not m:
        return None
    return next((g for g in m.groups() if g is not None), None)


_SEED_GLOB_REPLACE = re.compile(r"sd\d+|seed_\d+")


def _seed_collision_check(filename, search_dirs, expected_seed):
    """Glob for files matching ``filename`` with the seed component replaced
    by ``*``. If any hit exists with a non-matching seed, raise WorkflowError.

    This catches the case where the user has e.g. ``...sd24.full.ts`` in a
    search dir but configured seed=99 -- without this check, exact lookup
    would silently miss the file and run stage 1 from scratch.
    """
    if expected_seed is None:
        return
    if not _SEED_GLOB_REPLACE.search(filename):
        return
    glob_pat = _SEED_GLOB_REPLACE.sub("sd*", filename)
    expected_str = str(expected_seed)
    from pathlib import Path
    for d in search_dirs or ():
        d = os.path.expanduser(str(d))
        if not os.path.isdir(d):
            continue
        for hit in Path(d).glob(glob_pat):
            got = _extract_seed(hit.name)
            if got is not None and got != expected_str:
                raise WorkflowError(
                    f"Search dir {d} contains {hit.name} with seed {got}, "
                    f"but config expects seed {expected_str} (would look for "
                    f"{filename}). Update the config or remove the search dir."
                )


def find_in_search_dirs(filename, search_dirs,
                        expected_seed=None, seed_pattern=DEFAULT_SEED_PATTERN,
                        verbose=False):
    """Return the first search-dir path containing ``filename``, or None.

    If a hit is found whose embedded seed differs from ``expected_seed``,
    raise ``WorkflowError``. If no exact hit is found but a same-shape file
    with a different seed exists in any search dir, also raise.
    """
    expected_str = None if expected_seed is None else str(expected_seed)

    for d in search_dirs or ():
        candidate = os.path.join(os.path.expanduser(str(d)), filename)
        if os.path.exists(candidate):
            if expected_str is not None:
                got = _extract_seed(filename, seed_pattern)
                if got is not None and got != expected_str:
                    raise WorkflowError(
                        f"Search hit {candidate} has seed {got}, "
                        f"but config expects seed {expected_str}. "
                        "Update the config or remove the search dir."
                    )
            if verbose:
                print(f"[search_dirs] {filename}: found at {candidate}")
            return candidate

    # Exact filename not found anywhere -- check for seed collisions before
    # giving up, so the user gets a clear error rather than a silent fresh run.
    _seed_collision_check(filename, search_dirs, expected_seed)

    if verbose:
        print(f"[search_dirs] {filename}: not found in any search dir")
    return None


def find_or_empty(filename, search_dirs, expected_seed=None,
                  seed_pattern=DEFAULT_SEED_PATTERN):
    """Snakemake-friendly wrapper: returns the search-dir path or ``[]``."""
    hit = find_in_search_dirs(filename, search_dirs,
                              expected_seed=expected_seed,
                              seed_pattern=seed_pattern)
    return hit if hit else []


def find_or_canonical(filename, search_dirs, canonical_dir,
                      expected_seed=None, seed_pattern=DEFAULT_SEED_PATTERN,
                      verbose=False):
    """Backwards-compatible wrapper: returns canonical path when no hit.

    Prefer ``find_in_search_dirs`` / ``find_or_empty`` in new code. This
    helper can introduce cyclic-dependency errors if used in a rule whose
    output equals the canonical path it would return.
    """
    hit = find_in_search_dirs(filename, search_dirs,
                              expected_seed=expected_seed,
                              seed_pattern=seed_pattern,
                              verbose=verbose)
    if hit is not None:
        return hit
    return os.path.join(os.path.expanduser(str(canonical_dir)), filename)


def find_first_match(pattern, search_dirs, canonical_dir, verbose=False):
    """Same as ``find_or_canonical`` but uses a glob pattern instead of an
    exact filename. Useful when legacy filenames vary (e.g. ``*_sd24.full.ts``
    matches ``revision_farm_selection_mult_100_gen_23_muts_26_bottlenecked_sd24.full.ts``).

    Returns (matched_path, is_real). If is_real is False, the canonical path
    is returned for the rule to produce.
    """
    for d in search_dirs or ():
        d = os.path.expanduser(str(d))
        if not os.path.isdir(d):
            continue
        matches = sorted(Path(d).glob(pattern))
        if matches:
            if verbose:
                print(f"[search_dirs] {pattern}: matched {matches[0]}")
            return str(matches[0]), True

    canonical = os.path.join(os.path.expanduser(str(canonical_dir)), pattern)
    return canonical, False
