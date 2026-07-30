"""Per-simulation parameter records.

Nothing on disk used to say what a simulation actually ran with, so
``helpers/summarize_coloc.py`` had to reverse-engineer it: read the ``.fmmaf``
sidecar, grep ``--ld-ctrl`` out of ``.snakemake/metadata/``, regex the
multipliers out of a directory name. This module writes the parameters down
instead.

Three files, one format:

``{workdir}/params/{run_tag}[.{output_tag}].params.txt``
    Written at Snakemake ``onstart``, so it exists as soon as the pipeline
    starts and survives a DAG that dies before stage 4.

``{stage4_dir}/{basename}[.{output_tag}].params.txt``
    Identical content, beside the final products and sharing their prefix. For
    any ``{prefix}.{gtex_cat}.enloc.{kind}.out``, strip
    ``.{gtex_cat}.enloc.{kind}.out`` and append ``.params.txt``.

``{stage2_inner}/stage2_params.txt``
    Written by ``create_gwas_files_and_phenotypes.py`` itself, and therefore the
    only one that can also record *outcomes* (pool size at each filter). Stage-2
    directories are symlink-adopted across output roots, so this is the only way
    an adopting run can know what the adopted phenotypes were -- which is why the
    Snakefile's stage-2 provenance guard reads it.

None of the three is a declared Snakemake output. Same reasoning as the
``.fmmaf`` sidecar in ``rules/common.smk``: runs that predate the file must stay
adoptable, and a metadata file should not append a job to every existing DAG.

The format is a flat YAML mapping saved with a ``.txt`` extension, so it is both
readable and ``yaml.safe_load``-able -- and re-feedable: ``snakemake --configfile
<params.txt>`` reproduces the run, since the Snakefile validates the keys it
requires and ignores the rest.
"""

import getpass
import hashlib
import json
import os
import socket
import subprocess
import sys
import tempfile
from datetime import datetime

import yaml

SIM_REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Config keys that determine stage-2 CONTENT -- the trait set and the
# phenotypes. Hashed into stage2_uid by both this module (from the resolved
# config) and create_gwas_files_and_phenotypes.py (from its own arguments, keyed
# by these same config names), so the stage-4 record and the stage-2 sidecar can
# be checked against each other after cross-root symlink adoption.
#
# pipeline, species, stage1_seed and recombination_rate also determine stage-2
# content, but only through the stage-1 tree sequence -- which the sidecar
# captures directly, by name and size, under `stage1`.
STAGE2_KEYS = (
    "causal_min_maf",
    "gwas_scaling",
    "gtex_scaling",
    "gtex_size",
    "n_samples",
    "n_central_traits",
    "n_flank_gtex_traits",
    "neutral_trait_vars",
    "already_includes_neutral",
    "L",
    "Q_scaling",
    "handoff_ticks",
    "deep_Q_scaling",
    "stage2_seed",
)

# Keys that change the ANALYSIS of a fixed stage-2 output: fine-mapping, the GWAS
# and colocalization. Two runs differing only in these share stage 2 and are what
# `output_tag` exists to keep apart.
ANALYSIS_KEYS = (
    "fm_min_maf",
    "min_maf",
    "ld_ctrl",
    "dapg_window",
    "fastenloc_prior",
    "fastenloc_total_snps",
    "output_tag",
    "skip_categories",
    "selection_cap",
    "selection_seed",
    "num_replicates",
)

_HEADER = """\
# Resolved simulation parameters -- written by helpers/params_record.py.
# NOT a Snakemake rule output (see the module docstring for why).
#
# Records the config AFTER the Snakefile's defaults and every --config override,
# so this is what ran, not what a YAML file says now. Valid YAML despite the
# .txt extension: yaml.safe_load() it, or feed it back with
# `snakemake --configfile <this file>`.
#
# THREE INDEPENDENT MAF FLOORS -- each governs exactly one thing, and none is
# derived from another:
#   causal_min_maf  which variants may be selected as CAUSATIVE      (stage 2)
#   fm_min_maf      which variants are FINE-MAPPED (the DAP-G SBAMS) (stage 3)
#   min_maf         which variants enter the GWAS (plink2 --maf)     (stage 5)
#
# Two things deliberately absent rather than guessed at:
#   fastenloc_total_snps: 'auto' is resolved inside scripts/run_fastenloc.sh at
#     job time; the number it used is in logs/stage4_fastenloc_{gtex_cat}.log.
#   Per-attempt SLURM resources: see .snakemake/log/ and `sacct`.
"""


# ---------------------------------------------------------------- serialisation

def _canon(v):
    """Canonical form for hashing, so int 2000000 and float 2000000.0 agree.

    The Snakefile sees ``L: 2000000`` from YAML while the stage-2 script parses
    ``--length`` as a float. Without this the two sides would hash the same
    simulation differently.
    """
    if isinstance(v, bool):
        return v
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, (list, tuple)):
        return [_canon(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _canon(x) for k, x in sorted(v.items())}
    if v is None:
        return None
    return str(v)


def _plain(obj):
    """Coerce to types yaml.safe_dump accepts.

    Snakemake's config is a dict subclass and may hold values that safe_dump
    refuses; a provenance file must never be the thing that breaks a run.
    """
    if isinstance(obj, bool) or obj is None or isinstance(obj, (int, float, str)):
        return obj
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_plain(v) for v in obj]
    return str(obj)


def _dump(mapping):
    return yaml.safe_dump(_plain(mapping), sort_keys=True, default_flow_style=False,
                          width=100, allow_unicode=True)


def uid(values, keys):
    """Short stable hash of ``values`` restricted to ``keys``.

    Missing keys hash as null, so adding a key to STAGE2_KEYS does not silently
    make two records incomparable -- it changes both sides' uid together.
    """
    payload = {k: _canon(values.get(k)) for k in keys}
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha1(blob.encode()).hexdigest()[:10]


# ------------------------------------------------------------------- the blocks

def _git(*args):
    try:
        out = subprocess.run(["git", "-C", SIM_REPO_DIR, *args],
                             capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def meta(argv=None):
    """Who/when/what-code. Never raises: the repo dir may be a read-only export."""
    try:
        import snakemake
        smk_version = getattr(snakemake, "__version__", None)
    except Exception:
        smk_version = None
    status = _git("status", "--porcelain")
    return {
        "written": datetime.now().astimezone().isoformat(timespec="seconds"),
        "host": socket.gethostname(),
        "user": getpass.getuser(),
        "python": sys.version.split()[0],
        "snakemake": smk_version,
        "sim_repo_dir": SIM_REPO_DIR,
        "git_commit": _git("rev-parse", "--short", "HEAD"),
        "git_dirty": None if status is None else bool(status),
        "argv": list(argv) if argv else None,
    }


def derived(cfg):
    """Values computed from the config rather than set in it."""
    from helpers import paths

    gtex_cats = paths.stage2_gtex_categories(cfg)
    return {
        "sim_id": os.path.basename(paths.workdir(cfg).rstrip("/")),
        "stage1_file": paths.stage1_full_filename(cfg),
        "stage2_run_tag": paths.stage2_run_tag(cfg),
        "stage2_inner": paths.stage2_inner(cfg),
        "stage2_params_file": paths.stage2_params_file(cfg),
        "stage2_categories": list(paths.stage2_categories(cfg)),
        "causal_maf_segment": paths.causal_maf_segment(cfg) or "",
        "output_tag": paths.output_tag(cfg),
        "stage3_outputs_subdir": paths.stage3_outputs_subdir_name(cfg),
        "stage4_prefixes": {g: paths.stage4_prefix(cfg, g) for g in gtex_cats},
        # The exact files this record describes.
        "stage4_outputs": [o for g in gtex_cats for o in paths.stage4_outputs(cfg, g)],
        "stage5_prefixes": {c: paths.stage5_prefix(cfg, c)
                            for c in paths.stage2_categories(cfg)},
        # Resolved from fastenloc_prior by scripts/run_fastenloc.sh.
        "fastenloc_flag": ("--coloc_default_prior"
                           if cfg.get("fastenloc_prior", "coloc_default") == "coloc_default"
                           else "-total_variants"),
        # Causal variants are confined to [5e5, L-5e5] by the stage-2 script.
        "trait_window": [5e5, float(cfg["L"]) - 5e5],
        "stage2_uid": uid(cfg, STAGE2_KEYS),
        "run_uid": uid(cfg, STAGE2_KEYS + ANALYSIS_KEYS),
    }


def render(cfg, extra_derived=None, argv=None):
    """The file body: header comment, ``_meta``, ``_derived``, then the config."""
    d = derived(cfg)
    if extra_derived:
        d.update(extra_derived)
    # Skip underscore keys so re-feeding a params file as --configfile (which
    # lands _meta/_derived in config) cannot emit duplicate top-level keys.
    plain_cfg = {k: cfg[k] for k in sorted(cfg) if not str(k).startswith("_")}
    return "".join([_HEADER,
                    _dump({"_meta": meta(argv)}),
                    _dump({"_derived": d}),
                    _dump(plain_cfg)])


# ---------------------------------------------------------------------- writing

def write_atomic(path, text):
    """Write via a temp file in the same directory, then ``os.replace``.

    The Snakefile's ``previous_workdirs`` adoption symlinks *every* file under
    ``stage4/<run_tag>/``, so a ``*.params.txt`` here can be a symlink pointing
    into another workdir. A plain ``open(path, 'w')`` would silently overwrite
    that other run's provenance; ``os.replace`` replaces the symlink itself.
    """
    d = os.path.dirname(path) or "."
    os.makedirs(d, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".params.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def read(path):
    """Parse a params file, or return None if it is missing or unreadable."""
    try:
        with open(path) as fh:
            data = yaml.safe_load(fh)
    except (OSError, yaml.YAMLError):
        return None
    return data if isinstance(data, dict) else None


def _disagreements(requested, on_disk, keys):
    out = []
    for k in keys:
        if k not in on_disk:
            continue                      # written by an older key set; not a conflict
        if _canon(on_disk.get(k)) != _canon(requested.get(k)):
            out.append((k, on_disk.get(k), requested.get(k)))
    return out


def stage2_disagreements(cfg, path):
    """Stage-2-determining keys where an existing stage-2 dir disagrees with cfg.

    Returns ``[(key, on_disk, requested), ...]`` -- empty when the stage-2 output
    at ``path``'s directory was produced with the parameters this run wants.
    """
    on_disk = read(path)
    if on_disk is None:
        return []
    return _disagreements(cfg, on_disk, STAGE2_KEYS)


def write_run_params(cfg, argv=None, strict=False):
    """Write both copies of the run record. Returns the paths written.

    If a record already exists at the stage-4 location and disagrees on any
    ANALYSIS key, the old one is archived as ``.params.<utc>.txt`` and a warning
    is printed -- that is the signature of two analysis variants sharing one
    output_tag. ``strict`` turns the warning into an error. Stage-2 disagreements
    are caught earlier and more cheaply by the Snakefile's parse-time guard.
    """
    from helpers import paths

    text = render(cfg, argv=argv)
    primary = paths.stage4_params_file(cfg)

    old = read(primary)
    if old is not None:
        clash = _disagreements(cfg, old, ANALYSIS_KEYS)
        if clash:
            detail = "; ".join(f"{k}: on disk {o!r} != requested {r!r}" for k, o, r in clash)
            msg = (f"{primary} already describes a run with different analysis "
                   f"parameters ({detail}). Two variants are sharing one "
                   f"output_tag, so their stage-3/4 outputs are colliding.")
            if strict:
                raise ValueError(msg)
            stamp = datetime.now().astimezone().strftime("%Y%m%dT%H%M%S%z")
            archive = primary[: -len(".txt")] + f".{stamp}.txt"
            try:
                write_atomic(archive, _dump({"_archived_from": primary}) + yaml.safe_dump(
                    _plain(old), sort_keys=True, default_flow_style=False))
                print(f"[params] WARNING: {msg}\n[params] previous record archived to {archive}")
            except OSError as e:
                print(f"[params] WARNING: {msg}\n[params] could not archive it: {e}")

    written = []
    for path in (paths.workdir_params_file(cfg), primary):
        write_atomic(path, text)
        written.append(path)
    return written


def write_stage2_params(path, values, pools, stage1=None):
    """Write the stage-2 sidecar. Called by create_gwas_files_and_phenotypes.py.

    ``values`` holds the stage-2-determining parameters keyed by their CONFIG
    names (so ``stage2_uid`` matches the run record's), ``pools`` the resulting
    pool sizes, ``stage1`` the identity of the input tree sequence.
    """
    body = dict(values)
    body["stage2_uid"] = uid(values, STAGE2_KEYS)
    body["_meta"] = meta(sys.argv)
    body["_pools"] = pools
    if stage1 is not None:
        body["stage1"] = stage1
    header = (
        "# Stage-2 parameters and outcomes, written by\n"
        "# create_gwas_files_and_phenotypes.py. Travels with this directory when it\n"
        "# is symlink-adopted into another output root, which is what lets the\n"
        "# Snakefile verify that adopted phenotypes were built under the causal MAF\n"
        "# floor the adopting run actually wants.\n"
        "#\n"
        "# _pools records how many variants survived each filter -- the quantity that\n"
        "# changes when causal_min_maf moves.\n"
    )
    write_atomic(path, header + _dump(body))
    return path
