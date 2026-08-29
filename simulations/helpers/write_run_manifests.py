"""Build ARMS.tsv and RUNS.tsv for the publication simulation set.

WHY THIS IS PYTHON AND NOT AWK

`RUNS.tsv` carries two columns that must be *exactly* what the pipeline will use,
not a bash reconstruction of them:

  stage1_file      the launcher's pre-flight checks that this file exists before
                   submitting. That check is the ONLY guard on human stage-1
                   adoption -- `search_dirs.DEFAULT_SEED_PATTERN` matches
                   `sd(\\d+)` and `seed_(\\d+)`, and `hts_11.ts` matches neither,
                   so `_extract_seed` returns None and the seed-collision check
                   never fires. A wrong human tree sequence is adopted silently.
  stage2_run_tag   the sole namespace for stages 2-5. summarize_coloc and the
                   publish script both locate stage3/4/5 through it.

Both come from `helpers/paths.py`, so they are produced by the same code that
will consume them. Rebuilding either by string concatenation would work until a
segment was added, at which point the manifest and the pipeline would disagree
and nothing would say so.

Usage:
    python helpers/write_run_manifests.py --out-dir <publication root> \\
        [--arms causal_maf001_paired ...] [--reps 1 2 3 4 5] [--cell g5t20] \\
        [--scratch-root ...] [--publish-root ...]

Writes <out-dir>/ARMS.tsv and <out-dir>/RUNS.tsv. Idempotent.
"""

import argparse
import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from helpers import paths  # noqa: E402


HERE = os.path.dirname(os.path.abspath(__file__))
SIM_DIR = os.path.dirname(HERE)
CATEGORIES_TSV = os.path.join(HERE, "publication_categories.tsv")
ARMS_TSV = os.path.join(HERE, "publication_arms.tsv")
DEEP_HISTORIES_TSV = os.path.join(HERE, "publication_deep_histories.tsv")

RUNS_COLUMNS = [
    "arm", "run_dir", "category", "species", "replicate", "letter", "seed",
    "pipeline", "config", "stage1_donor", "category_extra",
    "stage1_file", "stage2_run_tag", "workdir", "publishdir",
    # Appended, not inserted: every reader indexes by header name
    # (submit_publication.sh, publish_to_data2.sh, summarize_coloc.runs_from_manifest),
    # so appending is free and reordering would not be.
    "cattle_baseline_seed",   # "-" for human
    "deep_history_handoff",   # the ep7 this run's stage 1 resumed from; "-" for human
    # cattle_sel only. rules/stage1_cattle_sel.smk requires the marks file AND the
    # .full.ts; a present full file with a missing marks file does not fail, it
    # falls through to a multi-hour SLiM run inside a controller sized for
    # symlinks. The launcher cannot pre-flight what the manifest does not name.
    "stage1_marks_file",
]


def read_tsv(path):
    """Rows of a #-commented TSV as dicts, keyed by the first non-comment line."""
    header, rows = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                continue
            if len(fields) != len(header):
                raise SystemExit(
                    f"{path}: row has {len(fields)} fields, header has "
                    f"{len(header)}: {line!r}")
            rows.append(dict(zip(header, fields)))
    if header is None:
        raise SystemExit(f"{path}: no header row found")
    return header, rows


def coerce(value):
    """TSV text -> the type the config would hold, so stage2_uid can match.

    params_record._canon already normalises int/float, but True/False and the
    numeric floors have to arrive as the right Python type or the uid this
    manifest implies is not the uid the run computes.
    """
    if value in ("True", "False"):
        return value == "True"
    if value == "-":
        return None
    try:
        f = float(value)
    except ValueError:
        return value
    return int(f) if f.is_integer() and "." not in value and "e" not in value.lower() else f


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--cell", default="g5t20")
    ap.add_argument("--arms", nargs="*", default=None)
    # Default from the deep-history table, NOT a literal. publish_to_data2.sh
    # regenerates the manifests with no --reps, so a literal here silently
    # rewrote the published RUNS.tsv down to whatever it said.
    ap.add_argument("--reps", nargs="*", type=int, default=None)
    ap.add_argument("--categories", nargs="*", default=None,
                    help="Category LETTERS. Default: every row of the table.")
    ap.add_argument("--scratch-root", required=True,
                    help="Where the runs actually execute.")
    ap.add_argument("--publish-root", required=True,
                    help="Where the publish script copies them (= --out-dir on O2).")
    args = ap.parse_args()

    _, dh_rows = read_tsv(DEEP_HISTORIES_TSV)
    deep = {int(r["replicate"]): r for r in dh_rows}
    if args.reps is None:
        args.reps = sorted(deep)

    arm_header, arm_rows = read_tsv(ARMS_TSV)
    cat_header, cat_rows = read_tsv(CATEGORIES_TSV)

    arms = {r["arm"]: r for r in arm_rows}
    cats = {r["letter"]: r for r in cat_rows}

    # The pairing invariant, checked here rather than trusted: a category that
    # adopts a donor's stage 1 must sit in the donor's seed band, or the
    # seed-strict lookup will not find the donor's file and a fresh genetic
    # simulation starts instead of an adoption.
    for letter, c in cats.items():
        donor = c["stage1_donor"]
        if donor != "self":
            if donor not in cats:
                raise SystemExit(f"category {letter}: unknown stage1_donor {donor!r}")
            if c["seed_prefix"] != cats[donor]["seed_prefix"]:
                raise SystemExit(
                    f"category {letter} adopts {donor}'s stage 1 but sits in seed "
                    f"band {c['seed_prefix']}N against {donor}'s "
                    f"{cats[donor]['seed_prefix']}N. Stage-1 lookup is seed-strict, "
                    f"so this would not adopt -- it would start a fresh simulation.")

    want_arms = args.arms or [r["arm"] for r in arm_rows]
    want_cats = args.categories or [r["letter"] for r in cat_rows]
    for a in want_arms:
        if a not in arms:
            raise SystemExit(f"unknown arm {a!r}; known: {sorted(arms)}")
    for c in want_cats:
        if c not in cats:
            raise SystemExit(f"unknown category {c!r}; known: {sorted(cats)}")

    os.makedirs(args.out_dir, exist_ok=True)

    # ARMS.tsv is the arms table verbatim, minus its comments -- it is already
    # one row per arm and one column per config key.
    arms_out = os.path.join(args.out_dir, "ARMS.tsv")
    with open(arms_out, "w") as fh:
        fh.write("\t".join(arm_header) + "\n")
        for r in arm_rows:
            if r["arm"] in want_arms:
                fh.write("\t".join(r[k] for k in arm_header) + "\n")

    runs = []
    for arm in want_arms:
        arm_cfg = {k: coerce(v) for k, v in arms[arm].items()
                   if k not in ("arm", "description")}
        for letter in want_cats:
            c = cats[letter]
            cfg_rel = c["config"].replace("${CELL}", args.cell)
            cfg_path = os.path.join(SIM_DIR, cfg_rel)
            if not os.path.exists(cfg_path):
                raise SystemExit(f"category {letter}: config not found: {cfg_path}")
            with open(cfg_path) as fh:
                base = yaml.safe_load(fh)

            extra = c["extra_config"]
            extra_cfg = {}
            if extra and extra != "-":
                for tok in extra.split():
                    k, _, v = tok.partition("=")
                    extra_cfg[k] = coerce(v)

            for rep in args.reps:
                seed = int(f"{c['seed_prefix']}{rep}")
                run_dir = f"{c['dir_name']}_rep{rep}"
                cfg = {**base, **arm_cfg, **extra_cfg,
                       "stage1_seed": seed, "stage2_seed": seed}

                # MUST be set before stage1_full_filename/stage2_run_tag are
                # called: it is in E/L's name as cb_{cb} and, since the segment
                # landed, in F/G's as _cb{cb}. Without this, cfg inherits the
                # YAML's 20250303 and the manifest names a file no run produces.
                cb = dh = marks = "-"
                if c["species"] == "cattle":
                    row = deep.get(rep)
                    if row is None:
                        raise SystemExit(
                            f"replicate {rep} has no row in "
                            f"{os.path.basename(DEEP_HISTORIES_TSV)}; a cattle run "
                            f"without one falls back to the config's 20250303")
                    cb = int(row["cattle_baseline_seed"])
                    dh = row["handoff_file"]
                    cfg["cattle_baseline_seed"] = cb
                    if c["pipeline"] == "cattle_sel":
                        marks = paths.stage1_cattle_sel_marks(cfg)
                runs.append({
                    "arm": arm,
                    "run_dir": run_dir,
                    "category": c["dir_name"],
                    "species": c["species"],
                    "replicate": rep,
                    "letter": letter,
                    "seed": seed,
                    "pipeline": c["pipeline"],
                    "config": cfg_rel,
                    "stage1_donor": c["stage1_donor"],
                    "category_extra": extra if extra != "-" else "",
                    "stage1_file": paths.stage1_full_filename(cfg),
                    "stage2_run_tag": paths.stage2_run_tag(cfg),
                    "workdir": os.path.join(args.scratch_root, arm, run_dir),
                    "publishdir": os.path.join(args.publish_root, arm, run_dir),
                    "cattle_baseline_seed": cb,
                    "deep_history_handoff": dh,
                    "stage1_marks_file": marks,
                })

    # Checked BEFORE anything is written: a manifest that is about to be
    # rejected should not be left on disk for someone to launch from.
    #
    # Two runs on one workdir silently interleave two DAGs.
    seen = {}
    for r in runs:
        if r["workdir"] in seen:
            raise SystemExit(f"workdir collision: {r['workdir']} claimed by "
                             f"{seen[r['workdir']]} and {r['arm']}/{r['run_dir']}")
        seen[r["workdir"]] = f"{r['arm']}/{r['run_dir']}"

    # Sharing a stage2_run_tag within an arm is legal for EXACTLY one reason: a
    # category that adopts another's stage-1 tree sequence runs at the donor's
    # seed, so A and K (and E and L) produce the same tag wherever the arm's
    # other segments do not separate them -- which is both power arms, since
    # require_gtex_partner is False for both members there. They are safe because
    # the tag is namespaced by workdir and those differ.
    #
    # Any OTHER pair sharing a tag is the failure this check exists for: the tag
    # is the sole namespace for stages 2-5, so two such runs write the same
    # stage2/, stage3/, stage4/ and stage5/ outputs. Two cattle deep histories at
    # one stage1_seed did precisely that before paths.cattle_baseline_segment.
    donor_of = {r["letter"]: r["stage1_donor"] for r in runs}
    groups = {}
    for r in runs:
        groups.setdefault((r["arm"], r["stage2_run_tag"]), []).append(r)
    for (arm, tag), members in sorted(groups.items()):
        if len(members) == 1:
            continue
        letters = {m["letter"] for m in members}
        reps = {m["replicate"] for m in members}
        files = {m["stage1_file"] for m in members}
        linked = any(donor_of.get(a) == b or donor_of.get(b) == a
                     for a in letters for b in letters if a != b)
        if not (linked and len(reps) == 1 and len(files) == 1):
            raise SystemExit(
                f"stage2_run_tag collision in {arm}: {tag} claimed by "
                + ", ".join(sorted(m["run_dir"] for m in members))
                + " -- these are not a donor pair adopting one stage-1 file")

    runs_out = os.path.join(args.out_dir, "RUNS.tsv")
    with open(runs_out, "w") as fh:
        fh.write("\t".join(RUNS_COLUMNS) + "\n")
        for r in runs:
            fh.write("\t".join(str(r[k]) for k in RUNS_COLUMNS) + "\n")

    # The reader-facing seed record. RUNS.tsv carries these columns among 18
    # others; this is the one file to open to answer "which seeds produced this?"
    seeds_out = os.path.join(args.out_dir, "SEEDS.tsv")
    seed_cols = ["arm", "run_dir", "letter", "replicate", "seed",
                 "cattle_baseline_seed", "deep_history_handoff", "stage1_file"]
    with open(seeds_out, "w") as fh:
        fh.write("\t".join(seed_cols) + "\n")
        for r in runs:
            fh.write("\t".join(str(r[k]) for k in seed_cols) + "\n")

    print(f"wrote {arms_out} ({len(want_arms)} arms)")
    print(f"wrote {runs_out} ({len(runs)} runs: {len(want_arms)} arms x "
          f"{len(want_cats)} categories x {len(args.reps)} replicates)")
    print(f"wrote {seeds_out} ({len(runs)} runs)")


if __name__ == "__main__":
    main()
