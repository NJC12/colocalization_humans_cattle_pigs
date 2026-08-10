"""Relabel the power-sampled arms from causal_min_maf 0.01 to 0.

WHY THIS EXISTS

`submit_2Mb_r3_x20_psamp_fm001.sh` never passed `causal_min_maf`. It inherited
the old default of 0.01, and under `causal_sampling: power` the central pool
ignored the floor entirely -- the script's own header said so: "NOT a floor on
the causal pool under power sampling". The directories were nonetheless named
`maf_0.01`, because the floor named the directory whether or not it filtered
anything.

The floor now applies in both sampling schemes. Under that rule the content of
those directories is a `causal_min_maf: 0` run, and their names are wrong. This
script makes the names match the content. It does NOT recompute anything: the new
predicate at floor 0 is exactly the old no-floor predicate, which is what makes a
rename the right operation rather than a re-run.

WHAT MOVES

Four things, and they move together because `paths.causal_maf_segment()` stays
suppressed at 0.01 and therefore EMITS `.cmaf_0` at the new value:

  1. the stage-2 inner directory   gwas_G_gtex_T_maf_0.01  ->  ..._maf_0
  2. files whose names carry it    *_maf_0.01.tsv          ->  *_maf_0.tsv
  3. the run-tag directory under stage2..stage5, which gains `.cmaf_0`:
         hts_11.psamp_8000  ->  hts_11.cmaf_0.psamp_8000
  4. `causal_min_maf` inside stage2_params.txt and any run-level params file,
     with the uids recomputed -- otherwise the Snakefile's provenance guard
     rejects the directory it just renamed.

Two layouts are handled. The fetched local copies under `simulation_data/` are
FLAT (`<arm>/<ID>/<files>`) and have only 1, 2 and the arm directory's own
`cmaf01` token; the O2 workdirs are NESTED and have all four.

USAGE

Dry run is the default and prints every operation without performing any:

    python helpers/migrate_cmaf_psamp.py <root> [<root> ...]
    python helpers/migrate_cmaf_psamp.py <root> --apply

Run it on the local `simulation_data/` copies first, confirm the figure
notebook's build_manifest() still resolves those arms, then on O2.
"""

import argparse
import os
import re
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from helpers import params_record  # noqa: E402

OLD_MAF = "0.01"
NEW_MAF = "0"

# Only ever touch something that is demonstrably a power-sampled arm. Every guard
# below is a refusal to act on ambiguity: this script renames directories that
# other runs adopt by name, so acting on a near-miss is worse than doing nothing.
PSAMP_DIR = re.compile(r"(?:^|[._])psamp[_]?\d+(?:[._]|$)")
FILE_MAF = re.compile(r"_maf_" + re.escape(OLD_MAF) + r"(?=\.|$)")
STAGE2_INNER = re.compile(r"^(gwas_\d+_gtex_\d+)_maf_" + re.escape(OLD_MAF) + r"$")
RUN_TAG = re.compile(r"^(?P<base>.+?)(?P<samp>\.psamp_\d+)$")
ARM_CMAF = re.compile(r"(?:^|_)cmaf01(?=_|$)")


class Plan:
    """Collected operations, so a dry run prints exactly what --apply would do."""

    def __init__(self):
        self.renames = []      # (src, dst)
        self.rewrites = []     # (path, note)
        self.skipped = []      # (path, why)

    def rename(self, src, dst):
        self.renames.append((src, dst))

    def rewrite(self, path, note):
        self.rewrites.append((path, note))

    def skip(self, path, why):
        self.skipped.append((path, why))


def _rewrite_params(path, plan, apply_):
    """Set causal_min_maf to 0 in a params file and recompute its uids.

    Both file shapes are handled. The stage-2 sidecar is a flat mapping carrying
    `stage2_uid` at the top level; the run-level record carries a `_derived`
    block holding `stage2_uid`, `run_uid` and every resolved path. Without this
    the Snakefile's parse-time guard would compare the renamed directory against
    a record still claiming 0.01 and refuse it.
    """
    data = params_record.read(path)
    if data is None:
        plan.skip(path, "unreadable or not a mapping")
        return
    if str(data.get("causal_min_maf")) not in (OLD_MAF, "0.01"):
        plan.skip(path, f"causal_min_maf is {data.get('causal_min_maf')!r}, not {OLD_MAF}")
        return
    if data.get("causal_sampling") != "power":
        plan.skip(path, f"causal_sampling is {data.get('causal_sampling')!r}, not 'power'")
        return

    data["causal_min_maf"] = 0
    body = {k: v for k, v in data.items() if not k.startswith("_")}

    if "_derived" in data:
        # A run-level record: re-feedable as a config, so let params_record
        # recompute the whole derived block rather than patching strings.
        try:
            data["_derived"] = params_record.derived(body)
            note = "causal_min_maf -> 0, _derived recomputed"
        except Exception as exc:                                  # noqa: BLE001
            plan.skip(path, f"could not recompute _derived ({type(exc).__name__}: {exc})")
            return
    else:
        data["stage2_uid"] = params_record.uid(body, params_record.STAGE2_KEYS)
        note = f"causal_min_maf -> 0, stage2_uid -> {data['stage2_uid']}"

    plan.rewrite(path, note)
    if apply_:
        head = []
        with open(path) as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                head.append(line)
        with open(path, "w") as fh:
            fh.write("".join(head))
            fh.write("# causal_min_maf relabelled 0.01 -> 0 by helpers/migrate_cmaf_psamp.py.\n")
            fh.write("# The value never filtered this run's central pool; see that module.\n")
            yaml.safe_dump(data, fh, sort_keys=True, default_flow_style=False,
                           width=100, allow_unicode=True)


def _plan_dir_tree(root, plan, apply_):
    """Walk one root bottom-up, planning renames of files then directories.

    Bottom-up so a directory is renamed only after its contents have been dealt
    with; otherwise the paths collected for the contents go stale mid-walk.
    """
    for dirpath, dirnames, filenames in os.walk(root, topdown=False):
        for name in filenames:
            src = os.path.join(dirpath, name)
            if name.endswith("params.txt"):
                _rewrite_params(src, plan, apply_)
            if FILE_MAF.search(name):
                dst = os.path.join(dirpath, FILE_MAF.sub("_maf_" + NEW_MAF, name))
                plan.rename(src, dst)
                if apply_:
                    os.rename(src, dst)

        for name in dirnames:
            src = os.path.join(dirpath, name)
            m = STAGE2_INNER.match(name)
            if m:
                dst = os.path.join(dirpath, f"{m.group(1)}_maf_{NEW_MAF}")
            else:
                t = RUN_TAG.match(name)
                if not t or ".cmaf_" in name:
                    continue
                dst = os.path.join(dirpath, t.group("base") + ".cmaf_0" + t.group("samp"))
            plan.rename(src, dst)
            if apply_:
                os.rename(src, dst)


def _plan_arm_dir(root, plan, apply_):
    """The arm directory's own `cmaf01` token, for the flat local copies.

    figure2_revision2.ipynb's parse_arm() reads causal_min_maf out of this name
    (`cmaf<digits>` -> "0." + digits), so `cmaf0` decodes to 0. Renamed last, for
    the same bottom-up reason as above.
    """
    name = os.path.basename(root.rstrip("/"))
    if not ARM_CMAF.search(name):
        return root
    dst = os.path.join(os.path.dirname(root.rstrip("/")),
                       ARM_CMAF.sub(lambda m: m.group(0).replace("cmaf01", "cmaf0"), name))
    plan.rename(root, dst)
    if apply_:
        os.rename(root, dst)
    return dst


def _patch_provenance(root, plan, apply_):
    """PROVENANCE.txt states the old value in prose; leave a correction rather
    than silently disagreeing with the directory it describes."""
    path = os.path.join(root, "PROVENANCE.txt")
    if not os.path.exists(path):
        return
    text = open(path).read()
    if "migrate_cmaf_psamp" in text:
        plan.skip(path, "already carries the relabel note")
        return
    note = (
        "\n"
        "RELABELLED causal_min_maf 0.01 -> 0 (helpers/migrate_cmaf_psamp.py).\n"
        "Point 1 above was right about the behaviour and wrong about the name: the\n"
        "floor never gated this arm's central pool, so under the current rule --\n"
        "where a floor applies in both sampling schemes -- this arm IS a\n"
        "causal_min_maf = 0 run. Directory and filenames now say so. Nothing was\n"
        "recomputed: the floor-0 pool predicate is exactly the old no-floor one.\n"
        "To reproduce this arm, pass causal_min_maf=0 (now the script's default).\n"
    )
    plan.rewrite(path, "appended the relabel note")
    if apply_:
        with open(path, "a") as fh:
            fh.write(note)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("roots", nargs="+",
                    help="Arm directories (local flat copies) or run workdirs (O2).")
    ap.add_argument("--apply", action="store_true",
                    help="Perform the operations. Without it nothing is touched.")
    ap.add_argument("--force", action="store_true",
                    help="Skip the 'is this a power-sampled arm?' name check. Only for a "
                         "root whose own name lacks a psamp token but whose contents "
                         "carry .psamp_ run tags.")
    args = ap.parse_args()

    plan = Plan()
    for root in args.roots:
        root = root.rstrip("/")
        if not os.path.isdir(root):
            print(f"ERROR: not a directory: {root}", file=sys.stderr)
            return 1
        looks_psamp = bool(PSAMP_DIR.search(os.path.basename(root))) or any(
            ".psamp_" in d for _, dirs, _ in os.walk(root) for d in dirs)
        if not looks_psamp and not args.force:
            print(f"ERROR: {root} has no psamp token in its name and no .psamp_ run tag "
                  f"inside it. Refusing -- this renames directories other runs adopt by "
                  f"name. Pass --force if you are certain.", file=sys.stderr)
            return 1

        _patch_provenance(root, plan, args.apply)
        _plan_dir_tree(root, plan, args.apply)
        _plan_arm_dir(root, plan, args.apply)

    verb = "RENAMED" if args.apply else "would rename"
    for src, dst in plan.renames:
        print(f"{verb}: {src}\n      -> {dst}")
    verb = "REWROTE" if args.apply else "would rewrite"
    for path, note in plan.rewrites:
        print(f"{verb}: {path}  ({note})")
    for path, why in plan.skipped:
        print(f"skipped: {path}  ({why})")

    print(f"\n{len(plan.renames)} rename(s), {len(plan.rewrites)} rewrite(s), "
          f"{len(plan.skipped)} skipped.")
    if not args.apply:
        print("Dry run -- nothing was changed. Re-run with --apply.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
