#!/usr/bin/env python
"""Build the cross-cell colocalization / fine-mapping summary table.

One row per (sim_category, gtex_n, gwas_mult, gtex_mult), pooling all stage-1
replicates found for that category.  Everything is discovered from the output
tree rather than hard-coded, so new replicates (A2, A3, ... E2, ...) and new
multiplier cells are picked up automatically.

Definitions
-----------
RCP / LCP     SIGNAL-level RCP / LCP from `*.enloc.sig.out`.  A sig.out signal
              is named `<gene>:<cluster>(@)` and its cluster id is the DAP-G
              cluster id (verified: Num_SNP and CPIP_qtl match DAP-G's
              member_snp and cluster_pip exactly), so each colocalizing signal
              can be traced back to a specific credible set.
correct match the signal cluster that CONTAINS the gene's defined causal
              variant, for a gene that has a matching GWAS trait.  Flank genes
              have no matching GWAS trait and therefore no correct match, so
              every call they produce is wrong by construction.
true positive a gene with at least one signal cluster over the cutoff where one
              of those clusters is the correct match.  Denominator is the
              number of defined GWAS traits.
false positive a gene with at least one signal cluster over the cutoff, none of
              which is the correct match.  This includes CAUSATIVE genes that
              colocalize through the wrong credible set -- right gene, wrong
              signal -- not just flank genes.  Reported as a false DISCOVERY
              rate, fp / (tp + fp), so the denominator is the number of
              positive calls.
              Genes with no cluster over the cutoff are neither.
auc_enloc_N   single-operating-point AUC at RCP > 0.N: (sensitivity +
              specificity) / 2, with sensitivity = TP / defined GWAS traits and
              specificity = 1 - FP / all genes evaluated.  A per-threshold
              summary that pairs with the tp/fp columns, NOT a threshold-free
              ROC-AUC.  0.5 is chance.
CPIP          PIP of the DAP-G signal cluster CONTAINING the true causal
              variant, over the causative loci.  A locus fine-mapped
              confidently to the wrong cluster counts as a failure.
credible set  a DAP-G independent signal cluster.
              *_cs_all        every cluster at every locus fine-mapped in that
                              panel (GWAS: the GWAS traits; GTEx: all genes,
                              causative and flank).  Unconditional.
              *_cs_causative  the single cluster containing the defined causal
                              variant, at the causative loci only -- the same
                              locus set as the CPIP columns.  Conditional on
                              the causal variant having been clustered at all.
fm_filtered   whether the genotypes handed to DAP-G were MAF-filtered.
              Read from the run's params file when present; otherwise detected by
              comparing the stage-3 sbams row count against the stage-2 variant
              table, filtered and unfiltered.
causal_min_maf the MAF floor the CAUSATIVE variants were drawn from -- a distinct
              knob from the fine-mapping floor (fm_filtered/fm_r2) and from the
              GWAS floor.  Part of the grouping key: the same cell rerun at a
              different causal floor has different phenotypes and must not pool.
n_causal_tested of the `n_gwas_traits` defined GWAS loci, how many had their
              causal variant fine-mapped at all -- i.e. its MAF cleared the
              fine-mapping floor.  Below `causal_min_maf` < `fm_min_maf` this is
              less than n_gwas_traits, and it is the denominator of the CPIP and
              *_cs_causative columns.  The colocalization columns use ALL loci: a
              locus whose causal variant went untested can still colocalize
              correctly, because the GWAS and GTEx traits there share that one
              variant and both credible sets are built from tags of it.  See
              classify().  What such a locus cannot tell you is whether the
              credible set contains the causal variant -- a fine-mapping question,
              not a colocalization one.
params_source `params` when the parameters were read from the run's params file,
              `inferred` when they were reverse-engineered from paths, row counts
              and Snakemake metadata (every run predating the params file).
"""
import argparse
import glob
import gzip
import os
import re
import sys
from statistics import median

try:
    import yaml
except ImportError:   # no PyYAML: params reading degrades to the legacy inference
    yaml = None

DEMOGRAPHY = {"A": "human", "B": "human", "C": "human", "D": "human",
              "E": "cattle", "F": "cattle", "G": "cattle"}

COLUMNS = ["sim_category", "demography", "gtex_n", "replicates",
           "causal_min_maf", "n_gwas_traits", "n_causal_tested",
           "fm_filtered", "fm_r2", "params_source", "gwas_mult", "gtex_mult",
           "enloc_tp_rcp50", "enloc_fp_rcp50", "enloc_tp_rcp90", "enloc_fp_rcp90",
           "enloc_tp_lcp50", "enloc_fp_lcp50", "enloc_tp_lcp90", "enloc_fp_lcp90",
           "auc_enloc_50", "auc_enloc_90",
           "gwas_cpip_50", "gtex_cpip_50", "gwas_cpip_90", "gtex_cpip_90",
           "gwas_mean_cs_all", "gwas_mean_cs_causative",
           "gtex_mean_cs_all", "gtex_mean_cs_causative",
           "gwas_median_cs_all", "gwas_median_cs_causative",
           "gtex_median_cs_all", "gtex_median_cs_causative"]


def parse_dapg(path):
    """-> cluster_of_snp, cluster_pip, cluster_size"""
    clus, cpip, csize = {}, {}, {}
    sec = None
    with open(path) as fh:
        for line in fh:
            t = line.strip()
            if t.startswith("Posterior inclusion probability"):
                sec = "pip"; continue
            if t.startswith("Independent association signal clusters"):
                sec = "clu"; continue
            if sec == "pip" and t.startswith("(("):
                f = t.split()
                if len(f) >= 5:
                    clus[f[1]] = int(f[4])
            elif sec == "clu" and t.startswith("{"):
                f = t.replace("{", " ").replace("}", " ").split()
                if len(f) >= 3:
                    cid = int(f[0]); csize[cid] = int(f[1]); cpip[cid] = float(f[2])
    return clus, cpip, csize


def read_sig_out(path):
    """Signal-level enloc output -> {gene: {cluster_id: (RCP, LCP)}}.

    Signal names look like `tr88703:1(@)` -- gene before ':', DAP-G cluster id
    after.  Clusters absent from sig.out never cleared any threshold and are
    simply missing here.
    """
    out = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.split()
            if len(f) < 7:
                continue
            name = f[0].split("(")[0]
            gene, _, cid = name.partition(":")
            try:
                cid = int(cid)
            except ValueError:
                continue
            out.setdefault(gene, {})[cid] = (float(f[5]), float(f[6]))
    return out


def classify(genes, idx, thr):
    """genes: iterable of (is_shared, causal_cluster_id, {cid: (rcp, lcp)}, causal_tested).

    Returns (n_true_positive, n_false_positive).  A gene contributes nothing
    unless at least one of its signal clusters clears `thr` on metric `idx`.

    The correct-signal test is by causal-variant identity, which requires that
    variant to have been fine-mapped at all.  When it was not (`causal_tested`
    False -- its MAF fell below fm_min_maf, so stage3_index_geno dropped it from
    the SBAMS), the test falls back to the GENE level: a shared gene that
    colocalizes counts as a true positive.

    That fallback is sound because of how the simulation is built, not as an
    approximation.  A shared central locus gives the GWAS trait and the GTEx
    trait the SAME causal variant, so there is exactly one true signal there.
    With that variant untested, both credible sets are assembled from tags of
    that one signal, and fastEnloc colocalizing them is a genuine detection --
    what is lost is only the ability to confirm the credible set CONTAINS the
    causal variant, which is a fine-mapping question (the CPIP and
    *_cs_causative columns) rather than a colocalization one.  Scoring those
    loci by identity would count every correct colocalization as a false
    positive and inflate the reported FDR.
    """
    tp = fp = 0
    for is_shared, causal_cid, sigs, causal_tested in genes:
        above = [cid for cid, v in sigs.items() if v[idx] > thr]
        if not above:
            continue
        correct = (causal_cid in above) if causal_tested else True
        if is_shared and correct:
            tp += 1
        else:
            fp += 1
    return tp, fp


def auc_point(genes, n_pos, idx, thr):
    """Single-operating-point AUC at threshold `thr` -- balanced accuracy.

    AUC of the ROC polygon through one measured point:
        (sensitivity + specificity) / 2
    with
        sensitivity = TP / (defined GWAS traits)
        specificity = 1 - FP / (all genes evaluated)

    The FP denominator is EVERY gene, not just the flank ones, because under
    the signal-level definition a causative gene can also emit a false call by
    colocalizing through the wrong credible set -- so every gene has an
    opportunity to be wrong.  0.5 is chance; below 0.5 means the calls carry
    less information than a coin flip at that cutoff.

    Note this is a per-threshold summary, not a threshold-free ROC-AUC: it
    answers "how good is the classifier AT this cutoff", which is what pairs
    with the tp/fp columns beside it.
    """
    n = len(genes)
    if not n_pos or not n:
        return float("nan")
    tp, fp = classify(genes, idx, thr)
    return (tp / n_pos + (1.0 - fp / n)) / 2.0


def header_traits(path):
    with open(path) as fh:
        return fh.readline().rstrip("\n").split("\t")[2:]


def n_individuals(path):
    with open(path) as fh:
        return sum(1 for _ in fh) - 1


def gz_rows(path):
    with gzip.open(path, "rt") as fh:
        return sum(1 for _ in fh)


def detect_fm_filtered(stage2_dir, stage3_dir, cat, min_maf):
    """True/False if resolvable, else None.

    Prefers the `.fmmaf` sidecar written by stage3_index_geno.  Runs predating
    that sidecar fall back to comparing the SBAMS row count against the stage-2
    variant table, filtered and unfiltered.
    """
    geno = os.path.join(stage3_dir, "geno.sbams.gz")
    sidecar = geno + ".fmmaf"
    if os.path.exists(sidecar):
        try:
            with open(sidecar) as fh:
                return float(fh.read().strip()) > 0
        except (OSError, ValueError):
            pass
    vf = glob.glob(os.path.join(stage2_dir, f"{cat}_vars_*.tsv"))
    if not (os.path.exists(geno) and vf):
        return None
    total = kept = 0
    with open(vf[0]) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        try:
            mi = cols.index("maf")
        except ValueError:
            return None
        for line in fh:
            total += 1
            if float(line.split("\t")[mi]) >= min_maf:
                kept += 1
    n = gz_rows(geno)
    if n == total:
        return False
    if n == kept:
        return True
    return None


def detect_ld_ctrl(rep_dir, cap=400):
    """The --ld-ctrl actually passed to dap-g, read from Snakemake's metadata.

    `.snakemake/metadata/<b64 path>` records the shell command that produced
    each output, so this is what ran rather than what a config says now.
    Returns a set of observed values (empty if the metadata is gone).
    """
    md = os.path.join(rep_dir, ".snakemake", "metadata")
    if not os.path.isdir(md):
        return set()
    vals, seen = set(), 0
    try:
        names = os.listdir(md)
    except OSError:
        return set()
    for name in names:
        if seen >= cap:
            break
        try:
            with open(os.path.join(md, name)) as fh:
                txt = fh.read()
        except (OSError, UnicodeDecodeError):
            continue
        seen += 1
        for m in re.finditer(r"--ld-ctrl\s+([0-9.]+)", txt):
            vals.add(float(m.group(1)))
    return vals


def read_params(rep, stage1_tag):
    """The run's params file as a dict, or None for runs that predate it.

    Both copies written by helpers/params_record.py have identical content, but
    they are NOT equally safe to locate.  The stage-4 copy lives in a directory
    named for the run tag, so it is unambiguous; the workdir copy is named
    `{run_tag}[.{output_tag}].params.txt`, and since the causal-MAF segment is also
    dot-joined, a prefix glob for `hts_11` happily matches
    `hts_11.cmaf_0.001.params.txt` -- i.e. a baseline run would read the variant's
    parameters.  So: stage-4 first, then the workdir copy only when its recorded
    run tag confirms it.
    """
    if yaml is None:
        return None

    def load(path):
        # Skip the .params.<timestamp>.txt archives left by a conflicting rerun.
        if re.search(r"\.params\.\d{8}T\d{6}", os.path.basename(path)):
            return None
        try:
            with open(path) as fh:
                data = yaml.safe_load(fh)
        except Exception:
            return None
        return data if isinstance(data, dict) else None

    for path in sorted(glob.glob(os.path.join(rep, "stage4", stage1_tag, "*.params.txt"))):
        data = load(path)
        if data is not None:
            return data

    exact = os.path.join(rep, "params", f"{stage1_tag}.params.txt")
    data = load(exact) if os.path.exists(exact) else None
    if data is not None:
        return data
    for path in sorted(glob.glob(os.path.join(rep, "params", f"{stage1_tag}.*.params.txt"))):
        data = load(path)
        if data is None:
            continue
        # Only trust a prefix match that says which run it belongs to.
        if (data.get("_derived") or {}).get("stage2_run_tag") == stage1_tag:
            return data
    return None


def causal_maf_by_trait(stage2_dir, cat):
    """{trait_id: MAF of that trait's defined causal variant} from stage 2.

    The vars table carries every variant's MAF in this panel plus a non-zero
    `beta` on the causal one, and a trait is named `tr<position>` -- so the causal
    variant's frequency in the very sample being fine-mapped is recoverable
    without touching the genotypes.
    """
    vf = glob.glob(os.path.join(stage2_dir, f"{cat}_vars_*.tsv"))
    if not vf:
        return {}
    out = {}
    with open(vf[0]) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        try:
            pi, mi, bi = cols.index("position"), cols.index("maf"), cols.index("beta")
        except ValueError:
            return {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                if float(f[bi]) == 0.0:
                    continue
                out["tr%d" % int(float(f[pi]))] = float(f[mi])
            except (ValueError, IndexError):
                continue
    return out


def tested_flag(maf_by_trait, trait, fm_floor):
    """Was this trait's causal variant in the set handed to DAP-G?

    True when no filter was applied (fm_floor 0) or when its MAF cleared the
    floor.  Unknown -- no vars table, or the causal variant is missing from it --
    also reads True, which is the legacy assumption and keeps pre-params runs on
    one row instead of inventing a split.
    """
    if not fm_floor:
        return True
    maf = maf_by_trait.get(trait)
    return True if maf is None else maf >= fm_floor


def pct(hits, n):
    return round(100.0 * hits / n, 1) if n else float("nan")


def _new_row():
    return dict(reps=set(), genes=[], n_gwas=0, n_causal_tested=0,
                gw_cpip=[], gt_cpip=[],
                gw_cs_all=[], gw_cs_caus=[], gt_cs_all=[], gt_cs_caus=[])


def collect(root):
    rows = {}
    for rep in sorted(glob.glob(os.path.join(root, "[A-G][0-9]*"))):
        cat_letter = os.path.basename(rep)[0]
        # EVERY stage-2 directory, not just the first. One workdir can legitimately
        # hold several -- a run that changed gwas/gtex scaling or the causal MAF
        # floor gets its own stage-2 dir and its own stage-3/4/5 subtree -- and
        # taking [0] silently reported one variant's numbers under whichever
        # directory glob happened to return first.
        for s2 in sorted(d for d in glob.glob(os.path.join(rep, "stage2", "*", "*"))
                         if os.path.isdir(d)):
            _collect_one(rows, rep, cat_letter, s2)
    return rows


def _collect_one(rows, rep, cat_letter, s2):
    stage1 = os.path.basename(os.path.dirname(s2))
    # The maf_ component is the CAUSAL floor: it is the only one of the three MAF
    # knobs that changes what stage 2 produces, so it is what names the directory.
    m = re.search(r"gwas_(\d+)_gtex_(\d+)_maf_([0-9.eE+-]+)", os.path.basename(s2))
    if not m:
        return
    gwas_mult, gtex_mult = int(m.group(1)), int(m.group(2))
    causal_maf = float(m.group(3))

    s3root = os.path.join(rep, "stage3", stage1)
    if not os.path.isdir(s3root):
        return
    cats = sorted(os.listdir(s3root))
    gwas_cat = next((c for c in cats if c.endswith("gwas")), None)
    gtex_cats = [c for c in cats if "gtex" in c]
    if not gwas_cat:
        return
    gwf = glob.glob(os.path.join(s2, f"{gwas_cat}_traits_*.tsv"))
    if not gwf:
        return
    gwas_traits = header_traits(gwf[0])

    # Parameters: read them when the run wrote them down, infer them otherwise.
    # Everything already on disk predates the params file, so the inference paths
    # below stay load-bearing rather than becoming dead code.
    params = read_params(rep, stage1)
    if params is not None:
        causal_maf = float(params.get("causal_min_maf", causal_maf))
        fm_floor = float(params.get("fm_min_maf", 0) or 0)
        ldc = params.get("ld_ctrl")
        ldc = None if ldc is None else float(ldc)
        source = "params"
    else:
        # min_maf and causal_min_maf were the same 0.01 in every pre-params run, so
        # the parsed directory value is also the right fine-mapping floor to test
        # row counts against.
        fm_floor = causal_maf
        ld_ctrl_seen = detect_ld_ctrl(rep)
        ldc = next(iter(ld_ctrl_seen)) if len(ld_ctrl_seen) == 1 else None
        source = "inferred"

    gw_filt = detect_fm_filtered(s2, os.path.join(s3root, gwas_cat), gwas_cat, fm_floor)
    if params is not None:
        gw_filt = fm_floor > 0

    # Was each locus's defined causal variant actually in the tested set? This
    # gates the FINE-MAPPING columns only. Colocalization is scored over every
    # locus -- see classify() for why an untested causal variant still permits a
    # genuine true positive.
    gwas_causal_maf = causal_maf_by_trait(s2, gwas_cat)

    # ---- GWAS side: independent of the eQTL panel ----
    gw_out = os.path.join(s3root, gwas_cat, "outputs")
    gw = dict(cpip=[], cs_all=[], cs_caus=[])
    for tr in gwas_traits:
        f = os.path.join(gw_out, f"{tr}.dapg.out")
        if not os.path.exists(f):
            continue
        clus, cpip, csize = parse_dapg(f)
        cid = clus.get("snp" + tr[2:], -1)
        # cs_all is unconditional: it describes every credible set the panel
        # produced, which is meaningful whether or not the causal variant is in
        # one. CPIP and cs_caus are about the causal variant specifically, so a
        # locus where it was never a candidate would contribute a guaranteed 0
        # and silently drag the mean down.
        gw["cs_all"].extend(csize.values())
        if not tested_flag(gwas_causal_maf, tr, fm_floor):
            continue
        gw["cpip"].append(cpip.get(cid, 0.0) if cid > 0 else 0.0)
        if cid > 0 and cid in csize:
            gw["cs_caus"].append(csize[cid])

    for gcat in gtex_cats:
        tf = glob.glob(os.path.join(s2, f"{gcat}_traits_*.tsv"))
        outdir = os.path.join(s3root, gcat, "outputs")
        sig = glob.glob(os.path.join(rep, "stage4", stage1, f"*.{gcat}.enloc.sig.out"))
        if not tf or not os.path.isdir(outdir) or not sig:
            continue
        gtex_traits = header_traits(tf[0])
        gtex_n = n_individuals(tf[0])
        sigs = read_sig_out(sig[0])
        gtex_causal_maf = causal_maf_by_trait(s2, gcat)

        # The fine-mapping settings and the causal floor are part of the GROUPING
        # key, not just reported: the same cell rerun at a different ld_ctrl, MAF
        # floor or causal floor is a separate experimental condition and must not
        # be pooled with the original.
        gt_filt = (fm_floor > 0 if params is not None
                   else detect_fm_filtered(s2, os.path.join(s3root, gcat), gcat, fm_floor))
        fmf = gw_filt if gw_filt is not None else gt_filt

        def row():
            key = (cat_letter, gtex_n, gwas_mult, gtex_mult, causal_maf,
                   fmf, ldc, source)
            r = rows.setdefault(key, _new_row())
            r["reps"].add(os.path.basename(rep))
            return r

        # Whether the causal variant was fine-mapped on the GTEx side, which is the
        # side whose cluster ids the correct-match test uses.
        def tested_of(tr):
            return tested_flag(gtex_causal_maf, tr, fm_floor)

        r = row()
        # Denominator for the true-positive rate: every defined GWAS trait. Loci
        # whose causal variant went untested stay IN it -- they can still
        # colocalize, so excluding them would understate the achievable rate.
        r["n_gwas"] += len(gwas_traits)
        r["n_causal_tested"] += sum(1 for tr in gwas_traits if tested_of(tr))
        r["gw_cpip"].extend(gw["cpip"])
        r["gw_cs_all"].extend(gw["cs_all"])
        r["gw_cs_caus"].extend(gw["cs_caus"])

        # One pass over every GTEx gene -- causative and flank alike. The flank
        # genes are needed here, not just for the FP count: under the signal-level
        # definition a causative gene that colocalizes through the wrong credible
        # set is itself a false positive, so the classification needs each gene's
        # causal cluster id.
        for tr in gtex_traits:
            f = os.path.join(outdir, f"{tr}.dapg.out")
            if not os.path.exists(f):
                continue
            clus, cpip, csize = parse_dapg(f)
            cid = clus.get("snp" + tr[2:], -1)
            is_shared = tr in gwas_traits
            tested = tested_of(tr)
            r["genes"].append((is_shared, cid, sigs.get(tr, {}), tested))
            r["gt_cs_all"].extend(csize.values())
            if is_shared and tested:
                r["gt_cpip"].append(cpip.get(cid, 0.0) if cid > 0 else 0.0)
                if cid > 0 and cid in csize:
                    r["gt_cs_caus"].append(csize[cid])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roots", nargs="+", required=True)
    ap.add_argument("--fm-r2", type=float, default=0.75,
                    help="fallback -ld_control, used only when a run's "
                         "Snakemake metadata is unavailable")
    ap.add_argument("-o", "--out", default="-")
    a = ap.parse_args()

    allrows = {}
    for root in a.roots:
        for k, v in collect(root).items():
            if k in allrows:
                for f in ("genes", "gw_cpip", "gt_cpip",
                          "gw_cs_all", "gw_cs_caus", "gt_cs_all", "gt_cs_caus"):
                    allrows[k][f].extend(v[f])
                allrows[k]["reps"] |= v["reps"]
                allrows[k]["n_gwas"] += v["n_gwas"]
            else:
                allrows[k] = v

    def mean(xs):
        return round(sum(xs) / len(xs), 2) if xs else "NA"

    def med(xs):
        return median(xs) if xs else "NA"

    fh = sys.stdout if a.out == "-" else open(a.out, "w")
    print("\t".join(COLUMNS), file=fh)
    for key, r in sorted(
            allrows.items(),
            # category, then causal floor (descending, so the historical 0.01 sorts
            # first), then the fine-mapping settings, then panel size and multipliers.
            key=lambda kv: (kv[0][0], -kv[0][4], bool(kv[0][5]),
                            kv[0][6] if kv[0][6] else 0, -kv[0][1],
                            kv[0][2], kv[0][3])):
        cat, gtex_n, gwm, gtm, causal_maf, fmf, ldc, source = key
        out = []
        for idx, thr in ((0, 0.5), (0, 0.9), (1, 0.5), (1, 0.9)):
            ntp, nfp = classify(r["genes"], idx, thr)
            out.append((pct(ntp, r["n_gwas"]), pct(nfp, ntp + nfp)))
        (r50, r90, l50, l90) = out
        fm_filtered = "NA" if fmf is None else str(fmf)
        # Measured from Snakemake's record of the dap-g command line; the
        # --fm-r2 CLI value is only a fallback when that metadata is gone.
        fm_r2 = a.fm_r2 if ldc is None else ldc
        vals = [
            cat, DEMOGRAPHY.get(cat, "?"), gtex_n, len(r["reps"]),
            causal_maf, r["n_gwas"], r["n_causal_tested"],
            fm_filtered, fm_r2, source, gwm, gtm,
            r50[0], r50[1], r90[0], r90[1], l50[0], l50[1], l90[0], l90[1],
            round(auc_point(r["genes"], r["n_gwas"], 0, 0.5), 3),
            round(auc_point(r["genes"], r["n_gwas"], 0, 0.9), 3),
            pct(sum(1 for x in r["gw_cpip"] if x > 0.5), len(r["gw_cpip"])),
            pct(sum(1 for x in r["gt_cpip"] if x > 0.5), len(r["gt_cpip"])),
            pct(sum(1 for x in r["gw_cpip"] if x > 0.9), len(r["gw_cpip"])),
            pct(sum(1 for x in r["gt_cpip"] if x > 0.9), len(r["gt_cpip"])),
            mean(r["gw_cs_all"]), mean(r["gw_cs_caus"]),
            mean(r["gt_cs_all"]), mean(r["gt_cs_caus"]),
            med(r["gw_cs_all"]), med(r["gw_cs_caus"]),
            med(r["gt_cs_all"]), med(r["gt_cs_caus"]),
        ]
        print("\t".join(str(v) for v in vals), file=fh)
    if fh is not sys.stdout:
        fh.close()


if __name__ == "__main__":
    main()
