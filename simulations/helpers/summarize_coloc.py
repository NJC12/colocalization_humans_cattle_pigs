#!/usr/bin/env python
"""Build the cross-cell colocalization / fine-mapping summary table.

One row per (sim_category, gtex_n, gwas_mult, gtex_mult), pooling all stage-1
replicates found for that category.  Everything is discovered from the output
tree rather than hard-coded, so new replicates (A2, A3, ... E2, ...) and new
multiplier cells are picked up automatically.

Definitions
-----------
COLOCALIZATION is scored between PAIRS OF TRAITS, and never by which credible set
              carried the signal.  A GWAS trait and a GTEx trait share a name
              (`tr<position>`) exactly when they share a causal variant, so
              `tr887241` <-> `tr887241` is the one correct pairing for that GWAS
              trait.  Colocalizing through a tagging variant instead of the causal
              one is still a genuine detection -- and is the INTENDED route
              whenever `causal_min_maf` < the fine-mapping floor, which drops the
              causal variant from the tested set on purpose.

              The pairing comes from `*.enloc.sig.out`, whose signals are named
              `<eqtl_gene>:<eqtl_cluster>(@)<gwas_trait>:<gwas_cluster>`
              (controller.cc:1436).  Cluster ids on both sides are discarded and a
              GTEx trait's several clusters collapse to one entry by max, so every
              count below is over DISTINCT TRAITS.

              Per GWAS trait W, at one metric and cutoff:
                self_above   some signal with eqtl_gene == W over the cutoff
                other_traits distinct eqtl_gene != W over the cutoff
enloc_pow_*   COUNT of GWAS traits with self_above, however many other_traits they
              also hit.  Denominator: n_gwas_traits (50 per replicate).
enloc_fp_*    COUNT of GWAS traits with NOT self_above and other_traits non-empty
              -- it colocalized, just with the wrong trait.  Denominator for a
              false-discovery rate: enloc_pow_* + enloc_fp_*.  The two outcomes are
              mutually exclusive, and a trait clearing nothing is neither.
enloc_hit_med_gtex_sigs_*
              over the enloc_pow_* traits only, the MEDIAN number of distinct GTEx
              traits colocalized with, counting the correct one -- so 1 is the floor
              and a median of 1 means the typical true positive is specific.  NA
              when there are no power traits.
auc_enloc_N   single-operating-point AUC at RCP > 0.N: (sensitivity +
              specificity) / 2, with sensitivity = pow / n_gwas_traits and
              specificity = 1 - fp / n_gwas_traits (both outcomes are indexed on
              the GWAS trait, so n_gwas_traits bounds either).  A per-threshold
              summary that pairs with the pow/fp columns, NOT a threshold-free
              ROC-AUC, so it is sensitive to where the cutoff falls for a given
              condition.  0.5 is chance.

              _rcp* columns test RCP only and _lcp* columns LCP only.  LCP >= RCP
              always, so each lcp count is >= its rcp counterpart.

All enloc_pow_*, enloc_fp_*, gwas_cpip_* and gtex_cpip_* columns are RAW COUNTS,
not percentages; their denominators are the n_* columns described below.
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
              variant and both credible sets are built from tags of it.  What such
              a locus cannot tell you is whether the credible set contains the
              causal variant -- a fine-mapping question, not a colocalization one.
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
except ImportError:
    yaml = None


def _parse_params(text):
    """Minimal reader for the params files, for environments without PyYAML.

    Silently falling back to the legacy inference would be worse than useless
    here: it would assign the run a fine-mapping floor taken from the stage-2
    directory name, which for a lowered causal floor is the WRONG number, and the
    table would look fine while being wrong. The files are a flat mapping of
    scalars plus the `_meta` / `_derived` blocks, so parsing the part this script
    needs takes no real YAML support.

    Returns top-level scalars plus `_derived`'s scalar children.
    """
    def scalar(v):
        v = v.strip()
        if len(v) >= 2 and v[0] == v[-1] and v[0] in "'\"":
            return v[1:-1]
        low = v.lower()
        if low in ("true", "false"):
            return low == "true"
        if low in ("null", "~", ""):
            return None
        for cast in (int, float):
            try:
                return cast(v)
            except ValueError:
                pass
        return v

    out, block = {}, None
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indented = raw[0] in " \t"
        line = raw.strip()
        if ":" not in line:
            continue
        key, _, val = line.partition(":")
        key, val = key.strip(), val.strip()
        if not indented:
            # A key with no value opens a nested block; we only follow _derived.
            block = key if not val else None
            if val:
                out[key] = scalar(val)
            elif key == "_derived":
                out["_derived"] = {}
        elif block == "_derived" and val and not val.endswith(("{", "[")):
            out["_derived"][key] = scalar(val)
    return out

DEMOGRAPHY = {"A": "human", "B": "human", "C": "human", "D": "human",
              "E": "cattle", "F": "cattle", "G": "cattle",
              # H and I: neutral GENEALOGY (a pure coalescent, no selection in
              # stage 1 at all) with effect sizes drawn from the truncated DFE --
              # H under the human demography, I under the cattle one. They draw
              # from the SAME distribution, so H-vs-I isolates the demography.
              # B is the other neutral model -- a genome shaped by selection, with
              # the effects moved onto neutral variants. Keep them apart.
              # J is category A with the other branch of OutOfAfrica_2T12 sampled
              # (population: AFR). Same species, same model, same selection --
              # A vs J isolates ancestry, i.e. the LD and allele-frequency
              # structure colocalization power actually depends on.
              # K and L are the BACKGROUND SELECTION pair: A's and E's genomes
              # (forward runs under the DFE, so the genealogy carries background
              # selection) with H/I's effect model on top -- causal loci drawn
              # from the strictly NEUTRAL variants, betas from the truncated DFE.
              # K - H isolates the genealogy; A - K isolates the effect
              # assignment. They are the only arms where those two are separable.
              # M and N are the Wang-2014 Finnish founder pair -- a human
              # demography that is not in stdpopsim's catalog at all. M samples
              # the FIN deme, N the NFE one, and the two configs are otherwise
              # identical, so M - N is the founder event alone and N - A is the
              # model swap alone. Both run Q_scaling 3, not 10: the FIN deme is
              # too small at Q=10 to yield 9,000 individuals.
              "H": "human", "I": "cattle", "J": "human",
              "K": "human", "L": "cattle",
              "M": "human", "N": "human",
              # O and P are A's and E's genomes with the neutral class thinned
              # so pi halves -- same seeds, same causal loci, less background
              # variation. A - O and E - P isolate variant density.
              "O": "human", "P": "cattle"}

COLUMNS = [# "NA" for every root without a RUNS.tsv, i.e. everything written
           # before the publication tree.
           "arm",
           "sim_category", "demography", "gtex_n", "replicates",
           "causal_min_maf", "n_gwas_traits", "n_causal_tested",
           "n_gwas_cpip_tested", "n_gtex_cpip_tested",
           "fm_filtered", "fm_r2", "params_source", "gwas_mult", "gtex_mult",
           # How the causal set was drawn. Part of the grouping key: two arms
           # that differ only here are different experiments, and pooling them
           # halves nothing and doubles the denominator.
           "causal_sampling", "sampling_gwas_n", "sampling_power_plateau",
           "n_central_traits", "require_gtex_partner",
           # Raw COUNTS, not percentages. Denominators: power over
           # n_gwas_traits; the false-discovery rate over (pow + fp).
           "enloc_pow_rcp50", "enloc_fp_rcp50", "enloc_pow_rcp90", "enloc_fp_rcp90",
           "enloc_pow_lcp50", "enloc_fp_lcp50", "enloc_pow_lcp90", "enloc_fp_lcp90",
           "enloc_hit_med_gtex_sigs_rcp50", "enloc_hit_med_gtex_sigs_rcp90",
           "enloc_hit_med_gtex_sigs_lcp50", "enloc_hit_med_gtex_sigs_lcp90",
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


# Signals dropped because fastEnloc attributed them to no GWAS locus, plus the
# merged-locus tally. Reported on stderr by main() rather than silently swallowed:
# if the sig.out naming ever changes, every colocalization would land here and the
# table would quietly read as "nothing colocalized anywhere".
_SIG_STATS = {"rows": 0, "no_gwas": 0, "no_gwas_over_half": 0,
              "merged": 0, "merged_over_half": 0}

# Long-format per-GWAS-trait records, enabled by --dump-traits. None means off, and
# the extra work (a second gwas_partners pass per panel) is skipped entirely. The
# dump exists so questions the aggregate table cannot answer -- power split by
# whether the causal variant was testable, sensitivity to the merged-locus policy
# -- can be answered without another cluster round trip.
_TRAIT_DUMP = None


def read_sig_out(path):
    """Signal-level enloc output -> [(eqtl_gene, gwas_trait, RCP, LCP)].

    fastEnloc names a signal `<eqtl_gene>:<eqtl_cluster>(@)<gwas_locus>` --
    controller.cc:1436 builds it as `eqtl_vec[i].id + "(@)" + snp2gwas_locus[snp]`
    -- e.g. `tr544663:1(@)tr544663:1`.  The suffix names the GWAS locus this eQTL
    signal colocalized with, which is what makes trait-level scoring possible.

    Cluster ids on both sides are deliberately discarded: colocalization is scored
    per trait pair, so which credible set carried the signal does not matter.

    Two suffix forms beyond the simple one, both seen in real output:

      EMPTY -- `tr16568:2(@)` with nothing after.  The signal's SNPs fall in no
        GWAS locus at all, so it cannot be attributed to any GWAS trait and is
        dropped.  Counted in _SIG_STATS; harmless in practice because these carry
        RCP at the 1e-4 floor, but a loud problem if one ever clears a cutoff.

      MERGED -- `tr806894:2(@)tr634882:1_tr717151:1`, an underscore-joined list.
        fastEnloc merges GWAS loci that share SNPs, so one eQTL signal can overlap
        the signal regions of SEVERAL GWAS traits at once.  Every named trait is
        emitted as its own row, i.e. the signal is credited to all of them.

        Taking only the first, as an earlier version did, would be arbitrary --
        the order is SNP order, not evidence order -- and would silently convert a
        correct pairing into a miss whenever the matching trait happened to be
        listed second.  Dropping merged rows outright would understate power and
        false positives together.  Crediting every constituent is the honest
        reading of "this eQTL signal colocalized with a region containing these
        GWAS traits", and it is symmetric: a merged locus can just as easily hand
        a trait a false positive as a true one.  Trait names are `tr<position>`
        and never contain an underscore, so the split is unambiguous.
    """
    rows = []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.split()
            if len(f) < 7:
                continue
            eq, _, gw = f[0].partition("(@)")
            gene = eq.split(":")[0]
            rcp, lcp = float(f[5]), float(f[6])
            _SIG_STATS["rows"] += 1
            if not gw:
                _SIG_STATS["no_gwas"] += 1
                if max(rcp, lcp) > 0.5:
                    _SIG_STATS["no_gwas_over_half"] += 1
                continue
            loci = [x.split(":")[0] for x in gw.split("_") if x]
            if len(loci) > 1:
                _SIG_STATS["merged"] += 1
                if max(rcp, lcp) > 0.5:
                    _SIG_STATS["merged_over_half"] += 1
            for gwas in loci:
                rows.append((gene, gwas, rcp, lcp, len(loci) > 1))
    return rows


def gwas_partners(sig_rows, drop_merged=False):
    """-> {gwas_trait: {eqtl_gene: (max RCP, max LCP)}}

    Collapses a GTEx trait's several clusters into one entry by taking the max,
    which is what makes the downstream counts "distinct traits" rather than
    "signals".

    `drop_merged` discards every signal whose GWAS attribution spanned merged
    loci, giving the strict counterpart to the default credit-all policy.  Used to
    measure how much the reported numbers rest on the ambiguous attributions
    rather than to produce them -- see read_sig_out.
    """
    out = {}
    for gene, gwas, rcp, lcp, merged in sig_rows:
        if drop_merged and merged:
            continue
        d = out.setdefault(gwas, {})
        prev = d.get(gene, (0.0, 0.0))
        d[gene] = (max(prev[0], rcp), max(prev[1], lcp))
    return out


def _above(val_pair, idx, thr):
    return val_pair[idx] > thr


def power_and_fp(hits, idx, thr):
    """hits: iterable of (self_rcp, self_lcp, {other_gene: (rcp, lcp)}), one entry
    per defined GWAS trait -- including traits that colocalized with nothing.

    power = the GWAS trait colocalizes with its OWN GTEx trait, however many other
            GTEx traits it also hits.
    fp    = it does NOT colocalize with its own GTEx trait but DOES colocalize with
            at least one other.
    Neither when nothing clears the cutoff.  The `elif` is load-bearing: the two
    outcomes are mutually exclusive, so a correct hit is never also counted wrong.
    """
    pw = fp = 0
    for self_rcp, self_lcp, others in hits:
        if (self_rcp, self_lcp)[idx] > thr:
            pw += 1
        elif any(_above(v, idx, thr) for v in others.values()):
            fp += 1
    return pw, fp


def hit_med(hits, idx, thr):
    """Median number of DISTINCT GTEx traits colocalized with, over the power
    traits only, counting the correct trait itself -- so the floor is 1 and a
    median of 1 means the typical true positive is specific.
    """
    counts = []
    for self_rcp, self_lcp, others in hits:
        if (self_rcp, self_lcp)[idx] <= thr:
            continue
        counts.append(1 + sum(1 for v in others.values() if _above(v, idx, thr)))
    return median(counts) if counts else "NA"


def auc_point(hits, n_pos, idx, thr):
    """Single-operating-point AUC at threshold `thr` -- balanced accuracy.

    AUC of the ROC polygon through one measured point:
        (sensitivity + specificity) / 2
    with
        sensitivity = power / (defined GWAS traits)
        specificity = 1 - fp / (defined GWAS traits)

    Both denominators are the GWAS-trait count because both outcomes are now
    indexed on the GWAS trait: each of the n_pos traits can contribute at most
    one power or one fp, so n_pos is the ceiling for either.  0.5 is chance.

    Note this is a per-threshold summary, not a threshold-free ROC-AUC: it
    answers "how good is the classifier AT this cutoff", which is what pairs
    with the pow/fp columns beside it, and it is therefore sensitive to where
    the cutoff happens to fall for a given condition.
    """
    if not n_pos:
        return float("nan")
    pw, fp = power_and_fp(hits, idx, thr)
    return (pw / n_pos + (1.0 - fp / n_pos)) / 2.0


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
    def load(path):
        # Skip the .params.<timestamp>.txt archives left by a conflicting rerun.
        if re.search(r"\.params\.\d{8}T\d{6}", os.path.basename(path)):
            return None
        try:
            with open(path) as fh:
                text = fh.read()
        except OSError:
            return None
        try:
            data = yaml.safe_load(text) if yaml is not None else _parse_params(text)
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


def shared_partner_by_trait(stage2_dir, gwas_cat):
    """{gwas trait id: True if a GTEx trait was defined at the SAME causal variant}.

    Read from the `{cat}_trait_partners_*.tsv` sidecar stage 2 writes. Without it
    the pairing has to be inferred from trait-name equality, which cannot tell
    "this GWAS locus has no eQTL to find" apart from "colocalization failed" --
    and those are opposite results.

    The distinction only bites when the GTEx causal set is topped up rather than
    intersected: under `causal_sampling: power`, under the drawn-DFE arms (H, I
    and the background-selection pair K, L) in both sampling schemes, and under
    any arm that sets `require_gtex_partner: False` explicitly. A locus with no partner cannot produce a
    true colocalization, so any RCP it earns against some OTHER GTEx trait
    (`n_other_50` / `n_other_90` in the trait dump) is a FALSE POSITIVE, and the
    unpartnered loci are the only clean denominator for that rate.

    Returns {} when the sidecar is absent, which is every run predating it; the
    dump then carries `shared` as None rather than silently claiming False.
    """
    pf = glob.glob(os.path.join(stage2_dir, f"{gwas_cat}_trait_partners_*.tsv"))
    if not pf:
        return {}
    out = {}
    with open(pf[0]) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        try:
            ti, si = cols.index("gwas_trait"), cols.index("shared")
        except ValueError:
            return {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                out[f[ti]] = f[si].strip().lower() in ("true", "1")
            except IndexError:
                continue
    return out


def tested_flag(maf_by_trait, trait, fm_floor):
    """Was this trait's causal variant in the set handed to DAP-G?

    Three cases, and the difference between the last two matters:
      - no filter applied (fm_floor 0)            -> True
      - no vars table at all (empty mapping)      -> True, the legacy assumption,
        so runs whose stage-2 tables were never fetched behave as before
      - table present but this trait absent       -> False.  The vars table lists
        every variant segregating in the panel, so a causal variant missing from
        it was not a fine-mapping candidate here.  This is not hypothetical: in a
        500-individual GTEx subsample several causal variants drawn at
        causal_min_maf=0.001 are monomorphic and drop out entirely.  Counting such
        a locus in the CPIP denominator would contribute a guaranteed zero and drag
        the fine-mapping rates down for a reason that is not fine-mapping failure.

    Gates the FINE-MAPPING columns only.  The colocalization columns are scored
    per trait pair and ignore this flag entirely.
    """
    if not fm_floor:
        return True
    if not maf_by_trait:
        return True
    maf = maf_by_trait.get(trait)
    return False if maf is None else maf >= fm_floor


def _new_row():
    # gwas_hits: one (self_rcp, self_lcp, others) record per replicate x GWAS
    # trait. Records are appended rather than keyed by trait name so that pooling
    # replicates cannot collide when two stage-1 seeds happen to place a causal
    # variant at the same position.
    return dict(reps=set(), gwas_hits=[], n_gwas=0, n_causal_tested=0,
                gw_cpip=[], gt_cpip=[],
                gw_cs_all=[], gw_cs_caus=[], gt_cs_all=[], gt_cs_caus=[])


def runs_from_manifest(root):
    """[(rep_dir, letter, arm)] from RUNS.tsv, or None when there is no manifest.

    The publication tree names its runs interpretably --
    `human_background_selection_rep3`, not `K3` -- so the `[A-Z][0-9]*` glob below
    finds nothing there and the species letter is not in the path at all. RUNS.tsv
    carries both, plus the arm.

    Looked for beside the root AND one level up, because the layout is arm-major
    (`<publication root>/<arm>/<run>`) and `--roots` is normally given the arm
    directories, not their parent.

    Returns None rather than an empty list when absent, so the caller can tell
    "no manifest here, use the glob" from "manifest says this root is empty".
    """
    for cand in (os.path.join(root, "RUNS.tsv"),
                 os.path.join(os.path.dirname(os.path.abspath(root)), "RUNS.tsv")):
        if not os.path.isfile(cand):
            continue
        rows, header = [], None
        with open(cand) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                f = line.split("\t")
                if header is None:
                    header = {name: i for i, name in enumerate(f)}
                    continue
                rows.append({k: f[i] for k, i in header.items() if i < len(f)})

        # ONE arm's rows, when the root IS an arm directory.
        #
        # The run_dir column is arm-relative and the same 30 names repeat under
        # every arm, so `<root>/<run_dir>` resolves for rows belonging to OTHER
        # arms just as happily as for this one. Without this filter each run is
        # visited once per arm in the manifest and emitted under every arm's
        # label -- 4x the rows, three quarters of them mislabelled. The metrics
        # survive it (the grouping key carries the real parameters, read from
        # each run's own params file, so physically different arms still do not
        # pool) but the `arm` column becomes fiction, and that column is what
        # the figures select on.
        #
        # Keyed on the manifest's own arm values, not on the assumption that a
        # root is an arm: a root that is not one is left unfiltered, which is
        # what keeps a flat single-arm tree working.
        arms = {r.get("arm") for r in rows if r.get("arm")}
        here = os.path.basename(os.path.abspath(root))
        if here in arms:
            rows = [r for r in rows if r.get("arm") == here]

        out = []
        for row in rows:
            rep = os.path.join(root, row["run_dir"])
            if not os.path.isdir(rep):
                # A run that has not landed yet.
                continue
            out.append((rep, row.get("letter", "?"), row.get("arm")))
        if out:
            return sorted(out)
    return None


def collect(root):
    rows = {}
    found = runs_from_manifest(root)
    if found is None:
        # [A-Z], not [A-J]: the glob is the thing a new category is most likely to be
        # missed by, and it fails by silently reporting nothing rather than erroring.
        # A directory only matches if a digit follows the letter, so sibling dirs
        # like OLD_r2_75/ and simulation_summaries/ are still excluded.
        #
        # Every root written before the publication tree is discovered this way, so
        # this stays load-bearing rather than becoming dead code.
        found = [(rep, os.path.basename(rep)[0], None)
                 for rep in sorted(glob.glob(os.path.join(root, "[A-Z][0-9]*")))]
    for rep, cat_letter, arm in found:
        # EVERY stage-2 directory, not just the first. One workdir can legitimately
        # hold several -- a run that changed gwas/gtex scaling or the causal MAF
        # floor gets its own stage-2 dir and its own stage-3/4/5 subtree -- and
        # taking [0] silently reported one variant's numbers under whichever
        # directory glob happened to return first.
        for s2 in sorted(d for d in glob.glob(os.path.join(rep, "stage2", "*", "*"))
                         if os.path.isdir(d)):
            _collect_one(rows, rep, cat_letter, s2, arm=arm)
    return rows


def _collect_one(rows, rep, cat_letter, s2, arm=None):
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
    # locus: it is a trait-pair question, so an untested causal variant still
    # permits a genuine detection through a tagging variant.
    gwas_causal_maf = causal_maf_by_trait(s2, gwas_cat)

    # Which GWAS loci actually have a GTEx trait at the same causal variant.
    # Only meaningful where the GTEx set is topped up rather than intersected;
    # empty for every run predating the sidecar, and for those `shared` stays None.
    gwas_shared = shared_partner_by_trait(s2, gwas_cat)

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
        sig_rows = read_sig_out(sig[0])
        partners = gwas_partners(sig_rows)
        # Only built when a dump was asked for: it exists to measure the
        # merged-locus policy, not to feed the table.
        partners_strict = (gwas_partners(sig_rows, drop_merged=True)
                           if _TRAIT_DUMP is not None else None)
        gtex_causal_maf = causal_maf_by_trait(s2, gcat)

        # The fine-mapping settings and the causal floor are part of the GROUPING
        # key, not just reported: the same cell rerun at a different ld_ctrl, MAF
        # floor or causal floor is a separate experimental condition and must not
        # be pooled with the original.
        gt_filt = (fm_floor > 0 if params is not None
                   else detect_fm_filtered(s2, os.path.join(s3root, gcat), gcat, fm_floor))
        fmf = gw_filt if gw_filt is not None else gt_filt

        # How the causal set was drawn. Read from the params file only: these are
        # not recoverable from a path (n_central_traits and require_gtex_partner
        # appear in the run tag only when they differ from their legacy value, and
        # the run tag is not parsed here anyway), so a run without a params record
        # gets None for all of them and groups exactly as it did before.
        _p = params or {}
        draw = (_p.get("causal_sampling"), _p.get("sampling_gwas_n"),
                _p.get("sampling_power_plateau"), _p.get("n_central_traits"),
                _p.get("require_gtex_partner"))

        def row():
            # `draw` is every stage-2 knob that decides WHICH loci become causal
            # and is not already represented above. Without it, arms that differ
            # only in how the causal set was drawn pool into one row with a
            # doubled denominator -- silently, and looking like a well-powered
            # result. The publication set has two such pairs: paired vs unpaired
            # at identical floors (require_gtex_partner), and the two power arms
            # at identical floors (sampling_gwas_n).
            #
            # Every element is None for a run whose params file predates it, so
            # each is constant across any pre-existing root and the grouping there
            # is unchanged.
            key = (cat_letter, gtex_n, gwas_mult, gtex_mult, causal_maf,
                   fmf, ldc, source, draw, arm)
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

        # Colocalization is scored per GWAS trait, over EVERY defined GWAS trait --
        # iterating the declared list rather than `partners` keys keeps traits that
        # colocalized with nothing in the denominator instead of silently dropping
        # them and inflating power.
        for tr in gwas_traits:
            d = partners.get(tr, {})
            self_rcp, self_lcp = d.get(tr, (0.0, 0.0))
            others = {g: v for g, v in d.items() if g != tr}
            r["gwas_hits"].append((self_rcp, self_lcp, others))
            if _TRAIT_DUMP is not None:
                ds = partners_strict.get(tr, {})
                s_rcp, s_lcp = ds.get(tr, (0.0, 0.0))
                oth_s = {g: v for g, v in ds.items() if g != tr}
                _TRAIT_DUMP.append(dict(
                    root=os.path.basename(os.path.dirname(rep.rstrip("/"))),
                    rep=os.path.basename(rep), cat=cat_letter,
                    gwas_cat=gwas_cat, gtex_cat=gcat, gtex_n=gtex_n,
                    causal_min_maf=causal_maf, fm_min_maf=fm_floor,
                    gwas_mult=gwas_mult, gtex_mult=gtex_mult, trait=tr,
                    # Both sides: the GWAS panel and the GTEx panel filter
                    # independently, and a variant can clear the floor in one and
                    # not the other (the 500-person panel drops the most).
                    tested_gwas=int(tested_flag(gwas_causal_maf, tr, fm_floor)),
                    tested_gtex=int(tested_of(tr)),
                    # None when the sidecar is missing. False means the locus has
                    # NO GTEx trait at its causal variant, so self_rcp is 0 by
                    # construction and any n_other_* it carries is a false positive.
                    shared=gwas_shared.get(tr),
                    self_rcp=self_rcp, self_lcp=self_lcp,
                    n_other_50=sum(1 for v in others.values() if v[0] > 0.5),
                    n_other_90=sum(1 for v in others.values() if v[0] > 0.9),
                    self_rcp_strict=s_rcp, self_lcp_strict=s_lcp,
                    n_other_50_strict=sum(1 for v in oth_s.values() if v[0] > 0.5),
                    n_other_90_strict=sum(1 for v in oth_s.values() if v[0] > 0.9),
                ))

        # Separate pass over every GTEx gene for the FINE-MAPPING columns. Flank
        # genes contribute to cs_all (it describes every credible set the panel
        # produced) but not to the causal-variant columns, which are defined only
        # at the shared loci.
        for tr in gtex_traits:
            f = os.path.join(outdir, f"{tr}.dapg.out")
            if not os.path.exists(f):
                continue
            clus, cpip, csize = parse_dapg(f)
            cid = clus.get("snp" + tr[2:], -1)
            r["gt_cs_all"].extend(csize.values())
            if tr in gwas_traits and tested_of(tr):
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
    ap.add_argument("--dump-traits", metavar="PATH",
                    help="also write a long-format TSV with one row per GWAS "
                         "trait per GTEx panel, carrying the causal-variant "
                         "tested flags and the hit counts under both the default "
                         "merged-locus policy and the strict one")
    a = ap.parse_args()

    global _TRAIT_DUMP
    if a.dump_traits:
        _TRAIT_DUMP = []

    allrows = {}
    for root in a.roots:
        for k, v in collect(root).items():
            if k in allrows:
                for f in ("gwas_hits", "gw_cpip", "gt_cpip",
                          "gw_cs_all", "gw_cs_caus", "gt_cs_all", "gt_cs_caus"):
                    allrows[k][f].extend(v[f])
                allrows[k]["reps"] |= v["reps"]
                # Both counters, not just n_gwas: the earlier version omitted
                # n_causal_tested here, so a grouping key present in two roots
                # silently reported only the first root's tested count.
                for f in ("n_gwas", "n_causal_tested"):
                    allrows[k][f] += v[f]
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
            key=lambda kv: (kv[0][9] or "", kv[0][0], -kv[0][4], bool(kv[0][5]),
                            kv[0][6] if kv[0][6] else 0, -kv[0][1],
                            kv[0][2], kv[0][3],
                            tuple(str(x) for x in kv[0][8]))):
        cat, gtex_n, gwm, gtm, causal_maf, fmf, ldc, source, draw, arm = key
        sampling, samp_n, plateau, n_central, require_partner = draw
        hits = r["gwas_hits"]
        # (idx, thr) pairs in COLUMNS order: RCP at .5/.9 then LCP at .5/.9.
        # Each column family tests its own metric only -- LCP >= RCP always, so
        # the lcp counts are >= their rcp counterparts.
        counts = {k: power_and_fp(hits, *k)
                  for k in ((0, 0.5), (0, 0.9), (1, 0.5), (1, 0.9))}
        (r50, r90, l50, l90) = (counts[(0, 0.5)], counts[(0, 0.9)],
                                counts[(1, 0.5)], counts[(1, 0.9)])
        fm_filtered = "NA" if fmf is None else str(fmf)
        # Measured from Snakemake's record of the dap-g command line; the
        # --fm-r2 CLI value is only a fallback when that metadata is gone.
        fm_r2 = a.fm_r2 if ldc is None else ldc
        vals = [
            "NA" if arm is None else arm,
            cat, DEMOGRAPHY.get(cat, "?"), gtex_n, len(r["reps"]),
            causal_maf, r["n_gwas"], r["n_causal_tested"],
            len(r["gw_cpip"]), len(r["gt_cpip"]),
            fm_filtered, fm_r2, source, gwm, gtm,
            *("NA" if x is None else x
              for x in (sampling, samp_n, plateau, n_central, require_partner)),
            r50[0], r50[1], r90[0], r90[1], l50[0], l50[1], l90[0], l90[1],
            hit_med(hits, 0, 0.5), hit_med(hits, 0, 0.9),
            hit_med(hits, 1, 0.5), hit_med(hits, 1, 0.9),
            round(auc_point(hits, r["n_gwas"], 0, 0.5), 3),
            round(auc_point(hits, r["n_gwas"], 0, 0.9), 3),
            sum(1 for x in r["gw_cpip"] if x > 0.5),
            sum(1 for x in r["gt_cpip"] if x > 0.5),
            sum(1 for x in r["gw_cpip"] if x > 0.9),
            sum(1 for x in r["gt_cpip"] if x > 0.9),
            mean(r["gw_cs_all"]), mean(r["gw_cs_caus"]),
            mean(r["gt_cs_all"]), mean(r["gt_cs_caus"]),
            med(r["gw_cs_all"]), med(r["gw_cs_caus"]),
            med(r["gt_cs_all"]), med(r["gt_cs_caus"]),
        ]
        # A row that is one field short or long is not an error anywhere -- it is
        # a TSV whose columns have silently shifted from that row on.
        assert len(vals) == len(COLUMNS), (
            f"row has {len(vals)} values but COLUMNS declares {len(COLUMNS)}")
        print("\t".join(str(v) for v in vals), file=fh)
    if fh is not sys.stdout:
        fh.close()

    # Provenance check, on stderr so it never contaminates the TSV. A non-zero
    # over-cutoff count means colocalizing signals could not be attributed to a
    # GWAS trait and were excluded from every count -- i.e. the numbers are
    # understated and the sig.out naming needs re-checking.
    s = _SIG_STATS
    print(f"[sig.out] {s['rows']} signals read; {s['no_gwas']} with no GWAS locus "
          f"({s['no_gwas_over_half']} of those above 0.5)", file=sys.stderr)
    if s["no_gwas_over_half"]:
        print("[sig.out] WARNING: dropped signals that cleared 0.5 -- power and "
              "fp counts are understated.", file=sys.stderr)
    # Merged GWAS loci are credited to every constituent trait (see read_sig_out),
    # so a large over-cutoff count here means a meaningful share of the calls rest
    # on an attribution that fastEnloc itself left ambiguous.
    print(f"[sig.out] {s['merged']} signals spanned MERGED GWAS loci "
          f"({s['merged_over_half']} of those above 0.5), credited to every "
          f"constituent trait", file=sys.stderr)

    if _TRAIT_DUMP is not None:
        cols = list(_TRAIT_DUMP[0]) if _TRAIT_DUMP else []
        with open(a.dump_traits, "w") as dh:
            print("\t".join(cols), file=dh)
            for rec in _TRAIT_DUMP:
                print("\t".join(str(rec[c]) for c in cols), file=dh)
        print(f"[dump] {len(_TRAIT_DUMP)} trait records -> {a.dump_traits}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
