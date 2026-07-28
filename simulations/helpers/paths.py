"""Canonical filename + directory builders for the simulation Snakefile.

All paths derive from a small set of inputs so the four pipelines can share the
common rules. The convention:

    {basename}_{stagetag}_sd{seed}.{ext}

with stage-2 outputs additionally namespaced under
``gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}/`` (this layout is
hardcoded in ``create_gwas_files_and_phenotypes.py`` and we mirror it).
"""

import os


# ---------- output tag (for parallel result sets, e.g. r2=0.25 vs r2=0.75) ----------

def output_tag(cfg):
    """Optional tag inserted into stage-3 / stage-4 paths.

    When empty (default), stage-3 outputs live in ``outputs/`` / ``logs/`` and
    stage-4 prefix is ``{basename}.{gtex_cat}`` — i.e. paths are unchanged from
    the original layout. When set (e.g. ``"r2_0_25"``), stage-3 outputs move to
    ``outputs_<tag>/`` / ``logs_<tag>/`` with a ``.stage3.<tag>.done`` sentinel,
    and stage-4 prefix becomes ``{basename}.<tag>.{gtex_cat}``. Stages 1, 2,
    and 5 are tag-agnostic so they are shared across r2 values.
    """
    return cfg.get("output_tag", "") or ""


def _tagged_dir(name, tag):
    return name if not tag else f"{name}_{tag}"


# ---------- per-pipeline stage-1 filename builders ----------

def stage1_human_ts(cfg):
    """Top-level human tree sequence (one file feeds stage 2 directly).

    Matches the actual filename produced by ``human_simulation_o2.py``:
    ``hts_{seed}.ts``. The total individual count and SLiM scaling are
    now controlled by the ``--n_samples`` and ``--Q`` CLI flags rather
    than by selecting a separate script.
    """
    return f"hts_{cfg['stage1_seed']}.ts"


def stage1_cattle_baseline_orig(cfg):
    """Output of farm_create_orig_pop_e2.py."""
    Q = cfg["Q_scaling"]
    L = cfg["L"]
    s1 = cfg["stage1_seed"]
    return f"farm_orig_pop.Q_{Q}.L_{L}.seed_{s1}.ts"


def stage1_cattle_baseline_burnin(cfg):
    """Output of farm_burn_in_e2.slim (final tick checkpoint we actually use)."""
    Q = cfg["Q_scaling"]
    L = cfg["L"]
    s1 = cfg["stage1_seed"]
    cycle = cfg.get("burn_in_cycle", 25000)
    return f"farm_burn_in.Q_{Q}.L_{L}.seed_{s1}.cycle_{cycle}.ts"


def stage1_cattle_baseline_full(cfg):
    """Output of farm_selection.slim — the cattle-baseline full TS.

    Note the inconsistent separator: actual on-disk files use underscore after
    ``farm_selection`` (``farm_selection_Q_0.01...``) while ``farm_burn_in`` and
    ``farm_orig_pop`` use a period (``farm_burn_in.Q_0.01...``).
    """
    Q = cfg["Q_scaling"]
    L = cfg["L"]
    s1 = cfg["stage1_seed"]
    return f"farm_selection_Q_{Q}.L_{L}.seed_{s1}.full.ts"


def stage1_cattle_baseline_epoch_8(cfg):
    """END-of-epoch-8 checkpoint from farm_selection.slim.

    Retained for reference/back-compat only. Do NOT feed this to
    farm_selection_from_ep8.slim: that script begins by setting the epoch-8
    population size and running ep_lengths[8] worth of ticks, so starting it
    from a state that has already been through epoch 8 runs epoch 8 twice. Every
    cattle arm through round 2 therefore had 12 generations at Ne=1,000 instead
    of 6 -- 0.003 of spurious coalescent time, 5.7% of all drift between epoch 8
    and sampling. Use stage1_cattle_baseline_handoff() instead.
    """
    Q = cfg["Q_scaling"]
    L = cfg["L"]
    s1 = cfg["cattle_baseline_seed"]
    return f"farm_selection_Q_{Q}.L_{L}.seed_{s1}.epoch_8.ts"


def stage1_cattle_baseline_handoff(cfg):
    """END-of-epoch-7 checkpoint — the correct entry point for the from-ep8 runs.

    farm_selection.slim writes this at `startE8`, immediately before it applies
    the epoch-8 population size, so it is the state epochs 8-12 should start
    from. It is also where the Q=1 deep history naturally stops, which is why
    helpers/rescale_checkpoint.py writes its output under this name.
    """
    Q = cfg["Q_scaling"]
    L = cfg["L"]
    s1 = cfg["cattle_baseline_seed"]
    return f"farm_selection_Q_{Q}.L_{L}.seed_{s1}.ep7.ts"


def stage1_cattle_baseline_from_midpoint_full(cfg):
    """Output of farm_selection_from_ep8.slim — cattle baseline continued from epoch 8.

    Embeds both cattle_baseline_seed (predecessor whose epoch_8.ts is loaded)
    and stage1_seed (this continuation run) so adoption is seed-strict on
    both axes and the file cannot be confused with a full-trajectory baseline
    output.
    """
    Q = cfg["Q_scaling"]
    L = cfg["L"]
    cb = cfg["cattle_baseline_seed"]
    s1 = cfg["stage1_seed"]
    return f"farm_selection_from_ep8.Q_{Q}.L_{L}.cb_{cb}.seed_{s1}.full.ts"


def stage1_cattle_sel_full(cfg):
    """Output of farm_selection_from_ep8.slim (num_muts_selected > 0) — cattle+selection full TS."""
    base = cfg["basename"]
    mult = cfg["selection_multiplier"]
    gen = cfg["selected_generations"]
    muts = cfg["num_muts_selected"]
    s1 = cfg["stage1_seed"]
    return f"{base}_mult{mult}_gen{gen}_muts{muts}_sd{s1}.full.ts"


def stage1_cattle_sel_marks(cfg):
    base = cfg["basename"]
    mult = cfg["selection_multiplier"]
    gen = cfg["selected_generations"]
    muts = cfg["num_muts_selected"]
    s1 = cfg["stage1_seed"]
    return f"{base}_mult{mult}_gen{gen}_muts{muts}_sd{s1}.m4_marks.tsv"


def stage1_full_filename(cfg):
    """The single .full.ts that stage 2 will consume, by pipeline."""
    p = cfg["pipeline"]
    if p == "human":
        return stage1_human_ts(cfg)
    if p == "cattle_baseline":
        return stage1_cattle_baseline_full(cfg)
    if p == "cattle_baseline_from_midpoint":
        return stage1_cattle_baseline_from_midpoint_full(cfg)
    if p == "cattle_sel":
        return stage1_cattle_sel_full(cfg)
    raise ValueError(f"Unknown pipeline: {p!r}")


def stage1_marks_filename(cfg):
    """The .m4_marks.tsv that stage 2 needs (cattle_sel only)."""
    if cfg["pipeline"] != "cattle_sel":
        return None
    return stage1_cattle_sel_marks(cfg)


# ---------- workdir layout ----------

def workdir(cfg):
    return os.path.expanduser(cfg["workdir"])


def stage1_dir(cfg):
    return os.path.join(workdir(cfg), "stage1")


def stage2_run_tag(cfg):
    """Per-run tag that pins stage-2 outputs to the upstream stage-1 file
    (so changing stage1_seed gets a fresh stage-2 dir)."""
    return os.path.splitext(os.path.splitext(stage1_full_filename(cfg))[0])[0]


def stage2_dir(cfg):
    """Top-level stage-2 dir for this run."""
    return os.path.join(workdir(cfg), "stage2", stage2_run_tag(cfg))


def stage2_inner(cfg):
    """The directory that ``create_gwas_files_and_phenotypes.py`` actually
    writes to (one level deeper, matching its hardcoded layout)."""
    gw = cfg["gwas_scaling"]
    gt = cfg["gtex_scaling"]
    maf = cfg["min_maf"]
    return os.path.join(stage2_dir(cfg), f"gwas_{gw}_gtex_{gt}_maf_{maf}")


def stage2_marker(cfg):
    return os.path.join(stage2_inner(cfg), ".stage2.done")


def stage2_file(cfg, name):
    """Path to a specific file inside the stage-2 output dir."""
    return os.path.join(stage2_inner(cfg), name)


# ---------- stage 2 file names ----------

def species_prefix(cfg):
    return "c" if cfg["species"] == "cattle" else "h"


def stage2_categories(cfg):
    """The categories produced by stage 2.

    Default is two: ``{p}gwas`` and ``{p}gtex``. When ``gtex_size: -1`` is set
    in the config, the stage-2 script additionally produces two sub-sampled
    GTEx datasets (500 and 250 individuals) for the power-analysis sweep -- we
    expose them here so stages 3-5 fan out over all of them.

    A category listed in ``skip_categories`` (config list) is filtered out of
    the downstream fan-out. Stage-2 outputs for that category that already
    exist on disk are left intact; downstream stages just stop targeting it.

    NOTE (2026-05-26): ``{p}gtex_smallest`` is permanently disabled from the
    downstream fan-out. The trait-simulation Python script still writes the
    smallest pheno/geno files (we want those preserved), but no DAP-G /
    fastEnloc / plink-GLM is run on them. To re-enable, un-comment the line
    below and restart any affected controllers.
    """
    p = species_prefix(cfg)
    cats = [f"{p}gwas", f"{p}gtex"]
    if cfg.get("gtex_size") == -1:
        cats += [f"{p}gtex_smaller"]
        # cats += [f"{p}gtex_smaller", f"{p}gtex_smallest"]  # disabled 2026-05-26
    skip = set(cfg.get("skip_categories", []) or [])
    return [c for c in cats if c not in skip]


def stage2_gtex_categories(cfg):
    """Just the GTEx categories (one or three, depending on gtex_size)."""
    return [c for c in stage2_categories(cfg) if not c.endswith("gwas")]


def stage2_gwas_category(cfg):
    """The single GWAS category for this species."""
    return next(c for c in stage2_categories(cfg) if c.endswith("gwas"))


def stage2_pheno(cfg, cat):
    scaling = cfg["gwas_scaling"] if cat.endswith("gwas") else cfg["gtex_scaling"]
    return stage2_file(cfg, f"{cat}_scaling_{scaling}_pheno.sbams")


def stage2_geno(cfg, cat):
    scaling = cfg["gwas_scaling"] if cat.endswith("gwas") else cfg["gtex_scaling"]
    return stage2_file(cfg, f"{cat}_scaling_{scaling}_geno.sbams")


# stage2_pca() removed: stage 2 no longer writes {cat}.pca and stage 3 no longer
# consumes it. See scripts/run_plink_glm.sh for why the PCA covariates went away.


def stage2_vcf(cfg, cat):
    """Cattle writes .vcf, human writes .vcf.gz."""
    if cfg["species"] == "cattle":
        return stage2_file(cfg, f"{cat}.vcf")
    return stage2_file(cfg, f"{cat}.vcf.gz")


def stage2_pheno_for_plink(cfg, cat):
    scaling = cfg["gwas_scaling"] if cat.endswith("gwas") else cfg["gtex_scaling"]
    return stage2_file(cfg, f"{cat}_traits.scaling_{scaling}.pheno")


def stage2_dap_annot_vcf(cfg, cat):
    """Annotation VCF produced from the stage-2 VCF for fastEnloc.

    Matches the existing convention from create_gwas_files_and_phenotypes.py's
    helper text: ``cgwas.dap.vcf.gz`` produced by the awk recipe.
    """
    return stage2_file(cfg, f"{cat}.dap.vcf.gz")


# ---------- stage 3 (DAP-G) ----------
#
# Stage 3 paths intentionally do NOT embed scaling in filenames. The {cat}
# directory (cgwas/cgtex/hgwas/hgtex) is enough to disambiguate, and dropping
# scaling from stage-3 names lets Snakemake `output:` patterns use just one
# wildcard ({cat} or {cat}+{trait}). Two flavors of helper:
#   - concrete (takes cat) for input lambdas
#   - _tmpl (returns a string with {cat}/{trait} wildcards) for `output:`

def stage3_dir(cfg, cat):
    return os.path.join(workdir(cfg), "stage3", stage2_run_tag(cfg), cat)


def stage3_dir_tmpl(cfg):
    return os.path.join(workdir(cfg), "stage3", stage2_run_tag(cfg), "{cat}")


def stage3_geno_gz(cfg, cat):
    return os.path.join(stage3_dir(cfg, cat), "geno.sbams.gz")


def stage3_geno_gz_tmpl(cfg):
    return os.path.join(stage3_dir_tmpl(cfg), "geno.sbams.gz")


def stage3_geno_tbi(cfg, cat):
    return stage3_geno_gz(cfg, cat) + ".tbi"


def stage3_geno_tbi_tmpl(cfg):
    return stage3_geno_gz_tmpl(cfg) + ".tbi"


def stage3_manifest(cfg, cat):
    return os.path.join(stage3_dir(cfg, cat), "manifest.txt")


def stage3_manifest_tmpl(cfg):
    return os.path.join(stage3_dir_tmpl(cfg), "manifest.txt")


def stage3_outputs_subdir_name(cfg):
    """Subdir name (``outputs`` or ``outputs_<tag>``) where DAP-G writes loci.

    Passed to ``run_fastenloc.sh`` via ``--outputs-subdir`` so dap2enloc
    scans only the matching r2-tagged outputs.
    """
    return _tagged_dir("outputs", output_tag(cfg))


def _stage3_logs_subdir_name(cfg):
    return _tagged_dir("logs", output_tag(cfg))


def _stage3_done_basename(cfg):
    tag = output_tag(cfg)
    return ".stage3.done" if not tag else f".stage3.{tag}.done"


def stage3_locus_out(cfg, cat, trait):
    return os.path.join(stage3_dir(cfg, cat), stage3_outputs_subdir_name(cfg), f"{trait}.dapg.out")


def stage3_locus_out_tmpl(cfg):
    return os.path.join(stage3_dir_tmpl(cfg), stage3_outputs_subdir_name(cfg), "{trait}.dapg.out")


def stage3_locus_log_tmpl(cfg):
    return os.path.join(stage3_dir_tmpl(cfg), _stage3_logs_subdir_name(cfg), "{trait}.dapg.out")


def stage3_done(cfg, cat):
    return os.path.join(stage3_dir(cfg, cat), _stage3_done_basename(cfg))


def stage3_done_tmpl(cfg):
    return os.path.join(stage3_dir_tmpl(cfg), _stage3_done_basename(cfg))


# ---------- stage 4 (fastEnloc) ----------

def stage4_dir(cfg):
    return os.path.join(workdir(cfg), "stage4", stage2_run_tag(cfg))


def _stage4_prefix_tag_segment(cfg):
    tag = output_tag(cfg)
    return "" if not tag else f"{tag}."


def stage4_prefix(cfg, gtex_cat):
    """Output prefix for a single fastEnloc run (one per GTEx variant).

    The gtex_cat is embedded in the filename so the multi-size power-analysis
    sweep gets one file set per gtex variant (e.g. ``{basename}.hgtex.*``,
    ``{basename}.hgtex_smaller.*``, ``{basename}.hgtex_smallest.*``). When
    ``output_tag`` is non-empty it is inserted between ``basename`` and
    ``gtex_cat`` so r2-tagged result sets coexist (e.g.
    ``{basename}.r2_0_25.hgtex.*``).
    """
    return os.path.join(stage4_dir(cfg), f"{cfg['basename']}.{_stage4_prefix_tag_segment(cfg)}{gtex_cat}")


def stage4_output_kinds(cfg):
    """Which .enloc.*.out files fastEnloc actually writes, given the prior mode.

    `enrich` and `mi` are produced inside fastEnloc's enrich_est(), which is
    skipped entirely when colocalization priors are supplied directly
    (--coloc_default_prior sets p1/p2/p12 and the binary prints "Applying
    user-specified colocalization priors, skipping enrichment analysis").
    Listing them as rule outputs in that mode makes every stage-4 job fail with
    MissingOutputException even though fastEnloc exited 0.
    """
    if cfg.get("fastenloc_prior", "coloc_default") == "estimated":
        return ("enrich", "gene", "mi", "sig", "snp")
    return ("gene", "sig", "snp")


def stage4_outputs(cfg, gtex_cat):
    p = stage4_prefix(cfg, gtex_cat)
    return [f"{p}.enloc.{kind}.out" for kind in stage4_output_kinds(cfg)]


def stage4_prefix_tmpl(cfg):
    """Output-pattern form with a ``{gtex_cat}`` wildcard (for use in Snakemake `output:`)."""
    return os.path.join(stage4_dir(cfg), f"{cfg['basename']}.{_stage4_prefix_tag_segment(cfg)}{{gtex_cat}}")


def stage4_outputs_tmpl(cfg):
    p = stage4_prefix_tmpl(cfg)
    return [f"{p}.enloc.{kind}.out" for kind in stage4_output_kinds(cfg)]


# ---------- stage 5 (plink GLM) ----------

def stage5_dir(cfg):
    return os.path.join(workdir(cfg), "stage5", stage2_run_tag(cfg))


def stage5_prefix(cfg, cat):
    return os.path.join(stage5_dir(cfg), f"{cat}_glm")


def stage5_done(cfg, cat):
    return os.path.join(stage5_dir(cfg), f".{cat}.stage5.done")


def stage5_done_tmpl(cfg):
    return os.path.join(stage5_dir(cfg), ".{cat}.stage5.done")


def stage5_prefix_tmpl(cfg):
    return os.path.join(stage5_dir(cfg), "{cat}_glm")
