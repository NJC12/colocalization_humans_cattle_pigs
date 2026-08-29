"""Canonical filename + directory builders for the simulation Snakefile.

All paths derive from a small set of inputs so the four pipelines can share the
common rules. The convention:

    {basename}_{stagetag}_sd{seed}.{ext}

with stage-2 outputs additionally namespaced under
``gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{causal_min_maf}/`` (this layout is
hardcoded in ``create_gwas_files_and_phenotypes.py`` and we mirror it).

Note which MAF floor appears in that name. There are three independent ones, and
each governs exactly one thing:

    causal_min_maf  which variants may be selected as CAUSATIVE (stage 2)
    fm_min_maf      which variants are FINE-MAPPED (the SBAMS DAP-G reads, stage 3)
    min_maf         which variants enter the GWAS (``plink2 --maf``, stage 5)

They are set separately and none is derived from another. Only ``causal_min_maf``
determines stage-2 content -- the trait set and the phenotypes -- so it is the one
in the stage-2 directory name, and it is also what ``stage2_run_tag`` carries so
stages 3-5 separate too.
"""

import os

from . import human_demography


# ---------- causal-variant MAF floor ----------

# The PATH-SUPPRESSION SENTINEL, and nothing else. Every simulation output that
# existed before ``causal_min_maf`` became a knob was produced at 0.01, and those
# paths must stay reachable, so that one historical value emits no path segment.
# This is a fact about what is already on disk, NOT a relationship between
# causal_min_maf and min_maf -- the two are independent knobs and are never
# compared.
#
# It used to double as the default. It no longer does: the default is
# DEFAULT_CAUSAL_MIN_MAF below. Keeping the two decoupled is what lets the
# default change without renaming every directory ever written. Do not re-couple
# them.
LEGACY_CAUSAL_MIN_MAF = 0.01

# No floor unless a config asks for one. The rule is "causal_min_maf defaults to
# 0, but always applies when provided" -- in BOTH sampling schemes, which is the
# part that changed: the power path used to ignore the floor on the central pool
# on the grounds that the detection-power weight subsumed it.
#
# 0 is not "no filter": ``maf > 0`` well-formedness still applies, since a site
# monomorphic within a panel is not a variant there and an effect assigned to it
# would produce a phenotype of pure noise. See create_gwas_files_and_phenotypes.py
# :causal_eligible.
DEFAULT_CAUSAL_MIN_MAF = 0


def causal_min_maf(cfg):
    """MAF floor a variant must clear to be eligible as causative.

    Governs the central causal pool (in both sampling schemes), the central
    neutral recipient pool used under ``neutral_trait_vars`` (categories B/D),
    and the flanking GTEx-only loci. Independent of ``min_maf`` and
    ``fm_min_maf``; see the module docstring.
    """
    return cfg.get("causal_min_maf", DEFAULT_CAUSAL_MIN_MAF)


def causal_maf_token(value):
    """Render a causal MAF floor the ONE way every path component must render it.

    Three producers have to agree on this string: ``stage2_inner`` (the output
    directory), ``causal_maf_segment`` (the stages-2..5 namespace), and the
    f-string in ``create_gwas_files_and_phenotypes.py`` that names the directory
    the script actually writes into. They agreed by accident while every value
    was 0.01 or 0.001, which formats identically however it arrives.

    Zero does not. The config carries YAML int ``0`` -> ``"0"``, but the stage-2
    script parses ``--causal_min_maf`` as a float and gets ``"0.0"``. Snakemake
    would then declare ``maf_0/`` as the rule's output while the script wrote
    ``maf_0.0/``, and the rule would fail with its outputs missing after stage 2
    had already done all the work.

    Integral values render without the trailing ``.0``; everything else renders
    as-is. Every value used before this existed is unchanged.
    """
    f = float(value)
    return str(int(f)) if f == int(f) else str(f)


def causal_maf_segment(cfg):
    """Path segment marking a non-historical causal MAF floor, else ``""``.

    Appended to ``stage2_run_tag`` so that stages 3, 4 and 5 -- which are
    otherwise namespaced by the stage-1 filename alone -- separate between runs
    that drew their causal variants from different pools. Empty at
    ``LEGACY_CAUSAL_MIN_MAF``, which keeps every pre-existing path unchanged.

    Note that the suppressed value is 0.01, not the default. A run at the
    default (0) therefore DOES emit ``.cmaf_0``. That is deliberate: 0 and 0.01
    draw from different pools and must not share a directory, and pinning the
    suppression to the historical value is the only way existing outputs stay
    adoptable.
    """
    if float(causal_min_maf(cfg)) == LEGACY_CAUSAL_MIN_MAF:
        return ""
    return f".cmaf_{causal_maf_token(causal_min_maf(cfg))}"


# ---------- causal-variant sampling scheme ----------

# The historical scheme: a uniform draw from the pool that clears causal_min_maf.
# Like LEGACY_CAUSAL_MIN_MAF this is a fact about what is already on disk -- it emits
# no path segment so every existing path stays reachable.
LEGACY_CAUSAL_SAMPLING = "uniform"

#: ``sampling_power_plateau`` value that leaves the weights alone. Power is
#: already bounded by 1, so this is a no-op and every pre-existing arm sits here.
NO_POWER_PLATEAU = 1.0


def power_plateau_token(value):
    """Render a plateau the one way every path component must render it.

    Drops the decimal point, matching how the arm names already write floors
    (``cmaf001`` for 0.001, ``fm01`` for 0.01): 0.5 -> ``05``, 0.25 -> ``025``.
    Keeping the dot out matters -- ``stage2_run_tag`` is dot-joined and
    ``summarize_coloc.read_params`` prefix-globs the result.
    """
    return f"{float(value):g}".replace(".", "")


def causal_sampling(cfg):
    """How the central causal loci are drawn: ``"uniform"`` or ``"power"``.

    Under ``"power"`` each central ``selco != 0`` variant is weighted by its
    probability of being detected in a GWAS of ``sampling_gwas_n`` individuals and
    the draw has inclusion probability proportional to that power (see
    ``helpers/causal_power.py``).

    ``causal_min_maf`` gates the central pool under BOTH schemes. It did not
    always: power sampling used to drop the floor on the grounds that the weight
    subsumed it, which meant the same ``causal_min_maf`` value described two
    different pools depending on the scheme. The two knobs are still not
    redundant -- the floor is a hard cut on frequency, the weight is a soft
    preference over (effect size, frequency) jointly -- they now just compose.
    """
    return cfg.get("causal_sampling", LEGACY_CAUSAL_SAMPLING)


def causal_sampling_segment(cfg):
    """Path segment marking a non-historical sampling scheme, else ``""``.

    Carries only the headline knob, ``sampling_gwas_n``: two power runs differing
    solely in ``sampling_sig_p`` land on the same path and are caught by the
    Snakefile's stage-2 provenance guard instead, which fails loudly rather than
    silently reusing the wrong phenotypes. Path segments are for the parameters a
    reader needs to tell runs apart at a glance; the params file is for the rest.
    """
    if causal_sampling(cfg) == LEGACY_CAUSAL_SAMPLING:
        return ""
    seg = f".psamp_{int(cfg['sampling_gwas_n'])}"
    # The plateau is the exception to "headline knob only". It changes WHICH
    # variants are drawn, not just how they are described, so a saturated and an
    # unsaturated run at the same sampling_gwas_n are different simulations and
    # must not share a directory. sampling_sig_p can share one because the
    # provenance guard catches it; this cannot, because the whole point of the
    # arm is to sit beside its unsaturated twin.
    plateau = float(cfg.get("sampling_power_plateau", NO_POWER_PLATEAU))
    if plateau != NO_POWER_PLATEAU:
        seg += f"_sat{power_plateau_token(plateau)}"
    return seg


# ---------- number of trait loci ----------

#: The central-trait count every run through round 3 used. Suppressed in the path
#: for the same reason LEGACY_CAUSAL_MIN_MAF is: it is a fact about what is
#: already on disk, so those paths stay reachable.
LEGACY_N_CENTRAL_TRAITS = 50


def trait_count_segment(cfg):
    """Path segment marking a non-historical central-trait count, else ``""``.

    ``n_central_traits`` determines stage-2 CONTENT -- it is how many causal loci
    the arm has -- but unlike causal_min_maf it appears in no filename. The
    stage-2 outputs are named ``gwas_G_gtex_T_maf_M/`` and the vars tables
    ``{panel}_vars_gwas_G_gtex_T_maf_M.tsv``, none of which carries it. So two
    runs differing only in the trait count would write the same filenames, and in
    the FLAT archive layout (``<arm>/<replicate>/<file>``) they would land on top
    of each other with nothing to tell them apart.

    Carries the central count only, matching causal_sampling_segment's rule of
    "the headline knob, and let the provenance guard catch the rest":
    ``n_flank_gtex_traits`` is in params_record.STAGE2_KEYS, so a run that
    disagrees on it is refused rather than silently reusing the wrong phenotypes.
    """
    n = cfg.get("n_central_traits")
    if n is None or int(n) == LEGACY_N_CENTRAL_TRAITS:
        return ""
    return f".ntr_{int(n)}"


# ---------- neutral-variant thinning (categories O and P) ----------

#: ``neutral_keep_fraction`` value that keeps every neutral site, i.e. the whole
#: neutral overlay as simulated. Every run before categories O and P sits here and
#: emits no path segment, so no pre-existing path moves.
NO_NEUTRAL_THINNING = 1.0


def neutral_keep_fraction(cfg):
    """Fraction of ``selco == 0`` sites the run keeps. 1.0 = the whole overlay.

    Neutral mutations are a Poisson process laid down on branches the forward
    simulation already fixed -- ``human_simulation_o2.py``'s 8.4e-9 overlay for the
    human arms, ``create_gwas_files_and_phenotypes.add_neutral``'s for the cattle
    ones. Dropping each neutral site independently with probability ``1 - k`` is
    therefore EXACTLY the same distribution as running that overlay at
    ``k * 8.4e-9`` (Poisson thinning), which is what lets categories O and P reuse
    A's and E's stage-1 tree sequences instead of re-simulating at a lower rate.

    Only the neutral class is touched, so the causal pools, the power weights, the
    drawn positions and the trait partner tables are identical to the un-thinned
    arm at the same seed -- PROVIDED the causal pool is the selected class. That is
    true for A/E and every arm that reads effect sizes off the tree sequence, and
    it is why O and P are defined against A and E. It is NOT true under
    ``synthetic_dfe_effects``, where ``causal_eligible`` inverts its predicate to
    ``selco == 0`` and the neutral class IS the causal pool -- thinning it would
    change the draw, not just the background. The Snakefile rejects that
    combination rather than relying on this note.
    """
    return float(cfg.get("neutral_keep_fraction", NO_NEUTRAL_THINNING))


def neutral_thin_segment(cfg):
    """Path segment marking a thinned neutral class, else ``""``.

    Thinning changes stage-2 CONTENT -- it is the variant set every panel, every
    SBAMS and every GWAS is built from -- but like ``n_central_traits`` it appears
    in no filename: the vars tables are named ``{panel}_vars_gwas_G_gtex_T_maf_M``
    whether or not the neutral class was thinned. Without a segment an O run and an
    A run could be pointed at the same stage 2 by ``stage2_search_dirs``, and in the
    flat archive layout their files would be indistinguishable.

    The token follows ``power_plateau_token``'s rule -- drop the decimal point, the
    way the arm names already write floors -- so 0.375 renders ``nkeep_0375``.
    """
    k = neutral_keep_fraction(cfg)
    if k == NO_NEUTRAL_THINNING:
        return ""
    return f".nkeep_{power_plateau_token(k)}"


# ---------- GWAS/GTEx pairing of the causal loci ----------

#: ``require_gtex_partner`` value that reproduces the historical DERIVED rule:
#: intersect the central causal pool with the GTEx panel iff the draw is uniform
#: AND the effect sizes are the variants' own. Every run before this key existed
#: sits here and emits no path segment, so no pre-existing path moves.
LEGACY_REQUIRE_GTEX_PARTNER = "auto"


def require_gtex_partner_auto(cfg):
    """The rule this knob replaces: intersect iff uniform and non-synthetic.

    Kept as a function rather than inlined because it is the definition of
    ``auto``, and ``require_gtex_partner_segment`` has to compare against it to
    know whether an explicit value is actually saying anything new.
    """
    return (causal_sampling(cfg) != "power"
            and not cfg.get("synthetic_dfe_effects", False))


def require_gtex_partner(cfg):
    """Must a central causal locus segregate in the GTEx panel to be drawn?

    ``True`` intersects the pool, so every GWAS causal locus has a GTEx partner
    by construction and the two panels' causal sets are the same positions.
    ``False`` does not, so the shared GTEx set is whichever drawn loci that panel
    happens to carry and the rest of its central slots are topped up uniformly --
    which means a GWAS locus can have no partner at all and cannot colocalize.
    ``"auto"`` is ``require_gtex_partner_auto``.

    This used to be derived rather than set, which made two experimental
    conditions unreachable: a uniform draw WITHOUT the intersection, and a
    synthetic-DFE or power draw WITH it. Both are arms of the publication run --
    the pairing is the thing being varied, so it has to be a knob.
    """
    return cfg.get("require_gtex_partner", LEGACY_REQUIRE_GTEX_PARTNER)


def require_gtex_partner_segment(cfg):
    """Path segment marking a pairing rule the derived one would not have given.

    Pairing changes stage-2 CONTENT -- which loci become causal, and whether a
    GWAS locus has a GTEx partner -- but like ``n_central_traits`` it appears in
    no filename. Two arms differing only in this key would write the same
    filenames, and in the flat archive layout (``<arm>/<replicate>/<file>``) they
    would land on top of each other with nothing to tell them apart.

    Emitted only when the resolved value DIFFERS from ``require_gtex_partner_auto``,
    so a run that merely states the rule it would have got anyway is
    byte-identical to one that says nothing -- which is what lets the launcher
    pass this key for all 120 publication runs without moving any existing path.
    """
    v = require_gtex_partner(cfg)
    if v == LEGACY_REQUIRE_GTEX_PARTNER or bool(v) == require_gtex_partner_auto(cfg):
        return ""
    return ".gtexreq_1" if v else ".gtexreq_0"


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

    Matches the actual filename produced by ``human_simulation_o2.py``. The
    total individual count and SLiM scaling are controlled by the
    ``--n_samples`` and ``--Q`` CLI flags rather than by selecting a separate
    script, and the sampled population by ``--population``.

    ``hts_{seed}.ts`` for EUR under the default demographic model -- every human
    tree sequence written before either key was a parameter, so A/B/C/D stay
    adoptable -- ``hts_{pop}_{seed}.ts`` for another population of that model
    (category J: ``hts_AFR_101.ts``), and ``hts_{tag}_{pop}_{seed}.ts`` for
    another model entirely (category M: ``hts_finwang_FIN_131.ts``). Both the
    model and the population are in the name because stage-1 adoption matches on
    filename plus seed -- and for human names the seed half never fires, so the
    filename is the whole guard. See helpers/human_demography.py.

    ``.get`` rather than ``[...]``: helpers/tests feeds raw config YAML that has
    not been through the Snakefile's ``setdefault``, and no config written
    before category J carries either key at all.
    """
    model = cfg.get("demographic_model") or human_demography.DEFAULT_MODEL
    population = (cfg.get("population")
                  or human_demography.default_population(model))
    return human_demography.stage1_ts_name(population, cfg["stage1_seed"], model)


def stage1_human_neutral_ts(cfg):
    """Top-level human_neutral tree sequence, from human_neutral_simulation.py.

    ``nts_`` rather than ``hts_``: the two arms share a seed band convention and a
    workdir layout, and stage-1 adoption via ``stage1_search_dirs`` matches on
    filename plus embedded seed. A shared prefix would let an A tree sequence
    satisfy an H run's input -- silently, and with a genome shaped by selection.
    """
    return f"nts_{cfg['stage1_seed']}.ts"


def stage1_cattle_neutral_ts(cfg):
    """Top-level cattle_neutral tree sequence, from cattle_neutral_simulation.py.

    ``cnts_`` continues the ``hts_``/``nts_`` convention and is distinct from both
    for the same reason they are distinct from each other: stage-1 adoption via
    ``stage1_search_dirs`` matches on filename plus embedded seed, and a shared
    prefix would let one arm's tree sequence satisfy another's input silently.
    Nothing else in the cattle family uses this shape -- E/F/G carry the long
    descriptive ``farm_*`` names because their filenames have to encode the SLiM
    parameters that produced them, and this one has no forward phase to describe.
    """
    return f"cnts_{cfg['stage1_seed']}.ts"


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


#: The one cattle deep history every publication run through the 5-replicate set
#: resumed from. Suppressed in the path for exactly the reason
#: LEGACY_CAUSAL_MIN_MAF is: it is a fact about what is already on disk, so those
#: paths stay reachable.
LEGACY_CATTLE_BASELINE_SEED = 20250303


def cattle_baseline_segment(cfg):
    """Path segment naming the cattle deep history, else ``""``.

    ``cattle_baseline_seed`` determines stage-1 CONTENT -- it is the entire
    ancestry the run resumes from, 33,130 of 33,154 generations -- but for the
    ``cattle_sel`` pipeline it appears in no filename at all. Two F runs at the
    same ``stage1_seed`` off different deep histories therefore write the same
    ``.full.ts``, the same ``.m4_marks.tsv``, the same ``stage2_run_tag``, and
    after publishing the same files in the same directory, with nothing to tell
    them apart. ``cattle_baseline_seed`` is not in ``params_record.STAGE2_KEYS``
    either, so the provenance guard cannot separate them; the sidecar records the
    stage-1 tree sequence by NAME, which is precisely what this segment fixes.

    (``stage1_cattle_baseline_from_midpoint_full`` has never had the problem: it
    carries ``cb_{cb}`` unconditionally, so E and L were always distinguishable.)

    Emitted only when the seed is not ``LEGACY_CATTLE_BASELINE_SEED``, so every F
    and G output written before four deep histories existed is byte-identical
    under this change and no finished run moves. Same construction as
    ``causal_maf_segment``, ``causal_sampling_segment``, ``trait_count_segment``,
    ``neutral_thin_segment`` and ``require_gtex_partner_segment``.

    Placed AFTER ``_sd{seed}`` so ``search_dirs._extract_seed``'s first match is
    still the stage-1 seed, and so the collision glob becomes
    ``..._sd*_cb{cb}...`` -- scoped to one deep history rather than reporting a
    sibling history's file as a seed conflict.
    """
    cb = cfg.get("cattle_baseline_seed")
    if cb is None or str(cb) == str(LEGACY_CATTLE_BASELINE_SEED):
        return ""
    return f"_cb{cb}"


def stage1_cattle_sel_full(cfg):
    """Output of farm_selection_from_ep8.slim (num_muts_selected > 0) — cattle+selection full TS."""
    base = cfg["basename"]
    mult = cfg["selection_multiplier"]
    gen = cfg["selected_generations"]
    muts = cfg["num_muts_selected"]
    s1 = cfg["stage1_seed"]
    return (f"{base}_mult{mult}_gen{gen}_muts{muts}_sd{s1}"
            f"{cattle_baseline_segment(cfg)}.full.ts")


def stage1_cattle_sel_marks(cfg):
    base = cfg["basename"]
    mult = cfg["selection_multiplier"]
    gen = cfg["selected_generations"]
    muts = cfg["num_muts_selected"]
    s1 = cfg["stage1_seed"]
    return (f"{base}_mult{mult}_gen{gen}_muts{muts}_sd{s1}"
            f"{cattle_baseline_segment(cfg)}.m4_marks.tsv")


def stage1_full_filename(cfg):
    """The single .full.ts that stage 2 will consume, by pipeline."""
    p = cfg["pipeline"]
    if p == "human":
        return stage1_human_ts(cfg)
    if p == "human_neutral":
        return stage1_human_neutral_ts(cfg)
    if p == "cattle_neutral":
        return stage1_cattle_neutral_ts(cfg)
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
    (so changing stage1_seed gets a fresh stage-2 dir), plus the causal MAF
    floor when it is not the historical 0.01, plus the causal sampling scheme
    when it is not the historical uniform draw, plus the central-trait count
    when it is not the historical 50, plus the neutral keep fraction when the
    neutral class was thinned, plus the GWAS/GTEx pairing rule when it is not the
    one the derived rule would have chosen.

    Segments are appended in that order and each is empty at its legacy value, so
    every tag ever written stays a strict prefix of the tag a later, more
    parameterised run produces -- which is what ``summarize_coloc.read_params``
    relies on when it prefix-globs ``{tag}.*.params.txt``.

    This is the sole namespace for stage2_dir, stage3_dir, stage4_dir and
    stage5_dir, and for the ``previous_workdirs`` adoption in the Snakefile, so
    appending causal_maf_segment() here separates the whole downstream tree at
    once -- including ``geno.sbams.gz``, the file every DAP-G job reads.
    """
    base = os.path.splitext(os.path.splitext(stage1_full_filename(cfg))[0])[0]
    return (base + causal_maf_segment(cfg) + causal_sampling_segment(cfg)
            + trait_count_segment(cfg) + neutral_thin_segment(cfg)
            + require_gtex_partner_segment(cfg))


def stage2_dir(cfg):
    """Top-level stage-2 dir for this run."""
    return os.path.join(workdir(cfg), "stage2", stage2_run_tag(cfg))


def stage2_inner(cfg):
    """The directory that ``create_gwas_files_and_phenotypes.py`` actually
    writes to (one level deeper, matching its hardcoded layout).

    The ``maf_`` component is causal_min_maf, NOT min_maf: it is the only one of
    the three floors that determines what lives in this directory (the trait set
    and the phenotypes). min_maf governs ``plink2 --maf`` in stage 5 and
    fm_min_maf the SBAMS handed to DAP-G, and neither touches stage-2 content.
    Both are 0.01 in every config through round 3, so this name is unchanged for
    every run that already exists.

    Naming it this way also means the all-or-nothing stage-2 adoption in the
    Snakefile -- which matches on this directory's basename -- cannot hand a run
    phenotypes that were built under a different causal floor.
    """
    gw = cfg["gwas_scaling"]
    gt = cfg["gtex_scaling"]
    return os.path.join(
        stage2_dir(cfg),
        f"gwas_{gw}_gtex_{gt}_maf_{causal_maf_token(causal_min_maf(cfg))}")


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


# ---------- per-simulation parameters record ----------
#
# Two copies of one file, written at Snakemake onstart by helpers/params_record.
# Neither is a declared rule output -- same reasoning as the ``.fmmaf`` sidecar in
# rules/common.smk: runs that predate the file must stay adoptable, and a
# metadata file should not add a job to every existing DAG.

def stage4_params_file(cfg):
    """Beside the stage-4 outputs, sharing their prefix.

    For any ``{prefix}.{gtex_cat}.enloc.{kind}.out``, strip
    ``.{gtex_cat}.enloc.{kind}.out`` and append ``.params.txt`` to get this file.
    Uniqueness comes from the directory: stage4_dir embeds stage2_run_tag (the
    stage-1 filename, hence stage1_seed, plus the causal MAF segment).
    ``basename`` alone is not unique -- every round-3 config says
    ``basename: human``.
    """
    return os.path.join(
        stage4_dir(cfg),
        f"{cfg['basename']}.{_stage4_prefix_tag_segment(cfg)}params.txt",
    )


def workdir_params_file(cfg):
    """Per-invocation copy at the workdir root, with identical content.

    Exists as soon as Snakemake starts, so the record survives a DAG that dies
    before stage 4 ever runs.
    """
    tag = output_tag(cfg)
    stem = stage2_run_tag(cfg) + (f".{tag}" if tag else "")
    return os.path.join(workdir(cfg), "params", f"{stem}.params.txt")


def stage2_params_file(cfg):
    """Sidecar written by create_gwas_files_and_phenotypes.py itself.

    Records the stage-2 parameters actually used plus the resulting pool sizes.
    Stage-2 directories are symlink-adopted across output roots, so this is the
    only way an adopting run can know what the adopted phenotypes were -- and it
    is the oracle for the Snakefile's stage-2 provenance guard.
    """
    return stage2_file(cfg, "stage2_params.txt")
