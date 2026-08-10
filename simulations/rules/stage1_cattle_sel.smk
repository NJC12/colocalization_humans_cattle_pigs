# Stage 1 -- cattle + selection (used for both bottlenecked + not_bottlenecked
# pipelines; the YAML supplies continue_bottlenecking=1 or 0).
#
# Reads an epoch-8 checkpoint from a baseline cattle simulation and runs
# farm_selection_from_ep8.slim with num_muts_selected > 0, so positive
# selection is applied to a chosen set of mutations at sample_tick.
# Produces:
#   {basename}_mult{mult}_gen{gen}_muts{muts}_sd{s1}.full.ts
#   {basename}_mult{mult}_gen{gen}_muts{muts}_sd{s1}.m4_marks.tsv
#
# Conditional input function: if both products exist in a search dir, the
# rule's only inputs are those hits (shell symlinks them). Otherwise we need
# the epoch_8 checkpoint, also looked up via search dirs (since the
# cattle_baseline rules are not included in this pipeline).


def _stage1_cattle_sel_inputs(wc):
    full_hit = search_dirs.find_in_search_dirs(
        paths.stage1_cattle_sel_full(config),
        config.get("stage1_search_dirs", []),
        expected_seed=config["stage1_seed"],
    )
    marks_hit = search_dirs.find_in_search_dirs(
        paths.stage1_cattle_sel_marks(config),
        config.get("stage1_search_dirs", []),
        expected_seed=config["stage1_seed"],
    )
    if full_hit and marks_hit:
        return {"existing_full": full_hit, "existing_marks": marks_hit, "epoch_8": []}
    epoch_8_hit = search_dirs.find_in_search_dirs(
        paths.stage1_cattle_baseline_handoff(config),
        config.get("cattle_baseline_search_dirs", []),
        expected_seed=config["cattle_baseline_seed"],
    )
    if not epoch_8_hit:
        # Defer the error to job runtime (Snakemake at DAG-build time would
        # crash with a confusing MissingInputException). Empty input here
        # produces a clearer error inside the shell.
        epoch_8_hit = []
    return {
        "existing_full":  [],
        "existing_marks": [],
        "epoch_8": epoch_8_hit,
    }


# Attempt-scaled because both numbers are guesses -- see the note in the rule's
# resources block. Named functions rather than lambdas, matching _dapg_runtime /
# _dapg_mem_mb in common.smk.
def _stage1_cattle_sel_runtime(wc, attempt):
    return 2 * 60 * attempt        # 2 h -> 4 h -> 6 h; `short` caps at 12 h


def _stage1_cattle_sel_mem_mb(wc, attempt):
    return 8000 * attempt          # 8 -> 16 -> 24 GB


rule stage1_cattle_sel:
    output:
        full  = paths.stage1_dir(config) + "/" + paths.stage1_cattle_sel_full(config),
        marks = paths.stage1_dir(config) + "/" + paths.stage1_cattle_sel_marks(config),
    input: unpack(_stage1_cattle_sel_inputs)
    params:
        seed        = config["stage1_seed"],
        Q           = config["Q_scaling"],
        L           = config["L"],
        recomb      = config["recombination_rate"],
        mult        = config["selection_multiplier"],
        gen         = config["selected_generations"],
        muts        = config["num_muts_selected"],
        bottleneck  = config["continue_bottlenecking"],
        slim        = config.get("slim_binary", "/home/njc12/bin/slim/build/slim"),
        slim_script = os.path.join(SIM_REPO_DIR, "farm_selection_from_ep8.slim"),
        out_dir     = paths.stage1_dir(config),
        out_base    = paths.stage1_cattle_sel_full(config).replace(".full.ts", ""),
    log: os.path.join(paths.workdir(config), "logs", "stage1_cattle_sel.log")
    resources:
        # Loads the handoff checkpoint and runs epochs 8-12 (populations 1,000
        # decreasing to 90, census 1/Q times that) -- the same SLiM script as
        # stage1_cattle_baseline_from_midpoint, plus positive selection.
        #
        # Kept at 2h/8 GB rather than that rule's measured 1h/4 GB: categories
        # F and G have NOT been run under round-3 code, so there is no direct
        # measurement here. E1 took 4:29 / 600 MB at L=2 Mb, and adding selected
        # mutations should not change that much, but G in particular behaved
        # very differently from E in round 2 (see _dapg_mem_mb_base). Tighten
        # once F/G have actually run.
        #
        # ATTEMPT-SCALED because both numbers are guesses. They were static, so a
        # G stage-1 that wanted 9 GB would fail identically on all three attempts
        # (--restart-times is 2) and take the replicate down with it. G is the
        # unmeasured case: continue_bottlenecking=0 leaves a ~100k population, an
        # order of magnitude more individuals than E/F carry through epochs 8-12.
        # Scaling to 16/24 GB and 4/6 h makes a wrong first guess recoverable
        # instead of fatal. `short` caps at 12 h, so attempt 3 still fits.
        slurm_partition = "short",
        runtime = _stage1_cattle_sel_runtime,
        mem_mb  = _stage1_cattle_sel_mem_mb,
        cpus_per_task = 4,
    conda: "../envs/coloc_sims.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.full})" "$(dirname {log})"
        if [ -n "{input.existing_full}" ] && [ -n "{input.existing_marks}" ]; then
            ln -sf "{input.existing_full}"  "{output.full}"
            ln -sf "{input.existing_marks}" "{output.marks}"
            exit 0
        fi
        if [ -z "{input.epoch_8}" ]; then
            echo "stage1_cattle_sel: no end-of-epoch-7 handoff checkpoint found in cattle_baseline_search_dirs" >&2
            echo "Add the dir containing farm_selection_Q_*.L_*.seed_*.ep7.ts to cattle_baseline_search_dirs in your config." >&2
            exit 1
        fi
        {params.slim} -m -t -l \
            -s {params.seed} \
            -d Q_scaling={params.Q} \
            -d genome_length={params.L} \
            -d recomb_rate={params.recomb} \
            -d selection_multiplier={params.mult} \
            -d selected_generations={params.gen} \
            -d num_muts_selected={params.muts} \
            -d continue_bottlenecking={params.bottleneck} \
            -d file_in='"{input.epoch_8}"' \
            -d dir_out='"{params.out_dir}/"' \
            -d file_out='"{params.out_base}"' \
            "{params.slim_script}" \
            > {log} 2>&1
        """
