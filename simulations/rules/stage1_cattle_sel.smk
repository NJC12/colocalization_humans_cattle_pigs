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
        paths.stage1_cattle_baseline_epoch_8(config),
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


rule stage1_cattle_sel:
    output:
        full  = paths.stage1_dir(config) + "/" + paths.stage1_cattle_sel_full(config),
        marks = paths.stage1_dir(config) + "/" + paths.stage1_cattle_sel_marks(config),
    input: unpack(_stage1_cattle_sel_inputs)
    params:
        seed        = config["stage1_seed"],
        Q           = config["Q_scaling"],
        L           = config["L"],
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
        # Loads epoch_8.ts and runs ~24 scaled ticks (epochs 8-12, populations
        # capped at 1,000 individuals decreasing to 90) on a 10 Mb genome.
        # Hours of compute, not days; 30 days/150 GB was a copy-paste from
        # run_farm_selection.sh which runs the full epoch 1-12 chain.
        slurm_partition = "short",
        runtime = 8 * 60,
        mem_mb  = 32000,
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
            echo "stage1_cattle_sel: no epoch_8 checkpoint found in cattle_baseline_search_dirs" >&2
            echo "Add the dir containing farm_selection_Q_*.L_*.seed_*.epoch_8.ts to cattle_baseline_search_dirs in your config." >&2
            exit 1
        fi
        {params.slim} -m -t -l \
            -s {params.seed} \
            -d Q_scaling={params.Q} \
            -d genome_length={params.L} \
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
