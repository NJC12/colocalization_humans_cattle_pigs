# Stage 1 -- cattle baseline continued from epoch 8.
#
# Reads the epoch-8 checkpoint from a predecessor cattle_baseline run and runs
# farm_selection_from_ep8.slim through epochs 8-12 with no artificial selection.
# Produces a baseline-equivalent .full.ts in a fraction of the wall time of the
# full 3-stage cattle_baseline chain.
#
# Output:
#   farm_selection_from_ep8.Q_{Q}.L_{L}.cb_{cattle_baseline_seed}.seed_{stage1_seed}.full.ts
#
# Conditional input function: if the full output already exists in a search
# dir (seed-matched), the rule's only input is that hit (shell symlinks it).
# Otherwise we need the epoch_8 checkpoint, looked up via
# cattle_baseline_search_dirs (since the cattle_baseline rules are not
# included in this pipeline).


def _stage1_baseline_from_midpoint_inputs(wc):
    full_hit = search_dirs.find_in_search_dirs(
        paths.stage1_cattle_baseline_from_midpoint_full(config),
        config.get("stage1_search_dirs", []),
        expected_seed=config["stage1_seed"],
    )
    if full_hit:
        return {"existing_full": full_hit, "epoch_8": []}
    epoch_8_hit = search_dirs.find_in_search_dirs(
        paths.stage1_cattle_baseline_epoch_8(config),
        config.get("cattle_baseline_search_dirs", []),
        expected_seed=config["cattle_baseline_seed"],
    )
    if not epoch_8_hit:
        # Defer the error to job runtime; Snakemake at DAG-build time would
        # crash with a less actionable MissingInputException.
        epoch_8_hit = []
    return {"existing_full": [], "epoch_8": epoch_8_hit}


rule stage1_cattle_baseline_from_midpoint:
    output:
        full = paths.stage1_dir(config) + "/" + paths.stage1_cattle_baseline_from_midpoint_full(config),
    input: unpack(_stage1_baseline_from_midpoint_inputs)
    params:
        seed        = config["stage1_seed"],
        Q           = config["Q_scaling"],
        L           = config["L"],
        slim        = config.get("slim_binary", "/home/njc12/bin/slim/build/slim"),
        slim_script = os.path.join(SIM_REPO_DIR, "farm_selection_from_ep8.slim"),
        out_dir     = paths.stage1_dir(config),
        out_base    = paths.stage1_cattle_baseline_from_midpoint_full(config).replace(".full.ts", ""),
    log: os.path.join(paths.workdir(config), "logs", "stage1_cattle_baseline_from_midpoint.log")
    resources:
        # Epochs 8-12 on a 10 Mb genome, populations 1000 -> 90. Hours, not
        # days. Match cattle_sel's resource profile.
        slurm_partition = "short",
        runtime = 8 * 60,
        mem_mb  = 32000,
        cpus_per_task = 4,
    conda: "../envs/coloc_sims.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.full})" "$(dirname {log})"
        if [ -n "{input.existing_full}" ]; then
            ln -sf "{input.existing_full}" "{output.full}"
            exit 0
        fi
        if [ -z "{input.epoch_8}" ]; then
            echo "stage1_cattle_baseline_from_midpoint: no epoch_8 checkpoint found in cattle_baseline_search_dirs" >&2
            echo "Add the dir containing farm_selection_Q_*.L_*.seed_*.epoch_8.ts to cattle_baseline_search_dirs in your config." >&2
            exit 1
        fi
        {params.slim} -m -t -l \
            -s {params.seed} \
            -d Q_scaling={params.Q} \
            -d genome_length={params.L} \
            -d num_muts_selected=0 \
            -d file_in='"{input.epoch_8}"' \
            -d dir_out='"{params.out_dir}/"' \
            -d file_out='"{params.out_base}"' \
            "{params.slim_script}" \
            > {log} 2>&1
        """
