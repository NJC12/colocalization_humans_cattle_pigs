# Stage 1 -- cattle baseline (3-job chain).
#
# 1) farm_create_orig_pop_e2.py: msprime coalescent + selection coefficients
#    -> farm_orig_pop.Q_{Q}.L_{L}.seed_{s1}.ts
# 2) farm_burn_in_e2.slim: burn-in to mutation-selection-drift equilibrium
#    -> farm_burn_in.Q_{Q}.L_{L}.seed_{s1}.cycle_{tick}.ts (checkpoints)
# 3) farm_selection.slim: 12-epoch demographic history under selection
#    -> farm_selection_Q_{Q}.L_{L}.seed_{s1}.full.ts
#
# Each rule uses a conditional input function: when a search-dir hit exists
# the rule's only input is that hit (and the shell symlinks it); otherwise it
# pulls in the upstream canonical path (and the shell runs the producer).
# This way a single search-dir hit short-circuits the chain without forcing
# every upstream stage to also run.


def _stage1_create_pop_inputs(wc):
    hit = search_dirs.find_in_search_dirs(
        paths.stage1_cattle_baseline_orig(config),
        config.get("cattle_baseline_search_dirs", []),
        expected_seed=config["stage1_seed"],
    )
    return {"existing": hit or []}


rule stage1_cattle_create_pop:
    output:
        ts = paths.stage1_dir(config) + "/" + paths.stage1_cattle_baseline_orig(config),
    input: unpack(_stage1_create_pop_inputs)
    params:
        seed       = config["stage1_seed"],
        Q          = config["Q_scaling"],
        L          = config["L"],
        script_dir = SIM_REPO_DIR,
    log: os.path.join(paths.workdir(config), "logs", "stage1_cattle_create_pop.log")
    resources:
        slurm_partition = "short",
        runtime = 4 * 60,
        mem_mb  = 50000,
        cpus_per_task = 4,
    conda: "../envs/coloc_sims.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.ts})" "$(dirname {log})"
        if [ -n "{input.existing}" ]; then
            ln -sf "{input.existing}" "{output.ts}"
            exit 0
        fi
        cd "{params.script_dir}"
        python farm_create_orig_pop_e2.py \
            --seed {params.seed} --Q {params.Q} --length {params.L} \
            > {log} 2>&1
        SRC=/home/njc12/slim_simulations/farm_slim_outputs/$(basename {output.ts})
        cp -u "$SRC" "{output.ts}"
        """


def _stage1_burn_in_inputs(wc):
    hit = search_dirs.find_in_search_dirs(
        paths.stage1_cattle_baseline_burnin(config),
        config.get("cattle_baseline_search_dirs", []),
        expected_seed=config["stage1_seed"],
    )
    if hit:
        return {"existing": hit, "upstream": []}
    return {
        "existing": [],
        "upstream": paths.stage1_dir(config) + "/" + paths.stage1_cattle_baseline_orig(config),
    }


rule stage1_cattle_burn_in:
    output:
        ts = paths.stage1_dir(config) + "/" + paths.stage1_cattle_baseline_burnin(config),
    input: unpack(_stage1_burn_in_inputs)
    params:
        seed        = config["stage1_seed"],
        Q           = config["Q_scaling"],
        L           = config["L"],
        slim        = config.get("slim_binary", "/home/njc12/bin/slim/build/slim"),
        slim_script = os.path.join(SIM_REPO_DIR, "farm_burn_in_e2.slim"),
        out_dir     = paths.stage1_dir(config),
        out_base    = paths.stage1_cattle_baseline_burnin(config).rsplit(".cycle_", 1)[0],
    log: os.path.join(paths.workdir(config), "logs", "stage1_cattle_burn_in.log")
    resources:
        slurm_partition = "priority",
        runtime = 30 * 24 * 60,
        mem_mb  = 150000,
        cpus_per_task = 4,
    conda: "../envs/coloc_sims.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.ts})" "$(dirname {log})"
        if [ -n "{input.existing}" ]; then
            ln -sf "{input.existing}" "{output.ts}"
            exit 0
        fi
        {params.slim} -m -t -l \
            -s {params.seed} \
            -d Q_scaling={params.Q} \
            -d length={params.L} \
            -d file_in='"{input.upstream}"' \
            -d dir_out='"{params.out_dir}/"' \
            -d file_out='"{params.out_base}"' \
            "{params.slim_script}" \
            > {log} 2>&1
        """


def _stage1_selection_inputs(wc):
    hit = search_dirs.find_in_search_dirs(
        paths.stage1_cattle_baseline_full(config),
        config.get("cattle_baseline_search_dirs", []),
        expected_seed=config["stage1_seed"],
    )
    if hit:
        return {"existing": hit, "upstream": []}
    return {
        "existing": [],
        "upstream": paths.stage1_dir(config) + "/" + paths.stage1_cattle_baseline_burnin(config),
    }


rule stage1_cattle_selection:
    output:
        full = paths.stage1_dir(config) + "/" + paths.stage1_cattle_baseline_full(config),
    input: unpack(_stage1_selection_inputs)
    params:
        seed        = config["stage1_seed"],
        Q           = config["Q_scaling"],
        L           = config["L"],
        tick        = config.get("burn_in_cycle", 25000),
        slim        = config.get("slim_binary", "/home/njc12/bin/slim/build/slim"),
        slim_script = os.path.join(SIM_REPO_DIR, "farm_selection.slim"),
        out_dir     = paths.stage1_dir(config),
        out_base    = paths.stage1_cattle_baseline_full(config).replace(".full.ts", ""),
    log: os.path.join(paths.workdir(config), "logs", "stage1_cattle_selection.log")
    resources:
        slurm_partition = "priority",
        runtime = 30 * 24 * 60,
        mem_mb  = 150000,
        cpus_per_task = 4,
    conda: "../envs/coloc_sims.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.full})" "$(dirname {log})"
        if [ -n "{input.existing}" ]; then
            ln -sf "{input.existing}" "{output.full}"
            exit 0
        fi
        {params.slim} -m -t -l \
            -s {params.seed} \
            -d Q_scaling={params.Q} \
            -d length={params.L} \
            -d tick={params.tick} \
            -d file_in='"{input.upstream}"' \
            -d dir_out='"{params.out_dir}/"' \
            -d file_out='"{params.out_base}"' \
            "{params.slim_script}" \
            > {log} 2>&1
        """
