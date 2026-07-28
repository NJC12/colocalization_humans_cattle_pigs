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
        recomb     = config["recombination_rate"],
        script_dir = SIM_REPO_DIR,
        python     = config["python_binary"],
    log: os.path.join(paths.workdir(config), "logs", "stage1_cattle_create_pop.log")
    resources:
        slurm_partition = "priority,long",
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
        "{params.python}" farm_create_orig_pop_e2.py \
            --seed {params.seed} --Q {params.Q} --length {params.L} \
            --recomb_rate {params.recomb} \
            --out_dir "$(dirname {output.ts})" \
            --out_name "$(basename {output.ts})" \
            > {log} 2>&1
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
        recomb      = config["recombination_rate"],
        slim        = config.get("slim_binary", "/home/njc12/bin/slim/build/slim"),
        # farm_burn_in_continue.slim, not farm_burn_in_e2.slim. The latter
        # terminates at a hard-coded tick 25000 and checkpoints every 10/Q ticks,
        # neither of which survives contact with the Q=1 deep history: it would
        # stop 5000 ticks short of burn_in_cycle=30000 (so the declared output
        # never appears) and write 3000 checkpoint files on the way. The continue
        # script takes both as arguments and reads an orig-pop at tick 1 just as
        # happily as it resumes from a checkpoint.
        slim_script = os.path.join(SIM_REPO_DIR, "farm_burn_in_continue.slim"),
        end_tick    = config.get("burn_in_cycle", 25000),
        # Checkpoints feed helpers/sfs_equilibrium.py; ~25-30 time points is
        # plenty. Falls back to the old 10/Q cadence when unset.
        ckpt_every  = config.get("burn_in_checkpoint_every",
                                 max(1, int(10 / config["Q_scaling"]))),
        out_dir     = paths.stage1_dir(config),
        out_base    = paths.stage1_cattle_baseline_burnin(config).rsplit(".cycle_", 1)[0],
    log: os.path.join(paths.workdir(config), "logs", "stage1_cattle_burn_in.log")
    resources:
        slurm_partition = "priority,long",
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
            -d recomb_rate={params.recomb} \
            -d end_tick={params.end_tick} \
            -d checkpoint_every={params.ckpt_every} \
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
        recomb      = config["recombination_rate"],
        tick        = config.get("burn_in_cycle", 25000),
        slim        = config.get("slim_binary", "/home/njc12/bin/slim/build/slim"),
        slim_script = os.path.join(SIM_REPO_DIR, "farm_selection.slim"),
        out_dir     = paths.stage1_dir(config),
        out_base    = paths.stage1_cattle_baseline_full(config).replace(".full.ts", ""),
    log: os.path.join(paths.workdir(config), "logs", "stage1_cattle_selection.log")
    resources:
        slurm_partition = "priority,long",
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
            -d recomb_rate={params.recomb} \
            -d tick={params.tick} \
            -d file_in='"{input.upstream}"' \
            -d dir_out='"{params.out_dir}/"' \
            -d file_out='"{params.out_base}"' \
            "{params.slim_script}" \
            > {log} 2>&1
        """
