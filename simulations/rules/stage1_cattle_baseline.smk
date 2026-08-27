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


# ---------------------------------------------------------------------------
# Resource tiers, by Q_scaling.
#
# These three rules had ONE tier each, sized for Q=0.01, and asked for
# `priority,long` throughout. Two problems, both found the first time the chain
# was actually driven through Snakemake rather than by hand:
#
#   1. `stage1_cattle_create_pop` asks for 4 hours. O2 REFUSES a sub-12-hour job
#      on `long` -- "please submit jobs that are less than 12 hours long to the
#      short (or priority) partition" -- and rejects the whole comma-list as an
#      invalid partition name. So that rule could never submit at all. It went
#      unnoticed because every cattle deep history to date was launched as a
#      standalone sbatch (jobs `cattle_origpop_q001`, `cattle_burnin_q1`), not
#      through this rule.
#   2. The burn-in and selection rules ask for 30 days and 150 GB. That is right
#      at Q=0.01, where the population is 17,000/0.01 = 1.7M individuals -- the
#      `cattle_burnin_hedge` run has spent ~32 days to reach 6.6% of its target.
#      It is wrong by three orders of magnitude at Q=1, which is what the round-3
#      deep history actually runs: MEASURED 38m48s for the 30,000-tick burn-in
#      (sacct 48879426, 32 GB) and ~40 min for epochs 1-11. A 30-day request on
#      `long` for a 40-minute job queues for hours behind real long jobs.
#
# Cost scales as roughly 1/Q^2: the burn-in runs `end_tick` ticks over
# `N_ancestral / Q` individuals. So the tier is chosen on Q, the way
# stage1_human.smk's _human_stage1_partition already does.
_CATTLE_FULL_SCALE_Q = 1.0


def _cattle_deep_is_rescaled(cfg):
    """True when Q < 1, i.e. the population is inflated by 1/Q and this is the
    multi-week regime rather than the sub-hour one."""
    return float(cfg["Q_scaling"]) < _CATTLE_FULL_SCALE_Q


def _cattle_deep_partition(cfg):
    # `short` accepts <= 12 h; `long` REFUSES anything shorter. Never both.
    return "priority,long" if _cattle_deep_is_rescaled(cfg) else "short"


def _cattle_deep_runtime(cfg, hours_at_full_scale):
    if _cattle_deep_is_rescaled(cfg):
        return 30 * 24 * 60
    return int(hours_at_full_scale * 60)


def _cattle_deep_mem_mb(cfg, mb_at_full_scale):
    return 150000 if _cattle_deep_is_rescaled(cfg) else mb_at_full_scale


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
        # 4 h either way (measured 13-16 min at Q=0.01, sacct cattle_origpop_q001),
        # but the PARTITION has to change: `long` refuses a 4-hour job.
        slurm_partition = _cattle_deep_partition(config),
        runtime = 4 * 60,
        mem_mb  = _cattle_deep_mem_mb(config, 32000),
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
        slurm_partition = _cattle_deep_partition(config),
        runtime = _cattle_deep_runtime(config, 8),
        mem_mb  = _cattle_deep_mem_mb(config, 32000),
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
        slurm_partition = _cattle_deep_partition(config),
        runtime = _cattle_deep_runtime(config, 8),
        mem_mb  = _cattle_deep_mem_mb(config, 32000),
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
