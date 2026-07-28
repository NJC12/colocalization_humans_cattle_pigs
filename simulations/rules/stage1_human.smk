# Stage 1 -- human (single Python step that calls SLiM via stdpopsim).
#
# Invokes human_simulation_o2.py, which writes hts_{seed}.ts. Stage 2 reads
# that .ts directly. Q (SLiM scaling) and n_samples (total EUR individuals)
# are passed as CLI args; n_samples is computed from gwas_size + the largest
# requested GTEx size for the run.


def _human_stage1_n_samples(cfg):
    """n_samples = gwas_size + largest GTEx size needed downstream.

    When gtex_size == -1 (multi-GTEx fan-out: 1000/500/250) the largest is 1000.
    For single-size runs (C/D with gtex_size=50000) use that value directly.
    """
    gwas_size = int(cfg["gwas_size"])
    gtex_size = int(cfg.get("gtex_size", -1))
    largest_gtex = 1000 if gtex_size == -1 else gtex_size
    return gwas_size + largest_gtex


# The three functions below key on BOTH sample count and L. They used to key on
# sample count alone, which meant every A/B run inherited figures calibrated at
# L=10 Mb: 7 days and 100 GB. The measured round-3 A1 run at L=2 Mb took 1:05
# and 232 MB -- 0.01% of the runtime request and 0.2% of the memory. Since the
# SLiM working set and runtime both scale with sequence length, small-L configs
# get their own tier rather than a blanket cut, so the 10 Mb path is unchanged
# if it is ever revived.
_SMALL_L = 4_000_000


def _human_stage1_small(cfg):
    return _human_stage1_n_samples(cfg) <= 50000 and cfg["L"] <= _SMALL_L


def _human_stage1_mem_mb(cfg):
    """Scale memory request to the actual SLiM working-set size.

    Round-3 A1 measured 232 MB at L=2 Mb / n=9k; 8 GB is 35x that. A/B at
    L=10 Mb / n=9k need ~100 GB (was 32 GB at L=1 Mb -- the working set scales
    with L). C/D at L=2.5 Mb / n=100k keep 200 GB.
    """
    if _human_stage1_small(cfg):
        return 8000
    return 200000 if _human_stage1_n_samples(cfg) > 50000 else 100000


def _human_stage1_partition(cfg):
    """Small-L A/B fits `short` comfortably at 2h, which also removes the need
    for controller_2Mb*.sbatch to override this rule's partition. A/B at
    L=10 Mb overshoot short's 12h cap, so route to priority (30-day). C/D keep
    the long partition's 30-day allotment."""
    if _human_stage1_small(cfg):
        return "short"
    return "long" if _human_stage1_n_samples(cfg) > 50000 else "priority"


def _human_stage1_runtime(cfg):
    """Small-L A/B: 2h, against a measured 1:05 at L=2 Mb. A/B at L=10 Mb:
    7 days as a safe ceiling on priority. C/D keep the prior 30-day
    reservation on long."""
    if _human_stage1_small(cfg):
        return 2 * 60
    return 30 * 24 * 60 if _human_stage1_n_samples(cfg) > 50000 else 7 * 24 * 60


rule stage1_human:
    output:
        ts = paths.stage1_dir(config) + "/" + paths.stage1_human_ts(config),
    input:
        existing = lambda wc: search_dirs.find_or_empty(
            paths.stage1_human_ts(config),
            config.get("stage1_search_dirs", []),
            expected_seed=config["stage1_seed"],
        ),
    params:
        seed       = config["stage1_seed"],
        gwas_h2    = config.get("gwas_h2", 0.01),
        gtex_h2    = config.get("gtex_h2", 0.4),
        L          = config["L"],
        Q          = config["Q_scaling"],
        recomb     = config["recombination_rate"],
        n_samples  = _human_stage1_n_samples(config),
        script_dir = SIM_REPO_DIR,
        script     = "human_simulation_o2.py",
        python     = config["python_binary"],
        slim       = config["slim_binary"],
        # Per-run scratch dir keeps round-2 outputs isolated from old runs.
        tmp_dir    = os.path.join(paths.workdir(config), "human_tmp"),
        produced_ts = os.path.join(paths.workdir(config), "human_tmp", f"hts_{config['stage1_seed']}.ts"),
    log:
        os.path.join(paths.workdir(config), "logs", "stage1_human.log"),
    resources:
        slurm_partition = _human_stage1_partition(config),
        runtime         = _human_stage1_runtime(config),
        mem_mb          = _human_stage1_mem_mb(config),
        cpus_per_task   = 4,
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
        export PATH="$(dirname {params.slim}):$PATH"
        mkdir -p "{params.tmp_dir}"
        "{params.python}" {params.script} \
            --seed {params.seed} \
            --gwas_h2 {params.gwas_h2} \
            --gtex_h2 {params.gtex_h2} \
            --length {params.L} \
            --Q {params.Q} \
            --recomb_rate {params.recomb} \
            --n_samples {params.n_samples} \
            --tmp_dir "{params.tmp_dir}" \
            > {log} 2>&1
        cp -u "{params.produced_ts}" "{output.ts}"
        """
