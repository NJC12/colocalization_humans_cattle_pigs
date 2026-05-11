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
        n_samples  = _human_stage1_n_samples(config),
        script_dir = SIM_REPO_DIR,
        script     = "human_simulation_o2.py",
        # The script writes to a hardcoded scratch dir.
        produced_ts = f"/n/scratch/users/n/njc12/sims/tmp/hts_{config['stage1_seed']}.ts",
    log:
        os.path.join(paths.workdir(config), "logs", "stage1_human.log"),
    resources:
        slurm_partition = "priority",
        runtime = 30 * 24 * 60,   # 30 days, in minutes
        mem_mb  = 200000,
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
        python {params.script} \
            --seed {params.seed} \
            --gwas_h2 {params.gwas_h2} \
            --gtex_h2 {params.gtex_h2} \
            --length {params.L} \
            --Q {params.Q} \
            --n_samples {params.n_samples} \
            > {log} 2>&1
        cp -u "{params.produced_ts}" "{output.ts}"
        """
