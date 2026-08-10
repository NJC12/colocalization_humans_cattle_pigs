# Stage 1 -- human_neutral (category H). A single Python step, coalescent only.
#
# Invokes human_neutral_simulation.py, which writes nts_{seed}.ts. Stage 2 reads
# that .ts directly through the --human_ts_file branch, exactly as it does the A
# arm's hts_{seed}.ts: msprime times are already in generations, which is what
# that branch's times_already_unscaled=True asserts.
#
# There is no SLiM here and no Q_scaling to undo, so this rule has none of
# stage1_human.smk's resource tiering. The A arm's tiers exist because a forward
# simulation's runtime and working set scale with L and with the population size;
# a coalescent over 2 Mb and 9k diploids is minutes and a couple of GB.


def _human_neutral_n_samples(cfg):
    """n_samples = gwas_size + largest GTEx size needed downstream.

    Same rule as the A arm (rules/stage1_human.smk): when gtex_size == -1 the
    multi-GTEx fan-out tops out at 1000.
    """
    gwas_size = int(cfg["gwas_size"])
    gtex_size = int(cfg.get("gtex_size", -1))
    largest_gtex = 1000 if gtex_size == -1 else gtex_size
    return gwas_size + largest_gtex


rule stage1_human_neutral:
    output:
        ts = paths.stage1_dir(config) + "/" + paths.stage1_human_neutral_ts(config),
    input:
        existing = lambda wc: search_dirs.find_or_empty(
            paths.stage1_human_neutral_ts(config),
            config.get("stage1_search_dirs", []),
            expected_seed=config["stage1_seed"],
        ),
    params:
        seed       = config["stage1_seed"],
        L          = config["L"],
        recomb     = config["recombination_rate"],
        # Required by the Snakefile for this pipeline, so no .get() default: the
        # total-vs-neutral-component choice is the one number that decides whether
        # this arm's variant density matches the arm it is compared against.
        mut_rate   = config["mutation_rate"],
        n_samples  = _human_neutral_n_samples(config),
        script_dir = SIM_REPO_DIR,
        script     = "human_neutral_simulation.py",
        python     = config["python_binary"],
        tmp_dir    = os.path.join(paths.workdir(config), "human_neutral_tmp"),
        produced_ts = os.path.join(paths.workdir(config), "human_neutral_tmp",
                                   f"nts_{config['stage1_seed']}.ts"),
    log:
        os.path.join(paths.workdir(config), "logs", "stage1_human_neutral.log"),
    resources:
        slurm_partition = "short",
        runtime         = 60,
        mem_mb          = 16000,
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
        mkdir -p "{params.tmp_dir}"
        "{params.python}" {params.script} \
            --seed {params.seed} \
            --length {params.L} \
            --n_samples {params.n_samples} \
            --recomb_rate {params.recomb} \
            --mut_rate {params.mut_rate} \
            --tmp_dir "{params.tmp_dir}" \
            > {log} 2>&1
        cp -u "{params.produced_ts}" "{output.ts}"
        """
