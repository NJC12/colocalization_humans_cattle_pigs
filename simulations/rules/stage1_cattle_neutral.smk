# Stage 1 -- cattle_neutral (category I). A single Python step, coalescent only.
#
# Invokes cattle_neutral_simulation.py, which writes cnts_{seed}.ts. Stage 2 reads
# it through the --cattle_ts_file branch, exactly as it does E/F/G's .full.ts:
# branch lengths are in Q-rescaled ticks, so get_vars_df's `time * Q_scaling` and
# add_neutral's `8.4e-9 / (1/Q_scaling)` are both correct without a
# times_already_unscaled exemption.
#
# NO SHARED CHECKPOINT. E/F/G declare cattle_baseline_search_dirs and load the
# .ep7.ts that farm_selection.slim produced once for all of them, because 29,800
# ticks of burn-in plus epochs 2-7 is unaffordable per replicate in SLiM. This
# arm simulates all twelve epochs itself -- backwards, in seconds -- so it has no
# baseline dependency and its configs carry neither key. The pre-flight guard in
# submit_2Mb_r3_cmaf_replicates.sh that checks for that checkpoint is scoped to
# ^(E|F|G)$ and must stay that way.
#
# There is no SLiM here, so this rule has none of stage1_cattle_sel.smk's
# resource tiering; the shape is stage1_human_neutral.smk's instead.


def _cattle_neutral_n_samples(cfg):
    """n_samples = gwas_size + largest GTEx size needed downstream.

    Same rule as the human arms. At Q_scaling=0.01 this comes to 9,000, which is
    also the terminal population size (epoch 12 = 90 real / 0.01) -- so the
    sample is the entire population, exactly as it is for E/F, where stage 2
    takes all of SLiM's individuals. That is not a coincidence to preserve by
    hand: if gwas_size or the GTEx fan-out ever exceeds 9,000, the coalescent
    would be asked to over-sample its own terminal population, which Hudson
    permits but which no longer approximates the forward arms.
    """
    gwas_size = int(cfg["gwas_size"])
    gtex_size = int(cfg.get("gtex_size", -1))
    largest_gtex = 1000 if gtex_size == -1 else gtex_size
    n = gwas_size + largest_gtex
    terminal = cattle_demography.terminal_size(float(cfg["Q_scaling"]))
    if n > terminal:
        raise ValueError(
            f"cattle_neutral would sample {n} individuals from a terminal population "
            f"of {terminal:,.0f} (epoch 12 = {cattle_demography.EP_SIZES[12]} at "
            f"Q_scaling={cfg['Q_scaling']}). E/F sample their whole SLiM population and "
            f"no more; over-sampling here would stop approximating them."
        )
    return n


rule stage1_cattle_neutral:
    output:
        ts = paths.stage1_dir(config) + "/" + paths.stage1_cattle_neutral_ts(config),
    input:
        existing = lambda wc: search_dirs.find_or_empty(
            paths.stage1_cattle_neutral_ts(config),
            config.get("stage1_search_dirs", []),
            expected_seed=config["stage1_seed"],
        ),
    params:
        seed       = config["stage1_seed"],
        L          = config["L"],
        recomb     = config["recombination_rate"],
        # Required by the Snakefile for this pipeline, so no .get() default. NOTE
        # this is 5.6e-9, the stage-1 COMPONENT -- not the 1.4e-8 total category H
        # asks for. Stage 2's add_neutral overlays the remaining 8.4e-9 on the
        # cattle branch whether this arm wants it or not, so reproducing E/F/G's
        # split is what lands the total on 1.4e-8 without touching shared code.
        mut_rate   = config["mutation_rate"],
        Q          = config["Q_scaling"],
        n_samples  = _cattle_neutral_n_samples(config),
        script_dir = SIM_REPO_DIR,
        script     = "cattle_neutral_simulation.py",
        python     = config["python_binary"],
        tmp_dir    = os.path.join(paths.workdir(config), "cattle_neutral_tmp"),
        produced_ts = os.path.join(paths.workdir(config), "cattle_neutral_tmp",
                                   f"cnts_{config['stage1_seed']}.ts"),
    log:
        os.path.join(paths.workdir(config), "logs", "stage1_cattle_neutral.log"),
    resources:
        # Measured at the production size (L = 2 Mb, 9,000 individuals, Q = 0.01):
        # 1.9 s wall, 103 MB peak RSS. The ancestral population is 1.7e6 at this Q,
        # which sounds alarming and is not -- rho = 4*N*r*L is invariant under the
        # rescaling, so the work is the same as at N = 17,000 in generations. Left
        # at the human_neutral tier rather than trimmed to the measurement: the
        # margin costs nothing on `short` and absorbs a larger L.
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
            --Q {params.Q} \
            --tmp_dir "{params.tmp_dir}" \
            > {log} 2>&1
        cp -u "{params.produced_ts}" "{output.ts}"
        """
