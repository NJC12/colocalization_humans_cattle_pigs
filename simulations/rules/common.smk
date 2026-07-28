# Stages 2, 3, 4, 5 -- shared across all four pipelines.
#
# Snakemake rule contract: `output:` strings must use only `{wildcards}`
# (resolved by Snakemake itself), so we use the `*_tmpl()` helpers from
# helpers/paths.py for outputs. `input:`, `params:`, and `log:` accept
# lambdas, so we use the concrete (cat-aware) helpers there.


# ----------------------------------------------------------------------------
# Stage 2: GWAS/GTEx split, trait simulation, PCA
# ----------------------------------------------------------------------------
#
# Wraps create_gwas_files_and_phenotypes.py. The script writes its outputs
# under {out_dir}/gwas_{gwas_scaling}_gtex_{gtex_scaling}_maf_{min_maf}/ --
# we mirror that exactly. Output declared as a touch marker since the script
# writes ~20 files in a fixed naming pattern.

def _stage2_outputs(cfg):
    """Files create_gwas_files_and_phenotypes.py writes that downstream rules
    explicitly depend on. We don't enumerate every file (the script writes
    20-some), only the ones consumed by stage 3 / 4 / 5."""
    out = [paths.stage2_marker(cfg)]
    for cat in paths.stage2_categories(cfg):
        out += [
            paths.stage2_pheno(cfg, cat),
            paths.stage2_geno(cfg, cat),
            paths.stage2_vcf(cfg, cat),
            paths.stage2_pheno_for_plink(cfg, cat),
        ]
    return out


rule stage2_split_pheno:
    output:
        _stage2_outputs(config),
    input:
        ts    = lambda wc: search_dirs.find_or_canonical(
            paths.stage1_full_filename(config),
            config.get("stage1_search_dirs", []),
            paths.stage1_dir(config),
            expected_seed=config["stage1_seed"],
        ),
        marks = lambda wc: (
            search_dirs.find_or_canonical(
                paths.stage1_marks_filename(config),
                config.get("stage1_search_dirs", []),
                paths.stage1_dir(config),
                expected_seed=config["stage1_seed"],
            ) if paths.stage1_marks_filename(config) else []
        ),
    params:
        ts_flag         = "--cattle_ts_file" if config["species"] == "cattle" else "--human_ts_file",
        m4_flag         = lambda wc, input: (
            f"--cattle_m4_file {input.marks}" if config["pipeline"] == "cattle_sel" and input.marks else ""
        ),
        gwas_scaling       = config["gwas_scaling"],
        gtex_scaling       = config["gtex_scaling"],
        # Stage-1 SLiM rescaling factor. Both engines record selection
        # coefficients as Q_scaling * the real value; get_vars_df divides them
        # back, and beta = sqrt(|selco|) * scaling depends on it. No
        # config.setdefault() for this one on purpose -- a silent default of 1
        # would reintroduce the un-corrected effect sizes, so a config missing
        # the key must fail loudly.
        Q_scaling          = config["Q_scaling"],
        # Cattle split-Q deep history: when stage 1 ran the burn-in and epochs
        # 2-7 at Q=1 and only epochs 8-12 at Q_scaling, the tree sequence has a
        # piecewise time scale and both the neutral overlay and the time column
        # have to be applied in two segments. Omit both keys for a single-Q
        # stage 1 (the Q=0.01 hedge, and the human arm) and stage 2 behaves as
        # before.
        handoff_flag       = (
            f"--handoff_ticks {config['handoff_ticks']} "
            f"--deep_Q_scaling {config['deep_Q_scaling']}"
            if config.get("handoff_ticks") is not None else ""
        ),
        min_maf            = config["min_maf"],
        L                  = config["L"],
        gtex_size          = config["gtex_size"],
        already_neutral    = config.get("already_includes_neutral", True),
        neutral_trait_vars = config.get("neutral_trait_vars", False),
        # Number of central trait loci (GWAS + shared GTEx). Omit the flag when
        # unset so the script falls back to "use all eligible" (legacy).
        central_flag       = (
            f"--n_central_traits {config['n_central_traits']}"
            if config.get("n_central_traits") is not None else ""
        ),
        # GTEx-only flank loci; 50 for all runs unless a config overrides it.
        n_flank_gtex       = config.get("n_flank_gtex_traits", 50),
        # Cap the total sampled pool (GWAS + GTEx) at n_samples individuals.
        # Omitted when unset so the script defaults to the full simulated
        # population (legacy). Used to cap the no-bottleneck cattle category (G),
        # whose ~100k population would otherwise give a ~99k-individual GWAS
        # instead of the intended 8k (bottlenecked E/F have ~9k pop -> 8k GWAS).
        n_samples_flag     = (
            f"--n_samples {config['n_samples']}"
            if config.get("n_samples") is not None else ""
        ),
        out_dir            = paths.stage2_dir(config),
        marker             = paths.stage2_marker(config),
        seed               = config["stage2_seed"],
        python             = config["python_binary"],
        script             = os.path.join(SIM_REPO_DIR, "create_gwas_files_and_phenotypes.py"),
    log: os.path.join(paths.workdir(config), "logs", "stage2_split_pheno.log")
    resources:
        # Sized from 6 measured round-3 runs (A1/E1 x x35/x10/x20 at L=2 Mb):
        # max elapsed 2:17, max RSS 369 MB. The previous priority/24h/128 GB was
        # carried over from L=10 Mb, where SBAMS/VCFs are ~10x larger.
        #
        # The 24h request is why controller_2Mb*.sbatch had to carry
        # `--set-resources stage2_split_pheno:slurm_partition=medium`: 24h does
        # not fit `short`'s 12h ceiling. At 1h it does, so the rule declares
        # `short` directly and the overrides are gone.
        #
        # If L=10 Mb is ever revived, raise these again -- they are calibrated
        # for 2 Mb only.
        slurm_partition = "short",
        runtime = 60,
        mem_mb  = 8000,
        cpus_per_task = 4,
    conda: "../envs/coloc_sims.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.out_dir}" "$(dirname {log})"
        "{params.python}" "{params.script}" \
            {params.ts_flag} "{input.ts}" \
            {params.m4_flag} \
            --gwas_scaling {params.gwas_scaling} \
            --gtex_scaling {params.gtex_scaling} \
            --Q_scaling {params.Q_scaling} \
            {params.handoff_flag} \
            --min_maf {params.min_maf} \
            --length {params.L} \
            --r2_value 0.2 \
            --gtex_size {params.gtex_size} \
            --seed {params.seed} \
            --out_dir "{params.out_dir}" \
            --already_includes_neutral {params.already_neutral} \
            --neutral_trait_vars {params.neutral_trait_vars} \
            {params.central_flag} \
            --n_flank_gtex_traits {params.n_flank_gtex} \
            {params.n_samples_flag} \
            > {log} 2>&1
        touch "{params.marker}"
        """


# ----------------------------------------------------------------------------
# Stage 3: DAP-G fine-mapping (native rules; checkpoint fans out per trait)
# ----------------------------------------------------------------------------

rule stage3_index_geno:
    output:
        gz  = paths.stage3_geno_gz_tmpl(config),
        tbi = paths.stage3_geno_tbi_tmpl(config),
    input:
        marker = ancient(paths.stage2_marker(config)),
        geno   = lambda wc: ancient(paths.stage2_geno(config, wc.cat)),
    params:
        bgzip      = config.get("bgzip_binary", "bgzip"),
        tabix      = config.get("tabix_binary", "tabix"),
        htslib_lib = config.get("htslib_lib", ""),
    log: os.path.join(paths.workdir(config), "logs", "stage3_index_{cat}.log")
    resources:
        # Measured over 18 round-3 jobs: max 1:47, max RSS 110 MB.
        slurm_partition = "short",
        runtime = 30,
        mem_mb  = 4000,
        cpus_per_task = 8,
    conda: "../envs/coloc_sims.yaml"
    shell:
        r"""
        set -euo pipefail
        if [ -n "{params.htslib_lib}" ]; then
            export LD_LIBRARY_PATH="{params.htslib_lib}:${{LD_LIBRARY_PATH:-}}"
        fi
        mkdir -p "$(dirname {output.gz})" "$(dirname {log})"
        # Prepend chrom (1) and numeric position parsed from snp_id, sort by
        # position, bgzip, tabix-index. Mirrors submit_revision_dapg_o2.sh.
        awk 'BEGIN{{OFS="\t"}} {{pos=$2; sub(/^snp/,"",pos); print "1", pos, $0}}' "{input.geno}" \
            | sort -k2,2n -S 4G --parallel={resources.cpus_per_task} \
            | "{params.bgzip}" -@ {resources.cpus_per_task} > "{output.gz}.tmp" 2> {log}
        mv "{output.gz}.tmp" "{output.gz}"
        "{params.tabix}" -s 1 -b 2 -e 2 "{output.gz}"
        """


checkpoint stage3_manifest:
    output:
        manifest = paths.stage3_manifest_tmpl(config),
    input:
        marker = ancient(paths.stage2_marker(config)),
        pheno  = lambda wc: ancient(paths.stage2_pheno(config, wc.cat)),
    run:
        # Default: emit all trait IDs (column 2 of the pheno sbams file).
        # If `selection_cap` is set, each replicate independently picks a
        # ceil(selection_cap / num_replicates) prefix of a (selection_seed,
        # stage1_seed)-seeded permutation -- see helpers/selection.py.
        os.makedirs(os.path.dirname(output.manifest), exist_ok=True)
        with open(input.pheno) as f:
            all_traits = [ln.split()[1] for ln in f if ln.strip()]
        cap = config.get("selection_cap")
        if cap is not None:
            from helpers.selection import select_traits_for_replicate
            selected = select_traits_for_replicate(
                all_traits,
                cap=int(cap),
                seed=int(config["selection_seed"]),
                num_replicates=int(config["num_replicates"]),
                stage1_seed=int(config["stage1_seed"]),
            )
        else:
            selected = all_traits
        with open(output.manifest, "w") as f:
            for t in selected:
                f.write(t + "\n")


def _dapg_mem_mb_base(wc):
    # Per-category base memory.
    #
    # MEASURED over 1500 round-3 DAP-G jobs (A1 human + E1 cattle-baseline, at
    # x35/x10/x20, L=2 Mb, PCA removed, ±0.25 Mb window). Peak RSS by category:
    #     hgwas 2580 MB   hgtex 249 MB   hgtex_smaller 112 MB
    #     cgwas  633 MB   cgtex 113 MB   cgtex_smaller 111 MB
    # The old numbers came from L=1 Mb ±1 Mb-window runs that still carried PCA
    # covariates, where hgwas RSS was 16-40 GB. Both changes cut the SBAMS size.
    cat = wc.cat
    if config["species"] == "human":
        if cat.endswith("gwas"):
            return 8000        # 3.1x measured peak; retries give 16/24 GB
        if cat == "hgtex":
            return 4000        # 16x measured peak
        return 4000            # hgtex_smaller, hgtex_smallest: 36x measured peak
    # cattle cgwas memory. DELIBERATELY NOT REDUCED, despite E1 measuring only
    # 633 MB. The 1500-job round-3 sample covers category E only; F and G have
    # not been rerun. Round-2 history: F1 (cattle_sel_bottlenecked) peaked at
    # ~947 MB, but G1-G4 (cattle_sel_not_bottlenecked, continue_bottlenecking=0)
    # OOM-killed 100% of cgwas DAP-G jobs at BOTH 4 GB and 12 GB, peaking near
    # 30 GB -- G's higher Ne means more common variants in the 8000-sample
    # cgwas and a much larger sparse-LD matrix. E is not a safe proxy for G.
    # Revisit once F/G have run under round-3 code; until then 32 GB stands,
    # scaling 32 -> 64 -> 96 GB on retries.
    if cat.endswith("gwas"):
        return 32000
    return 4000                # cgtex, cgtex_smaller: 35x measured peak


def _dapg_mem_mb(wc, attempt):
    # attempt-scaled: 1x base, 2x base, 3x base on retries (paired with
    # --restart-times in profiles/o2/config.yaml). Lets the upper tail of
    # heavy loci succeed without over-provisioning the median job.
    return _dapg_mem_mb_base(wc) * attempt


rule stage3_dapg_locus:
    output:
        out = paths.stage3_locus_out_tmpl(config),
        log_out = paths.stage3_locus_log_tmpl(config),
    input:
        pheno = lambda wc: paths.stage2_pheno(config, wc.cat),
        geno  = lambda wc: paths.stage3_geno_gz(config, wc.cat),
        tbi   = lambda wc: paths.stage3_geno_tbi(config, wc.cat),
        # No `pca` input: DAP-G runs without covariates. See run_plink_glm.sh.
    params:
        ld_ctrl      = config.get("ld_ctrl", 0.75),
        window       = config.get("dapg_window", 1_000_000),
        dapg         = config.get("dapg_binary",  "/n/data2/hms/dbmi/sunyaev/lab/nconnally/bin/dap/dap_src/dap-g"),
        tabix        = config.get("tabix_binary", "/n/data2/hms/dbmi/sunyaev/lab/nconnally/bin/htslib-1.23.1/tabix"),
        htslib_lib   = config.get("htslib_lib",   "/n/data2/hms/dbmi/sunyaev/lab/nconnally/bin/htslib-1.23.1"),
        dapg_libpath = config.get("dapg_libpath", ""),
        script       = os.path.join(SIM_REPO_DIR, "scripts", "dapg_one_locus.sh"),
    resources:
        # DAP-G is single-threaded (`seff` CPU efficiency = 6.08-6.20% on
        # 16-CPU jobs, i.e. ~1 of 16 cores used). cpus_per_task=4 gives
        # mild headroom for transient bursts without burning fairshare on
        # idle cores or blocking on 16-CPU slot availability. Memory is
        # attempt-scaled (32 -> 64 -> 96 GB for hgwas) to handle the upper
        # tail without over-provisioning the median.
        # Runtime measured over 1500 round-3 jobs (A1/E1 x x35/x10/x20, L=2 Mb).
        # The tail is entirely hgwas -- median 19:19, max 23:50 -- because human
        # carries ~36k variants per region against cattle's ~1.4k. Every other
        # category maxes at 6:29:
        #     hgwas  med 19:19  max 23:50      cgwas         med 2:20  max 4:44
        #     hgtex  med  1:22  max  3:02      cgtex         med 1:37  max 5:37
        #     hgtex_smaller med 1:47 max 5:19  cgtex_smaller med 1:53  max 6:29
        # 60 min for *gwas is 2.5x the observed hgwas max; 30 min for *gtex is
        # 4.6x the observed gtex max. A walltime kill is retried by
        # --restart-times, so the failure mode here is recoverable.
        slurm_partition = "short",
        runtime         = lambda wc: 60 if wc.cat.endswith("gwas") else 30,
        mem_mb          = _dapg_mem_mb,
        cpus_per_task   = 4,
    shell:
        r"""
        bash "{params.script}" \
            --pheno "{input.pheno}" \
            --geno  "{input.geno}" \
            --trait "{wildcards.trait}" \
            --window {params.window} \
            --ld-ctrl {params.ld_ctrl} \
            --dapg "{params.dapg}" \
            --tabix "{params.tabix}" \
            --htslib-lib "{params.htslib_lib}" \
            --dapg-libpath "{params.dapg_libpath}" \
            --threads {resources.cpus_per_task} \
            --out "{output.out}" \
            --log "{output.log_out}"
        """


def _stage3_targets(wildcards):
    chk = checkpoints.stage3_manifest.get(cat=wildcards.cat).output.manifest
    with open(chk) as f:
        traits = [line.strip() for line in f if line.strip()]
    return [paths.stage3_locus_out(config, wildcards.cat, t) for t in traits]


rule stage3_done:
    output: touch(paths.stage3_done_tmpl(config))
    input: _stage3_targets


# ----------------------------------------------------------------------------
# Stage 4: fastEnloc
# ----------------------------------------------------------------------------
#
# Per-SNP annotation VCF derived from the stage-2 genotype VCF. The awk recipe
# is from create_gwas_files_and_phenotypes.py (the helper .txt it writes).

rule stage4_annot_vcf:
    output:
        annot = paths.stage2_inner(config) + "/{cat}.dap.vcf.gz",
    input:
        marker = ancient(paths.stage2_marker(config)),
        vcf    = lambda wc: ancient(paths.stage2_vcf(config, wc.cat)),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.annot})"
        if [[ "{input.vcf}" == *.gz ]]; then
            zcat "{input.vcf}"
        else
            cat "{input.vcf}"
        fi | awk 'BEGIN{{OFS="\t"}} $1 == 1 {{$3="snp"$2; print $0}}' | gzip > "{output.annot}"
        """


def _stage4_inputs(wildcards):
    gwas = paths.stage2_gwas_category(config)
    return {
        "gwas_done":  paths.stage3_done(config, gwas),
        "gtex_done":  paths.stage3_done(config, wildcards.gtex_cat),
        "annot_gwas": paths.stage2_dap_annot_vcf(config, gwas),
        "annot_gtex": paths.stage2_dap_annot_vcf(config, wildcards.gtex_cat),
    }


rule stage4_fastenloc:
    output:
        outs = paths.stage4_outputs_tmpl(config),
    input:
        unpack(_stage4_inputs),
    params:
        gwas_dir   = lambda wc: paths.stage3_dir(config, paths.stage2_gwas_category(config)),
        gtex_dir   = lambda wc: paths.stage3_dir(config, wc.gtex_cat),
        prefix     = lambda wc: paths.stage4_prefix(config, wc.gtex_cat),
        outputs_subdir = paths.stage3_outputs_subdir_name(config),
        # "auto" means: count IDs that appear in both annotation VCFs (the
        # intersection of GWAS and GTEx tested SNPs). Override by setting
        # fastenloc_total_snps to an integer in the config YAML. Read only when
        # fastenloc_prior == "estimated"; the default prior ignores it.
        total_snps = config.get("fastenloc_total_snps", "auto"),
        prior      = config.get("fastenloc_prior", "coloc_default"),
        dap2enloc  = config["dap2enloc_binary"],
        fastenloc  = config["fastenloc_binary"],
        fastenloc_libpath = config.get("fastenloc_libpath", ""),
        script     = os.path.join(SIM_REPO_DIR, "scripts", "run_fastenloc.sh"),
    log: os.path.join(paths.workdir(config), "logs", "stage4_fastenloc_{gtex_cat}.log")
    resources:
        # The 32 GB request existed for the "auto" -total_variants computation,
        # which sorts and intersects both gunzipped annotation VCFs (~1.2 M IDs
        # for human at L=10 Mb). The default prior skips that entirely, so the
        # job is now just dap2enloc plus fastEnloc itself.
        # Measured over 12 round-3 jobs under the default prior: max 1:50
        # elapsed, 115 MB RSS. The "estimated" branch still gets 32 GB because
        # its -total_variants intersection has not been measured under round-3
        # code, and it is the reason the 32 GB tier existed at all.
        slurm_partition = "short",
        runtime = (4 * 60 if config.get("fastenloc_prior", "coloc_default") == "estimated"
                   else 30),
        mem_mb  = (32000 if config.get("fastenloc_prior", "coloc_default") == "estimated"
                   else 4000),
        cpus_per_task = 4,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        bash "{params.script}" \
            --gwas-dir "{params.gwas_dir}" \
            --gtex-dir "{params.gtex_dir}" \
            --outputs-subdir "{params.outputs_subdir}" \
            --annot-gwas "{input.annot_gwas}" \
            --annot-gtex "{input.annot_gtex}" \
            --out-prefix "{params.prefix}" \
            --prior "{params.prior}" \
            --total-snps {params.total_snps} \
            --dap2enloc "{params.dap2enloc}" \
            --fastenloc "{params.fastenloc}" \
            --fastenloc-libpath "{params.fastenloc_libpath}" \
            > {log} 2>&1
        """


# ----------------------------------------------------------------------------
# Stage 5: plink GLM (one rule instance per category)
# ----------------------------------------------------------------------------

rule stage5_plink_glm:
    output:
        done = touch(paths.stage5_done_tmpl(config)),
    input:
        marker = ancient(paths.stage2_marker(config)),
        vcf    = lambda wc: ancient(paths.stage2_vcf(config, wc.cat)),
        pheno  = lambda wc: ancient(paths.stage2_pheno_for_plink(config, wc.cat)),
    params:
        prefix = paths.stage5_prefix_tmpl(config),
        maf    = config["min_maf"],
        plink2 = config.get("plink2_binary", "plink2"),
        script = os.path.join(SIM_REPO_DIR, "scripts", "run_plink_glm.sh"),
    log: os.path.join(paths.workdir(config), "logs", "stage5_plink_{cat}.log")
    resources:
        # Measured over 18 round-3 jobs at L=2 Mb: max 0:46 elapsed, 110 MB RSS.
        # The 8h/32 GB figures were for L=10 Mb / 1.2 M SNPs, and also predate
        # the PCA removal that dropped the --indep-pairwise and --pca steps this
        # rule used to run before the GLM.
        slurm_partition = "short",
        runtime = 30,
        mem_mb  = 4000,
        cpus_per_task = 4,
    conda: "../envs/coloc_sims.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.done})" "$(dirname {log})"
        bash "{params.script}" \
            --vcf "{input.vcf}" \
            --pheno "{input.pheno}" \
            --out-prefix "{params.prefix}" \
            --maf {params.maf} \
            --plink2 "{params.plink2}" \
            > {log} 2>&1
        """
