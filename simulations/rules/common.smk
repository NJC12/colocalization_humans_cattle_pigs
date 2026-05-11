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
            paths.stage2_pca(cfg, cat),
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
        min_maf            = config["min_maf"],
        gtex_size          = config["gtex_size"],
        already_neutral    = config.get("already_includes_neutral", True),
        neutral_trait_vars = config.get("neutral_trait_vars", False),
        out_dir            = paths.stage2_dir(config),
        marker             = paths.stage2_marker(config),
        seed               = config["stage2_seed"],
        python             = config["python_binary"],
        script             = os.path.join(SIM_REPO_DIR, "create_gwas_files_and_phenotypes.py"),
    log: os.path.join(paths.workdir(config), "logs", "stage2_split_pheno.log")
    resources:
        # On O2, jobs <12h must go to `short` (or `priority`); `medium` is
        # reserved for jobs longer than 12h. Set 8h here -- the trait sim on a
        # 10 Mb cattle TS finishes in well under that.
        slurm_partition = "short",
        runtime = 8 * 60,
        mem_mb  = 64000,
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
            --min_maf {params.min_maf} \
            --r2_value 0.2 \
            --gtex_size {params.gtex_size} \
            --seed {params.seed} \
            --out_dir "{params.out_dir}" \
            --already_includes_neutral {params.already_neutral} \
            --neutral_trait_vars {params.neutral_trait_vars} \
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
        slurm_partition = "short",
        runtime = 4 * 60,
        mem_mb  = 32000,
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
    shell:
        r"""
        mkdir -p "$(dirname {output.manifest})"
        awk '{{print $2}}' "{input.pheno}" > "{output.manifest}"
        """


rule stage3_dapg_locus:
    output:
        out = paths.stage3_locus_out_tmpl(config),
        log_out = paths.stage3_locus_log_tmpl(config),
    input:
        pheno = lambda wc: paths.stage2_pheno(config, wc.cat),
        geno  = lambda wc: paths.stage3_geno_gz(config, wc.cat),
        tbi   = lambda wc: paths.stage3_geno_tbi(config, wc.cat),
        pca   = lambda wc: ancient(paths.stage2_pca(config, wc.cat)),
    params:
        ld_ctrl      = config.get("ld_ctrl", 0.75),
        window       = config.get("dapg_window", 1_000_000),
        dapg         = config.get("dapg_binary",  "/n/data2/hms/dbmi/sunyaev/lab/nconnally/bin/dap/dap_src/dap-g"),
        tabix        = config.get("tabix_binary", "/n/data2/hms/dbmi/sunyaev/lab/nconnally/bin/htslib-1.23.1/tabix"),
        htslib_lib   = config.get("htslib_lib",   "/n/data2/hms/dbmi/sunyaev/lab/nconnally/bin/htslib-1.23.1"),
        dapg_libpath = config.get("dapg_libpath", ""),
        script       = os.path.join(SIM_REPO_DIR, "scripts", "dapg_one_locus.sh"),
    resources:
        # Per-category sizing -- cgtex (~1k samples) is light; cgwas (~99k
        # samples) needs 200 GB / 16 CPUs / 2-day runtime per the legacy
        # submit_revision_dapg_o2.sh "large" profile (sample size >80k).
        # Mixing them in one rule with one resource block OOM-killed every
        # cgwas job at 32 GB (peak RSS during dap-g matrix work for a 99k x
        # ~5k window blew past the limit). Using callables that branch on
        # wildcards.cat: gwas-category → large; gtex-category → small.
        slurm_partition = lambda wc: "medium" if wc.cat.endswith("gwas") else "short",
        runtime         = lambda wc: 48 * 60  if wc.cat.endswith("gwas") else 4 * 60,
        mem_mb          = lambda wc: 200000   if wc.cat.endswith("gwas") else 8000,
        cpus_per_task   = lambda wc: 16       if wc.cat.endswith("gwas") else 4,
    shell:
        r"""
        bash "{params.script}" \
            --pheno "{input.pheno}" \
            --geno  "{input.geno}" \
            --pca   "{input.pca}" \
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
        # "auto" means: count IDs that appear in both annotation VCFs (the
        # intersection of GWAS and GTEx tested SNPs). Override by setting
        # fastenloc_total_snps to an integer in the config YAML.
        total_snps = config.get("fastenloc_total_snps", "auto"),
        dap2enloc  = config["dap2enloc_binary"],
        fastenloc  = config["fastenloc_binary"],
        fastenloc_libpath = config.get("fastenloc_libpath", ""),
        script     = os.path.join(SIM_REPO_DIR, "scripts", "run_fastenloc.sh"),
    log: os.path.join(paths.workdir(config), "logs", "stage4_fastenloc_{gtex_cat}.log")
    resources:
        slurm_partition = "short",
        runtime = 2 * 60,
        mem_mb  = 16000,
        cpus_per_task = 4,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        bash "{params.script}" \
            --gwas-dir "{params.gwas_dir}" \
            --gtex-dir "{params.gtex_dir}" \
            --annot-gwas "{input.annot_gwas}" \
            --annot-gtex "{input.annot_gtex}" \
            --out-prefix "{params.prefix}" \
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
        slurm_partition = "short",
        runtime = 4 * 60,
        mem_mb  = 16000,
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
