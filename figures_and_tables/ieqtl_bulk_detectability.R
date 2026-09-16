#!/usr/bin/env Rscript
# =============================================================================
# ieqtl_bulk_detectability.R
#
# WHY THIS FILE EXISTS
# --------------------
# Cell-type interaction eQTLs (ieQTLs) find variants whose effect on expression
# depends on tissue cell composition. The obvious question is how much NEW signal
# that buys: what fraction of significant ieQTLs would an ordinary bulk eQTL scan
# have found anyway?
#
# This script answers that for human (GTEx v8) and pig (PigGTEx v0), and asks
# whether the two species answer it the same way -- overall, and stratified by
# distance to the TSS and by minor allele frequency.
#
#   Q1  Are human and pig ieQTLs detectable in bulk tissue at similar rates?
#   Q2  Does the human-vs-pig similarity change with TSS distance?
#   Q3  Does it change with MAF?
#
# THE ESTIMAND, STATED PRECISELY
# ------------------------------
# For each variant-gene pair that is a significant cell-type interaction eQTL,
# does THAT SAME PAIR also reach statistical significance as an ordinary bulk
# cis-eQTL? Two things this deliberately is NOT:
#
#   * NOT "is this gene a bulk eGene". A gene can be an eGene through a
#     completely independent signal elsewhere in the cis-window; that says
#     nothing about whether the cell-type-specific eQTL is visible in bulk.
#     (Scoring it that way gives 86.0% human / 53.1% pig -- a different and,
#     for this purpose, wrong question.)
#   * NOT the fraction that TRULY have a non-null bulk effect. That is a
#     statement about truth rather than detection, and would be estimated with
#     Storey pi1 (~85% human / ~47% pig). The question here is explicitly about
#     detectability at a significance threshold, so power differences between
#     the two studies are part of what is being measured, not a confound to be
#     removed. Interpret accordingly: pig tissues run n = 73-501 against human
#     n = 65-670, and the median |z| of b_g is 7.42 in human vs 1.37 in pig.
#
# WHY pval_g IS THE BULK EFFECT (this is the load-bearing fact)
# ------------------------------------------------------------
# Both projects fit, per variant-gene pair,
#
#       expression ~ g + i + g:i + covariates
#
# where i is the cell-type score. Centering design columns does not change OLS
# slopes (Frisch-Waugh-Lovell), so b_g is the genotype effect AT i = 0. Both
# projects inverse-normal transform the cell-type score before fitting, which
# puts i = 0 at the median cell composition. So b_g / pval_g IS the bulk-tissue
# effect at average cell composition -- no external bulk eQTL file is needed, and
# the two species are measured on the same footing.
#
# Confirmed empirically: |b_g| vs GTEx v8 bulk slope, r = 0.9985 over 800 shared
# variant-gene pairs. See notebook/entries/2026-09-16-gtex-ieqtl-file-format.md.
#
# THE DATA
# --------
#   HUMAN  43 files, flat, 21 columns (GTEx post-processed "eigenMT.annotated")
#          <Tissue>.<CellType>.ieQTL.eigenMT.annotated.txt.gz
#          xCell enrichment scores, 7 broad cell types, MAF floor 0.05
#   PIG    47 files, one subdirectory per tissue, 18 columns (raw tensorQTL)
#          <Tissue>/<Tissue>.interactions_pro_<CellType>.cis_qtl_top_assoc.txt.gz
#          single-cell deconvolved proportions, 20 fine types, MAF floor 0.10
#
# Both are one row per gene (top association by |t_gi| within the cis-window),
# both tensorQTL + eigenMT + Benjamini-Hochberg, both +/- 1 Mb cis-windows.
# Fifteen statistics columns are identically named and directly comparable.
#
# WHAT IS COMPARABLE, AND WHY (all verified, and re-asserted at runtime)
# ---------------------------------------------------------------------
#   * pval_emt = min(pval_gi * tests_emt, 1) holds exactly in both species.
#   * pval_adj_bh is a BH step-up over pval_emt WITHIN each file in both species
#     -- including the rows where pval_emt saturates at 1. The FDR columns are
#     therefore constructed identically and can be thresholded the same way.
#   * median tests_emt is 1265 (human) vs 1289 (pig), so the eigenMT correction
#     used for the bulk call does NOT systematically favour either species. This
#     mattered enough to check: if pig had far more effective tests, the
#     eigenMT-corrected bulk rule would have penalised pig by construction.
#
# GOTCHAS THAT SILENTLY CORRUPT THE RESULT
# ----------------------------------------
#   * PIG FILENAMES CONTAIN DOTS INSIDE THE CELL TYPE NAME. "L2.3_IT" and
#     "CD8a+_ab_T.NK_cells" both do. Never split a pig filename on ".". Take the
#     tissue from the directory name and strip the literal prefix
#     "<Tissue>.interactions_pro_". Human tissue names contain "-" but never ".",
#     so human splits on the LAST dot safely.
#   * PIG CALLS ITS MAF COLUMN "af", AND IT IS ALREADY FOLDED (observed max 0.5,
#     min 0.101). It is not a raw ALT frequency despite the name. Asserted below,
#     because if a future release ships unfolded values the MAF axis silently
#     inverts for half the variants.
#   * tests_emt / pval_emt SIT IN THE OPPOSITE ORDER in the two schemas: human is
#     (..., pval_emt, tests_emt, pval_adj_bh) at 19/20/21, pig is
#     (..., tests_emt, pval_emt, pval_adj_bh) at 16/17/18. Read by NAME, never by
#     position.
#   * tss_distance IS SIGNED BY GENOME COORDINATE, NOT BY STRAND, in both
#     projects (variant_pos - TSS_pos, no strand flip). For a minus-strand gene
#     the biological reading inverts. Everything here uses abs_tss_dist so the
#     result never depends on remembering this.
#   * ieQTL EFFECT SIZES ARE MINOR-ALLELE REFERENCED while standard GTEx v8 bulk
#     slopes are ALT-referenced (~27% of pairs differ in sign; plink1.9
#     --make-bed sets A1 = minor). This analysis uses only p-values so it is
#     unaffected -- but b_g/b_gi are carried into the output table, so do not
#     reuse their SIGNS against another file without harmonizing alleles first.
#   * PIG GENE IDS ARE UNVERSIONED, HUMAN ARE VERSIONED. Human ENSG must have its
#     ".N" stripped before it is used as a key.
#
# DESIGN DECISIONS (all chosen by the user, see the plan)
# ------------------------------------------------------
#   Significant ieQTL   each species' own pval_adj_bh < 0.05 (primary),
#                       plus a rank-matched top-K sensitivity
#   Significant bulk    pval_g < the gene's own permutation-derived
#                       pval_nominal_threshold  (see below)
#   Genes               protein-coding only
#   MAF                 harmonized to >= 0.10 for every cross-species number;
#                       the human 0.05-0.10 bin is kept as a human-only stratum
#   Cell types          strict matched pairs only, plus a pooled all-pairs view
#   Orthologs           not used -- this compares RATES, not genes
#
# MATCHED CELL-TYPE PAIRS
#   lung_epithelium  human Lung.Epithelial_cells
#                <-> pig  Lung.Alveolar_Epithelial_Type_1 / _2
#   brain_neurons    human Brain_*.Neurons (13 regions pooled)
#                <-> pig  {Brain,Frontal_cortex,Hypothalamus}.L2.3_IT
# No other pair is defensible. Human Adipocytes/Hepatocytes/Keratinocytes/
# Myocytes have no pig counterpart (pig Liver is profiled only for infiltrating
# immune cells, not hepatocytes), and pig Blood carries no neutrophil term.
#
# THE POWER CAVEAT, WHICH IS CENTRAL AND NOT A FOOTNOTE
# -----------------------------------------------------
# Pig sample sizes run 73-501 and human 65-670, and they are NOT matched within a
# comparison: pig Lung is n=149 against human Lung n~515. Pig interaction
# p-values are also far more liberal (pig Blood/Monocytes has 26.5% of genes at
# FDR < 0.05; human tops out at 4.3%). A raw rate difference between species is
# therefore part biology and part study design. Three things in this script push
# back on that, and none of them fully solves it:
#   1. n_samples is reported next to every rate
#   2. the rank-matched sensitivity equalises how deep into each ranking we go
#   3. the adjusted model carries log10(n_samples) as a covariate
# Read the headline numbers with all three in view.
#
# HOW "SIGNIFICANT BULK EFFECT" IS CALLED, AND WHY THIS WAY
# ---------------------------------------------------------
# The ieQTL files ship a multiple-testing-corrected p-value for the INTERACTION
# term only (pval_adj_bh, which is BH over pval_emt = pval_gi * tests_emt).
# There is no corrected counterpart for pval_g, so one has to be supplied.
#
# An earlier version of this script used pval_g * tests_emt < 0.05 -- the same
# eigenMT burden, but with no BH step, so it was the more liberal of the two
# rules and not symmetric with the ieQTL call. It is now retained only as the
# `bulk_emt` diagnostic column.
#
# The primary rule is instead each gene's OWN permutation-derived nominal
# threshold, taken from the matching bulk cis-eQTL release:
#
#   human  GTEx v8      <Tissue>.v8.egenes.txt.gz        -> pval_nominal_threshold
#   pig    PigGTEx v0   <Tissue>.cis_qtl_fdr0.05.txt.gz  -> pval_nominal_threshold
#
# This is the threshold both consortia use to declare a variant-gene pair a
# significant cis-eQTL, so "detectable in bulk" here means exactly what it means
# in the parent projects. It is gene-specific (it absorbs cis-window size and
# variant density per gene) and it is derived by permutation rather than by a
# closed-form correction, so it is calibrated rather than approximate.
#
# Three properties of these bulk releases were checked rather than assumed:
#   * SAME DONORS. Bulk and ieQTL sample sizes match within each species and
#     tissue (pig Lung 149 in both, human Lung 515 in both), so within a species
#     the bulk test and the interaction test have the same n. The species
#     difference in power is between studies, not between the two tests.
#   * FULL COVERAGE. Every tissue in the ieQTL set has a bulk counterpart (35/35
#     human, 7/7 pig), and every ieGene row joins to a bulk record (asserted).
#   * ALL GENES, not just eGenes. Both files carry every tested gene, so the
#     threshold exists for genes that are not eGenes. This matters: the question
#     is about the VARIANT, so a gene need not be an eGene for its
#     cell-type-specific eQTL to be asked about.
#
# Switching from the eigenMT rule to this one moves the pooled rate from 68.6%
# to 69.8% in human and 9.3% to 12.0% in pig.
#
# PIG PROTEIN-CODING ANNOTATION
# -----------------------------
# The pig files carry no biotype column. pig_gene_biotype.tsv next to this script
# was built from Ensembl release 100, which resolves 21,626/21,626 (100%) of the
# pig ieQTL gene ids -- that exact coverage is what identifies release 100 as the
# annotation PigGTEx v0 used. Rebuild it with:
#
#   curl -sL -o pig100.gtf.gz \
#     https://ftp.ensembl.org/pub/release-100/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.100.gtf.gz
#   { printf 'gene_id\tgene_biotype\n'
#     gzcat pig100.gtf.gz | awk -F'\t' '$3=="gene"{
#       if (match($9,/gene_id "[^"]+"/))      g=substr($9,RSTART+9,RLENGTH-10)
#       if (match($9,/gene_biotype "[^"]+"/)) b=substr($9,RSTART+14,RLENGTH-15)
#       if (g!="" && b!="") print g"\t"b }' | sort -u
#   } > pig_gene_biotype.tsv
#
# Human protein-coding comes free from the files' own biotype column (GENCODE
# v26, the GTEx v8 annotation).
#
# OUTPUTS (written next to this script, or to [output_dir] if given)
#   ieqtl_bulk_detectability.tsv.gz   one row per significant ieGene, long,
#                                     species-keyed, with every statistic and
#                                     every stratum label
#   ieqtl_bulk_rates.tsv              per-stratum aggregate rates with Wilson CIs
#   ieqtl_bulk_models.tsv             the species x stratifier logistic models
#
# RUNTIME  ~2-3 minutes, ~2 GB peak. Reads 90 gzipped files.
#
# USAGE
#   Rscript ieqtl_bulk_detectability.R [output_dir]
# =============================================================================

suppressPackageStartupMessages(library(data.table))

options(width = 200)
hr <- function() cat(strrep("-", 104), "\n")

script_dir <- function() {
    ca  <- commandArgs(trailingOnly = FALSE)
    hit <- grep("^--file=", ca, value = TRUE)
    if (length(hit) == 0) return(normalizePath("."))
    normalizePath(dirname(sub("^--file=", "", hit[1])))
}

args   <- commandArgs(trailingOnly = TRUE)
HERE   <- script_dir()
OUTDIR <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else HERE

# ---- configuration ----------------------------------------------------------

HUMAN_DIR    <- path.expand("~/data/gtex_ieQTL/GTEx_Analysis_v8_ieQTL")
PIG_DIR      <- path.expand("~/data/gtex_ieQTL/PigGTEx_v0.Celltype_interaction_eQTL")
PIG_BIOTYPE  <- file.path(HERE, "pig_gene_biotype.tsv")

# Bulk cis-eQTL releases, for the per-gene permutation thresholds.
HUMAN_BULK_TAR <- path.expand("~/data/gtex_v8/GTEx_Analysis_v8_eQTL.tar")
PIG_BULK_DIR   <- path.expand(
    "~/eqtl_selection/data/gtex/pgtex/eqtls_permutation/PigGTEx_v0.permutations_eQTL")

N_HUMAN_FILES <- 43L
N_PIG_FILES   <- 47L   # 50 runs launched; 3 died upstream and shipped no result:
                       # Hypothalamus.Micro, Liver.CD8ab+_ab_T_cells,
                       # Spleen.CD8ab+_ab_T_cells

IE_FDR        <- 0.05   # pval_adj_bh threshold defining a significant ieQTL
BULK_ALPHA    <- 0.05   # threshold for pval_g * tests_emt
MAF_FLOOR     <- 0.10   # harmonized floor; pig was mapped at 0.10, human at 0.05
CIS_WINDOW    <- 1e6

TSS_BREAKS <- c(0, 1e4, 5e4, 1e5, 2.5e5, 5e5, 1e6)
TSS_LABELS <- c("<10kb", "10-50kb", "50-100kb", "100-250kb", "250-500kb", "500kb-1Mb")

MAF_BREAKS_FINE  <- seq(0.10, 0.50, by = 0.05)          # pooled comparison
MAF_BREAKS_COARSE <- c(0.10, 0.20, 0.30, 0.40, 0.50)    # matched pairs (fewer genes)

# Matched cell-type pairs. Human brain regions are pooled; pig L2.3_IT is the
# only neuronal pig type (Astro/Micro/Oligo are glia).
HUMAN_MATCH <- list(
    lung_epithelium = "Lung.Epithelial_cells",
    brain_neurons   = NULL   # filled below: every Brain_*.Neurons pair
)
PIG_MATCH <- list(
    lung_epithelium = c("Lung.Alveolar_Epithelial_Type_1",
                        "Lung.Alveolar_Epithelial_Type_2"),
    brain_neurons   = c("Brain.L2.3_IT",
                        "Frontal_cortex.L2.3_IT",
                        "Hypothalamus.L2.3_IT")
)

need_file <- function(p, what) {
    if (!file.exists(p)) stop(sprintf("%s not found: %s", what, p), call. = FALSE)
    p
}
need_dir <- function(p, what) {
    if (!dir.exists(p)) stop(sprintf("%s not found: %s", what, p), call. = FALSE)
    p
}

need_dir(HUMAN_DIR, "human ieQTL directory")
need_dir(PIG_DIR,   "pig ieQTL directory")
need_file(PIG_BIOTYPE, "pig biotype table")
need_file(HUMAN_BULK_TAR, "GTEx v8 bulk eQTL tar")
need_dir(PIG_BULK_DIR, "PigGTEx bulk eQTL permutation directory")

# ---- bulk permutation thresholds --------------------------------------------
# Returns one row per (species, tissue, gene) with that gene's permutation-
# derived nominal significance threshold for the bulk cis-eQTL scan.
#
# The human release lives inside a 1.5 GB tar. Extract every needed member in a
# SINGLE tar pass -- the same trick red_herring_eqtls.R uses -- because one pass
# per tissue would re-scan the archive 35 times.

read_bulk_thresholds <- function(human_tissues, pig_tissues) {
    tmp <- file.path(tempdir(), "gtex_v8_egenes")
    dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
    members <- sprintf("GTEx_Analysis_v8_eQTL/%s.v8.egenes.txt.gz", human_tissues)
    rc <- system2("tar", c("-xf", shQuote(HUMAN_BULK_TAR), "-C", shQuote(tmp),
                           shQuote(members)), stdout = NULL, stderr = NULL)
    if (rc != 0) stop("failed to extract egenes from ", HUMAN_BULK_TAR, call. = FALSE)

    h <- rbindlist(lapply(human_tissues, function(t) {
        f <- file.path(tmp, "GTEx_Analysis_v8_eQTL",
                       sprintf("%s.v8.egenes.txt.gz", t))
        need_file(f, sprintf("GTEx v8 egenes for %s", t))
        d <- fread(cmd = sprintf("gzcat %s", shQuote(f)), showProgress = FALSE,
                   select = c("gene_id", "qval", "pval_nominal_threshold"))
        # strip the version only -- these ids have no other suffix
        d[, .(species = "human", tissue = t, gene = sub("\\.[0-9]+$", "", gene_id),
              bulk_nom_thresh = pval_nominal_threshold, bulk_egene = qval < 0.05)]
    }))

    p <- rbindlist(lapply(pig_tissues, function(t) {
        f <- file.path(PIG_BULK_DIR, sprintf("%s.cis_qtl_fdr0.05.txt.gz", t))
        need_file(f, sprintf("PigGTEx permutation results for %s", t))
        d <- fread(cmd = sprintf("gzcat %s", shQuote(f)), showProgress = FALSE,
                   select = c("phenotype_id", "is_eGene", "pval_nominal_threshold"))
        d[, .(species = "pig", tissue = t, gene = phenotype_id,
              bulk_nom_thresh = pval_nominal_threshold,
              bulk_egene = is_eGene %in% c(TRUE, "TRUE"))]
    }))

    b <- rbind(h, p)
    if (anyDuplicated(b[, .(species, tissue, gene)]))
        stop("bulk threshold table has duplicate (species, tissue, gene) keys",
             call. = FALSE)
    b
}

# ---- shared column contract -------------------------------------------------
# Both readers must return exactly this, in this order. Nothing downstream is
# allowed to branch on species.

SHAPE <- c("species", "tissue", "cell_type", "pair", "gene", "variant_id",
           "tss_distance", "maf", "ma_samples", "ma_count",
           "pval_g", "b_g", "b_g_se", "pval_i", "b_i", "b_i_se",
           "pval_gi", "b_gi", "b_gi_se", "tests_emt", "pval_emt", "pval_adj_bh",
           "protein_coding")

# ---- integrity checks applied identically to both species -------------------

check_file_integrity <- function(d, label) {
    # one row per gene (these are top-association files)
    if (anyDuplicated(d$gene)) {
        stop(sprintf("%s: %d duplicated gene ids -- expected one row per gene",
                     label, sum(duplicated(d$gene))), call. = FALSE)
    }
    # pval_emt = min(pval_gi * tests_emt, 1)
    expect <- pmin(d$pval_gi * d$tests_emt, 1)
    bad <- sum(abs(expect - d$pval_emt) / pmax(d$pval_emt, 1e-300) > 1e-4)
    if (bad > 0) {
        stop(sprintf("%s: pval_emt != min(pval_gi*tests_emt,1) in %d/%d rows",
                     label, bad, nrow(d)), call. = FALSE)
    }
    # pval_adj_bh = BH step-up over pval_emt, over ALL rows including saturated
    o  <- order(d$pval_emt)
    n  <- nrow(d)
    bh <- rev(cummin(rev(d$pval_emt[o] * n / seq_len(n))))
    bh <- pmin(bh, 1)
    bad <- sum(abs(bh - d$pval_adj_bh[o]) / pmax(d$pval_adj_bh[o], 1e-300) > 1e-3)
    if (bad > 0) {
        stop(sprintf("%s: pval_adj_bh is not BH over pval_emt (%d/%d rows differ)",
                     label, bad, n), call. = FALSE)
    }
    # cis-window
    if (max(abs(d$tss_distance)) > CIS_WINDOW) {
        stop(sprintf("%s: |tss_distance| exceeds %g", label, CIS_WINDOW), call. = FALSE)
    }
    # MAF must be folded
    if (max(d$maf) > 0.5 + 1e-9) {
        stop(sprintf("%s: maf > 0.5 -- allele frequency is not folded", label),
             call. = FALSE)
    }
    invisible(TRUE)
}

# ---- human reader -----------------------------------------------------------
# Filenames: <Tissue>.<CellType>.ieQTL.eigenMT.annotated.txt.gz
# Human tissue names contain "-" (Brain_Spinal_cord_cervical_c-1) but never ".",
# so splitting on the last dot of the stem is safe.

read_human_ieqtl <- function() {
    files <- list.files(HUMAN_DIR, pattern = "\\.ieQTL\\.eigenMT\\.annotated\\.txt\\.gz$",
                        full.names = TRUE)
    if (length(files) != N_HUMAN_FILES) {
        stop(sprintf("expected %d human files, found %d in %s",
                     N_HUMAN_FILES, length(files), HUMAN_DIR), call. = FALSE)
    }
    out <- vector("list", length(files))
    for (k in seq_along(files)) {
        f    <- files[k]
        stem <- sub("\\.ieQTL\\.eigenMT\\.annotated\\.txt\\.gz$", "", basename(f))
        cut  <- regexpr("\\.[^.]+$", stem)          # last dot
        if (cut < 0) stop("cannot parse human filename: ", basename(f), call. = FALSE)
        tissue    <- substr(stem, 1, cut - 1)
        cell_type <- substr(stem, cut + 1, nchar(stem))

        d <- fread(cmd = sprintf("gzcat %s", shQuote(f)), showProgress = FALSE)
        check_file_integrity(d, paste("human", stem))

        d[, `:=`(species        = "human",
                 tissue         = tissue,
                 cell_type      = cell_type,
                 pair           = stem,
                 gene           = sub("\\.[0-9]+$", "", gene_id),   # strip version only
                 protein_coding = biotype == "protein_coding")]
        out[[k]] <- d[, ..SHAPE]
    }
    rbindlist(out)
}

# ---- pig reader -------------------------------------------------------------
# Filenames: <Tissue>/<Tissue>.interactions_pro_<CellType>.cis_qtl_top_assoc.txt.gz
# Cell type names contain "+", "-" AND "." -- take the tissue from the directory
# and strip the literal prefix. Never split on ".".

read_pig_ieqtl <- function() {
    files <- list.files(PIG_DIR, pattern = "\\.cis_qtl_top_assoc\\.txt\\.gz$",
                        full.names = TRUE, recursive = TRUE)
    if (length(files) != N_PIG_FILES) {
        stop(sprintf("expected %d pig files, found %d in %s",
                     N_PIG_FILES, length(files), PIG_DIR), call. = FALSE)
    }
    bt <- fread(PIG_BIOTYPE, showProgress = FALSE)
    setnames(bt, c("gene_id", "gene_biotype"))
    setkey(bt, gene_id)

    out <- vector("list", length(files))
    for (k in seq_along(files)) {
        f      <- files[k]
        tissue <- basename(dirname(f))
        stem   <- sub("\\.cis_qtl_top_assoc\\.txt\\.gz$", "", basename(f))
        prefix <- paste0(tissue, ".interactions_pro_")
        if (substr(stem, 1, nchar(prefix)) != prefix) {
            stop(sprintf("pig filename does not match <Tissue>.interactions_pro_<CellType>: %s",
                         basename(f)), call. = FALSE)
        }
        cell_type <- substr(stem, nchar(prefix) + 1, nchar(stem))

        d <- fread(cmd = sprintf("gzcat %s", shQuote(f)), showProgress = FALSE)
        # pig calls the folded MAF "af"; rename before the shared integrity check
        setnames(d, "af", "maf")
        setnames(d, "phenotype_id", "gene")
        check_file_integrity(d, paste("pig", tissue, cell_type))

        d[, `:=`(species   = "pig",
                 tissue    = tissue,
                 cell_type = cell_type,
                 pair      = paste0(tissue, ".", cell_type))]
        d[, protein_coding := bt[.(d$gene), gene_biotype] == "protein_coding"]
        if (anyNA(d$protein_coding)) {
            miss <- sum(is.na(d$protein_coding))
            stop(sprintf("pig %s.%s: %d/%d gene ids absent from %s",
                         tissue, cell_type, miss, nrow(d), basename(PIG_BIOTYPE)),
                 call. = FALSE)
        }
        out[[k]] <- d[, ..SHAPE]
    }
    rbindlist(out)
}

cat("\n"); hr()
cat("READING ieQTL FILES\n"); hr()

human <- read_human_ieqtl()
cat(sprintf("  human : %s rows across %d tissue x cell-type pairs\n",
            format(nrow(human), big.mark = ","), uniqueN(human$pair)))
pig <- read_pig_ieqtl()
cat(sprintf("  pig   : %s rows across %d tissue x cell-type pairs\n",
            format(nrow(pig), big.mark = ","), uniqueN(pig$pair)))

ie <- rbindlist(list(human, pig))
rm(human, pig); invisible(gc())

# ---- derived columns --------------------------------------------------------

ie[, abs_tss_dist := abs(tss_distance)]

# Sample size per pair, derived from the file itself: maf = ma_count / (2N).
# Validated against the published GTEx v8 n (Adipose_Subcutaneous -> 581).
ie[, n_samples := round(median(ma_count / (2 * maf), na.rm = TRUE)), by = .(species, pair)]

ie[, is_ieqtl := pval_adj_bh < IE_FDR]

# Attach each gene's permutation-derived bulk threshold, then call the bulk
# effect. The join is on (species, tissue, gene) and must be complete -- a
# missing bulk record would silently become a non-significant call.
bulk <- read_bulk_thresholds(
    human_tissues = sort(unique(ie[species == "human", tissue])),
    pig_tissues   = sort(unique(ie[species == "pig",   tissue])))
ie <- merge(ie, bulk, by = c("species", "tissue", "gene"), all.x = TRUE, sort = FALSE)

# The bulk scan tested FEWER genes than the ieQTL scan (e.g. pig Spleen: 14,909
# vs 18,951), so some ieQTL rows legitimately have no bulk record and get NA.
# That is fine for rows we never analyse, but a significant ieGene with no bulk
# counterpart could not be assessed at all and must not be silently scored as
# "not detectable". Coverage on the analysis set is 100% in both species today;
# this asserts it stays that way, and names the gap if it ever does not.
cov <- ie[protein_coding == TRUE & is_ieqtl == TRUE,
          .(n = .N, matched = sum(!is.na(bulk_nom_thresh))), by = species]
cov[, frac := matched / n]
if (any(cov$frac < 1)) {
    bad <- ie[protein_coding == TRUE & is_ieqtl == TRUE & is.na(bulk_nom_thresh),
              .N, by = .(species, tissue)][order(-N)]
    stop(sprintf("significant ieGenes with no bulk counterpart: %s. Worst: %s %s (%d). %s",
                 paste(sprintf("%s %d/%d", cov$species, cov$n - cov$matched, cov$n),
                       collapse = ", "),
                 bad$species[1], bad$tissue[1], bad$N[1],
                 "These cannot be assessed and must be excluded, not counted as negatives."),
         call. = FALSE)
}

ie[, bulk_sig := pval_g < bulk_nom_thresh]           # PRIMARY: permutation threshold
ie[, bulk_emt := pval_g * tests_emt < BULK_ALPHA]    # previous rule, diagnostic
ie[, bulk_nom := pval_g < BULK_ALPHA]                # uncorrected, diagnostic
ie[, bulk_pair := bulk_sig & bulk_egene]             # consortia's own signif-pair rule
ie[, z_g := abs(b_g) / b_g_se]                       # threshold-free power diagnostic

# rank within the interaction test, for the rank-matched sensitivity
ie[, ie_rank := frank(pval_gi, ties.method = "first"), by = .(species, pair)]

# matched cell-type groups
HUMAN_MATCH$brain_neurons <- sort(unique(
    ie[species == "human" & cell_type == "Neurons" &
       grepl("^Brain_", tissue), pair]))

ie[, matched_group := NA_character_]
for (g in names(PIG_MATCH)) {
    ie[species == "human" & pair %in% HUMAN_MATCH[[g]], matched_group := g]
    ie[species == "pig"   & pair %in% PIG_MATCH[[g]],   matched_group := g]
}

# strata
ie[, tss_bin := cut(abs_tss_dist, breaks = TSS_BREAKS, labels = TSS_LABELS,
                    include.lowest = TRUE, right = FALSE)]
ie[, maf_bin_fine   := cut(maf, breaks = MAF_BREAKS_FINE,   include.lowest = TRUE)]
ie[, maf_bin_coarse := cut(maf, breaks = MAF_BREAKS_COARSE, include.lowest = TRUE)]
ie[species == "human" & maf < MAF_FLOOR,
   `:=`(maf_bin_fine = NA, maf_bin_coarse = NA)]

# ---- assertions on the assembled table --------------------------------------

stopifnot(uniqueN(ie$species) == 2L)
if (uniqueN(ie[species == "human", pair]) != N_HUMAN_FILES)
    stop("human pair count changed", call. = FALSE)
if (uniqueN(ie[species == "pig", pair]) != N_PIG_FILES)
    stop("pig pair count changed", call. = FALSE)
for (g in names(PIG_MATCH)) {
    for (sp in c("human", "pig")) {
        if (!any(ie$species == sp & ie$matched_group == g, na.rm = TRUE))
            stop(sprintf("matched group '%s' has no %s pairs", g, sp), call. = FALSE)
    }
}
if (min(ie[species == "pig", maf]) < MAF_FLOOR - 1e-9)
    stop("pig MAF below the stated 0.10 floor -- harmonization assumption broken",
         call. = FALSE)

# =============================================================================
# ANALYSIS SETS
# =============================================================================
# Every cross-species number uses protein-coding genes at MAF >= 0.10.
# The human 0.05-0.10 bin is kept separately as a human-only stratum.

pc      <- ie[protein_coding == TRUE]
sig     <- pc[is_ieqtl == TRUE]
sig_h10 <- sig[maf >= MAF_FLOOR]                       # harmonized, both species
sig_hlo <- sig[species == "human" & maf < MAF_FLOOR]   # human-only low-MAF bin

# ---- Wilson binomial CI -----------------------------------------------------

wilson <- function(k, n, conf = 0.95) {
    z  <- qnorm(1 - (1 - conf) / 2)
    p  <- k / n
    d  <- 1 + z^2 / n
    ctr <- (p + z^2 / (2 * n)) / d
    hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
    list(lo = pmax(0, ctr - hw), hi = pmin(1, ctr + hw))
}

rate_table <- function(d, comparison, stratifier, stratum_col) {
    if (nrow(d) == 0) return(NULL)
    r <- d[, .(n_ieqtl     = .N,
               n_bulk_sig  = sum(bulk_sig),
               n_bulk_emt  = sum(bulk_emt),
               n_bulk_nom  = sum(bulk_nom),
               n_bulk_pair = sum(bulk_pair),
               median_z_g  = median(z_g),
               n_samples   = round(mean(n_samples)),
               n_pairs     = uniqueN(pair)),
           by = c("species", stratum_col)]
    setnames(r, stratum_col, "stratum")
    r[, rate := n_bulk_sig / n_ieqtl]
    r[, rate_emt := n_bulk_emt / n_ieqtl]
    r[, rate_nom := n_bulk_nom / n_ieqtl]
    r[, rate_pair := n_bulk_pair / n_ieqtl]
    ci <- wilson(r$n_bulk_sig, r$n_ieqtl)
    r[, `:=`(ci_lo = ci$lo, ci_hi = ci$hi,
             comparison = comparison, stratifier = stratifier)]
    # Preserve the factor's own ordering. Strata are binned distances and
    # frequencies, so alphabetical order is wrong ("10-50kb" would precede
    # "<10kb"); carry an explicit index so every downstream print and plot
    # inherits the intended order.
    r[, stratum_idx := if (is.factor(stratum)) as.integer(stratum) else 1L]
    r[, stratum := as.character(stratum)]
    setcolorder(r, c("comparison", "stratifier", "species", "stratum", "stratum_idx"))
    r[order(species, stratum_idx)]
}

# ---- rank-matched sensitivity ------------------------------------------------
# Within a comparison, go equally deep into each species' interaction ranking:
# K = min over species of the mean number of significant ieGenes per pair. This
# removes "pig simply calls more hits" as an explanation without discarding the
# published thresholds from the primary analysis.

rank_matched <- function(d) {
    if (nrow(d) == 0) return(d[0])
    per_pair <- d[, .N, by = .(species, pair)][, .(m = mean(N)), by = species]
    if (nrow(per_pair) < 2) return(d[0])
    K <- floor(min(per_pair$m))
    if (K < 1) return(d[0])
    d[ie_rank <= K]
}

# =============================================================================
# RESULTS
# =============================================================================

rates <- list()

# --- overall, pooled across all tissue x cell-type pairs ---
sig_h10[, .dummy := "all"]
rates$pooled_overall <- rate_table(sig_h10, "pooled", "overall", ".dummy")
rates$pooled_tss     <- rate_table(sig_h10, "pooled", "tss_distance", "tss_bin")
rates$pooled_maf     <- rate_table(sig_h10, "pooled", "maf", "maf_bin_fine")

# --- matched cell-type pairs ---
for (g in names(PIG_MATCH)) {
    m <- sig_h10[matched_group == g]
    m[, .dummy := "all"]
    rates[[paste0(g, "_overall")]] <- rate_table(m, g, "overall", ".dummy")
    rates[[paste0(g, "_tss")]]     <- rate_table(m, g, "tss_distance", "tss_bin")
    rates[[paste0(g, "_maf")]]     <- rate_table(m, g, "maf", "maf_bin_coarse")
}

# --- rank-matched sensitivity ---
rm_pooled <- rank_matched(sig_h10)
rm_pooled[, .dummy := "all"]
rates$rankmatched_overall <- rate_table(rm_pooled, "pooled_rankmatched", "overall", ".dummy")
rates$rankmatched_tss     <- rate_table(rm_pooled, "pooled_rankmatched", "tss_distance", "tss_bin")
rates$rankmatched_maf     <- rate_table(rm_pooled, "pooled_rankmatched", "maf", "maf_bin_fine")

for (g in names(PIG_MATCH)) {
    m <- rank_matched(sig_h10[matched_group == g])
    if (nrow(m) == 0) next
    m[, .dummy := "all"]
    rates[[paste0(g, "_rankmatched")]] <-
        rate_table(m, paste0(g, "_rankmatched"), "overall", ".dummy")
}

# --- human-only low-MAF stratum ---
if (nrow(sig_hlo) > 0) {
    sig_hlo[, .dummy := "0.05-0.10 (human only)"]
    rates$human_lowmaf <- rate_table(sig_hlo, "human_lowmaf", "maf", ".dummy")
}

rate_dt <- rbindlist(Filter(Negate(is.null), rates), use.names = TRUE)

# ---- shared-cutoff sweep ----------------------------------------------------
# A single uncorrected cutoff applied identically to every gene and both
# species. Sweeping it answers the question a single number cannot: is the
# species gap a property of the data, or of where the threshold was put? If the
# two curves were merely shifted versions of each other, one cutoff would make
# them agree.

SWEEP <- c(0.5, 0.1, 0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10)
sweep_dt <- rbindlist(lapply(SWEEP, function(a)
    sig_h10[, .(alpha = a, n_ieqtl = .N, n_hit = sum(pval_g < a),
                rate = mean(pval_g < a)), by = species]))
sweep_w <- dcast(sweep_dt, alpha ~ species, value.var = c("n_hit", "rate"))
sweep_w[, ratio := rate_human / rate_pig]

# ---- Fisher exact, human vs pig, within each stratum ------------------------

fisher_by_stratum <- function(r) {
    keys <- unique(r[, .(comparison, stratifier, stratum)])
    res  <- vector("list", nrow(keys))
    for (i in seq_len(nrow(keys))) {
        sub <- r[comparison == keys$comparison[i] &
                 stratifier == keys$stratifier[i] &
                 stratum    == keys$stratum[i]]
        if (uniqueN(sub$species) < 2) next
        h <- sub[species == "human"]; p <- sub[species == "pig"]
        tab <- matrix(c(h$n_bulk_sig, h$n_ieqtl - h$n_bulk_sig,
                        p$n_bulk_sig, p$n_ieqtl - p$n_bulk_sig), nrow = 2)
        ft <- fisher.test(tab)
        res[[i]] <- data.table(comparison = keys$comparison[i],
                               stratifier = keys$stratifier[i],
                               stratum    = keys$stratum[i],
                               human_rate = h$rate, pig_rate = p$rate,
                               odds_ratio = unname(ft$estimate),
                               p_fisher   = ft$p.value)
    }
    rbindlist(Filter(Negate(is.null), res))
}
fisher_dt <- fisher_by_stratum(rate_dt)

# ---- the models that actually answer Q2 and Q3 ------------------------------
# Q2/Q3 ask whether the species DIFFERENCE changes across a stratifier. That is
# a species x stratifier interaction, so the interaction coefficient is the
# answer -- the binned rates are the descriptive display of the same thing.

fit_models <- function(d, label, response = "bulk_sig") {
    d <- copy(d)
    d[, species := factor(species, levels = c("human", "pig"))]
    d[, .y := get(response)]
    out <- list()
    spec <- list(
        tss          = .y ~ species * log10(abs_tss_dist),
        tss_adj      = .y ~ species * log10(abs_tss_dist) + log10(n_samples),
        maf          = .y ~ species * maf,
        maf_adj      = .y ~ species * maf + log10(n_samples)
    )
    # The adjusted models add log10(n_samples) to separate the species effect
    # from study size. Inside a single matched pair each species contributes
    # only one sample size, so n_samples is then perfectly collinear with
    # species and the adjustment is not identifiable -- glm still returns a
    # fit, but the extra coefficient is aliased and its CI comes back NA.
    # Drop those models rather than print a meaningless row.
    identifiable <- d[, uniqueN(n_samples), by = species][, all(V1 > 1)]
    if (!identifiable) spec <- spec[!grepl("_adj$", names(spec))]

    for (nm in names(spec)) {
        m  <- glm(spec[[nm]], data = d, family = binomial)
        cf <- coef(summary(m))
        ix <- grep(":", rownames(cf))
        if (length(ix) == 0) next
        ci <- suppressMessages(confint.default(m))
        out[[nm]] <- data.table(
            comparison = label, rule = response, model = nm,
            term       = rownames(cf)[ix],
            estimate   = cf[ix, 1], se = cf[ix, 2],
            ci_lo      = ci[ix, 1], ci_hi = ci[ix, 2],
            p_value    = cf[ix, 4],
            n          = nrow(d))
    }
    rbindlist(out)
}

# Fit under both bulk rules: the gene-specific permutation threshold, and the
# single shared uncorrected cutoff. The second applies the same number to every
# gene in both species -- no per-gene correction, and no rescaling of p-values
# to a common sample size -- so any species difference under it is exactly the
# raw difference in how significant the bulk effects are.
model_dt <- rbindlist(c(
    lapply(c("bulk_sig", "bulk_nom"), function(rl) rbindlist(list(
        fit_models(sig_h10, "pooled", rl),
        fit_models(sig_h10[matched_group == "lung_epithelium"], "lung_epithelium", rl),
        fit_models(sig_h10[matched_group == "brain_neurons"],   "brain_neurons",   rl)
    ), use.names = TRUE))
), use.names = TRUE)

# =============================================================================
# OUTPUTS
# =============================================================================

OUT_GENES  <- file.path(OUTDIR, "ieqtl_bulk_detectability.tsv.gz")
OUT_RATES  <- file.path(OUTDIR, "ieqtl_bulk_rates.tsv")
OUT_MODELS <- file.path(OUTDIR, "ieqtl_bulk_models.tsv")
OUT_SWEEP  <- file.path(OUTDIR, "ieqtl_bulk_cutoff_sweep.tsv")

keep <- c(SHAPE, "abs_tss_dist", "n_samples", "is_ieqtl",
          "bulk_nom_thresh", "bulk_egene", "bulk_sig", "bulk_emt", "bulk_nom",
          "bulk_pair", "z_g",
          "ie_rank", "matched_group", "tss_bin", "maf_bin_fine", "maf_bin_coarse")
fwrite(sig[, ..keep], OUT_GENES, sep = "\t", quote = FALSE,
       compress = "gzip", na = "NA")
fwrite(rate_dt,  OUT_RATES,  sep = "\t", quote = FALSE, na = "NA")
fwrite(model_dt, OUT_MODELS, sep = "\t", quote = FALSE, na = "NA")
fwrite(sweep_dt, OUT_SWEEP,  sep = "\t", quote = FALSE, na = "NA")

# =============================================================================
# REPORT
# =============================================================================

pct <- function(x) sprintf("%.1f%%", 100 * x)

cat("\n"); hr()
cat("REGRESSION CHECKS  (values established when this analysis was designed)\n"); hr()

chk <- function(label, got, want, fmt = "%s") {
    ok <- isTRUE(all.equal(got, want, tolerance = 1e-3))
    cat(sprintf("  [%s] %-58s %s (expected %s)\n",
                if (ok) "ok" else "XX", label,
                sprintf(fmt, got), sprintf(fmt, want)))
    ok
}
all_ok <- c(
  chk("human protein-coding rows, all pairs",
      nrow(pc[species == "human"]), 702923L, "%d"),
  chk("human significant protein-coding ieGenes",
      nrow(sig[species == "human"]), 3116L, "%d"),
  chk("human pooled bulk rate, eigenMT rule, no MAF floor",
      mean(sig[species == "human", bulk_emt]), 0.6486, "%.4f"),
  chk("human pooled bulk rate, nominal, no MAF floor",
      mean(sig[species == "human", bulk_nom]), 0.7776, "%.4f"),
  chk("human background bulk rate, eigenMT rule, all PC rows",
      mean(pc[species == "human", bulk_emt]), 0.021644, "%.4f"),
  # The primary rule, at MAF >= 0.10. Established when the permutation
  # threshold replaced the ad-hoc eigenMT-on-pval_g rule.
  chk("human pooled bulk rate, PERMUTATION threshold, MAF>=0.10",
      mean(sig_h10[species == "human", bulk_sig]), 0.6982, "%.4f"),
  chk("pig pooled bulk rate, PERMUTATION threshold, MAF>=0.10",
      mean(sig_h10[species == "pig", bulk_sig]), 0.1197, "%.4f"),
  chk("median |z| of b_g, human / pig (power diagnostic)",
      median(sig_h10[species == "human", z_g]), 7.42, "%.2f"),
  chk("human Lung.Epithelial_cells ieGenes, MAF>=0.10",
      nrow(sig_h10[pair == "Lung.Epithelial_cells"]), 66L, "%d"),
  # 183, not the 195 quoted during planning: that figure predated the
  # protein-coding filter, which the pig files cannot apply on their own.
  chk("pig Lung.Alveolar_Epithelial_Type_1 ieGenes, PC, MAF>=0.10",
      nrow(sig_h10[pair == "Lung.Alveolar_Epithelial_Type_1"]), 183L, "%d"),
  chk("pig Blood.Monocytes ieGenes (all biotypes, no MAF floor)",
      nrow(ie[pair == "Blood.Monocytes" & is_ieqtl]), 3778L, "%d")
)
cat(sprintf("\n  %d/%d checks passed\n", sum(all_ok), length(all_ok)))

cat("\n"); hr()
cat("SAMPLE SIZES AND ieGENE COUNTS BY PAIR  (protein-coding, MAF >= 0.10)\n"); hr()
smry <- sig_h10[, .(n_ieqtl = .N, bulk = sum(bulk_sig), n = mean(n_samples)),
                by = .(species, pair, matched_group)][order(species, -n_ieqtl)]
smry[, rate := pct(bulk / n_ieqtl)]
print(head(smry[species == "human"], 10)); cat("  ...\n")
print(head(smry[species == "pig"], 10));   cat("  ...\n")

cat("\n"); hr()
cat("Q1  BULK-DETECTABLE RATE  (pval_g < the gene's permutation threshold)\n"); hr()
q1 <- rate_dt[stratifier == "overall"]
print(q1[, .(comparison, species, n_pairs, n_samples, n_ieqtl, n_bulk_sig,
             rate = pct(rate), ci = sprintf("[%s, %s]", pct(ci_lo), pct(ci_hi)),
             med_z = sprintf("%.2f", median_z_g))])
cat("\n  Alternative bulk rules on the same rows, for comparison:\n")
print(q1[, .(comparison, species,
             `permutation (primary)` = pct(rate),
             `+ gene is an eGene`    = pct(rate_pair),
             `eigenMT on pval_g`     = pct(rate_emt),
             `uncorrected p<0.05`    = pct(rate_nom))])
cat("\n  Fisher, human vs pig:\n")
print(fisher_dt[stratifier == "overall",
                .(comparison, human = pct(human_rate), pig = pct(pig_rate),
                  OR = sprintf("%.2f", odds_ratio),
                  p = format.pval(p_fisher, digits = 3))])

cat("\n"); hr()
cat(sprintf("Q1b  SHARED UNCORRECTED CUTOFF  (pval_g < %g, same number for every gene\n",
            BULK_ALPHA))
cat("     and both species; no per-gene correction, no rescaling to a common n)\n"); hr()
q1b <- rate_dt[stratifier == "overall"]
print(q1b[, .(comparison, species, n_ieqtl,
              `shared cutoff`   = pct(rate_nom),
              `per-gene perm`   = pct(rate),
              `ratio shared/perm` = sprintf("%.2f", rate_nom / rate))])
cat("\n  Human/pig ratio under each rule:\n")
rr <- dcast(q1b, comparison ~ species, value.var = c("rate", "rate_nom"))
print(rr[, .(comparison,
             `perm: human/pig`   = sprintf("%.1fx", rate_human / rate_pig),
             `shared: human/pig` = sprintf("%.1fx", rate_nom_human / rate_nom_pig))])

cat("\n"); hr()
cat("     Cutoff sweep, pooled. A single number applied to both species.\n"); hr()
print(sweep_w[, .(alpha = sprintf("%g", alpha),
                  human = pct(rate_human), pig = pct(rate_pig),
                  `human n` = n_hit_human, `pig n` = n_hit_pig,
                  `human/pig` = sprintf("%.1fx", ratio))])

cat("\n"); hr()
cat("Q2  RATE BY |TSS DISTANCE|\n"); hr()
q2 <- dcast(rate_dt[stratifier == "tss_distance"],
            comparison + stratum + stratum_idx ~ species,
            value.var = c("n_ieqtl", "rate"))
q2[, `:=`(human = pct(rate_human), pig = pct(rate_pig))]
print(q2[order(comparison, stratum_idx),
         .(comparison, stratum, n_human = n_ieqtl_human, human,
           n_pig = n_ieqtl_pig, pig)])
cat("\n  Same bins under the shared uncorrected cutoff:\n")
q2b <- dcast(rate_dt[stratifier == "tss_distance"],
             comparison + stratum + stratum_idx ~ species, value.var = "rate_nom")
print(q2b[order(comparison, stratum_idx),
          .(comparison, stratum, human = pct(human), pig = pct(pig))])

cat("\n  species x log10(|TSS distance|) interaction, under both rules:\n")
print(model_dt[model %in% c("tss", "tss_adj"),
               .(comparison, rule = ifelse(rule == "bulk_sig", "per-gene perm", "shared cutoff"),
                 model, estimate = sprintf("%+.3f", estimate),
                 ci = sprintf("[%+.3f, %+.3f]", ci_lo, ci_hi),
                 p = format.pval(p_value, digits = 3), n)][order(comparison, rule, model)])

cat("\n"); hr()
cat("Q3  RATE BY MAF\n"); hr()
q3 <- dcast(rate_dt[stratifier == "maf" & comparison != "human_lowmaf"],
            comparison + stratum + stratum_idx ~ species,
            value.var = c("n_ieqtl", "rate"))
q3[, `:=`(human = pct(rate_human), pig = pct(rate_pig))]
print(q3[order(comparison, stratum_idx),
         .(comparison, stratum, n_human = n_ieqtl_human, human,
           n_pig = n_ieqtl_pig, pig)])
if (!is.null(rates$human_lowmaf)) {
    cat("\n  human-only 0.05-0.10 stratum (no pig counterpart):\n")
    print(rates$human_lowmaf[, .(species, stratum, n_ieqtl, n_bulk_sig,
                                 rate = pct(rate))])
}
cat("\n  Same bins under the shared uncorrected cutoff:\n")
q3b <- dcast(rate_dt[stratifier == "maf" & comparison != "human_lowmaf"],
             comparison + stratum + stratum_idx ~ species, value.var = "rate_nom")
print(q3b[order(comparison, stratum_idx),
          .(comparison, stratum, human = pct(human), pig = pct(pig))])

cat("\n  species x MAF interaction, under both rules:\n")
print(model_dt[model %in% c("maf", "maf_adj"),
               .(comparison, rule = ifelse(rule == "bulk_sig", "per-gene perm", "shared cutoff"),
                 model, estimate = sprintf("%+.3f", estimate),
                 ci = sprintf("[%+.3f, %+.3f]", ci_lo, ci_hi),
                 p = format.pval(p_value, digits = 3), n)][order(comparison, rule, model)])

cat("\n"); hr()
cat("SENSITIVITY  rank-matched (equal depth into each species' ranking)\n"); hr()
print(rate_dt[grepl("rankmatched", comparison) & stratifier == "overall",
              .(comparison, species, n_ieqtl, n_bulk_sig, rate = pct(rate),
                ci = sprintf("[%s, %s]", pct(ci_lo), pct(ci_hi)))])

cat("\n"); hr()
cat("WROTE\n"); hr()
for (p in c(OUT_GENES, OUT_RATES, OUT_MODELS, OUT_SWEEP)) {
    cat(sprintf("  %-46s %s\n", basename(p),
                format(structure(file.size(p), class = "object_size"),
                       units = "auto")))
}
cat("\n")
if (!all(all_ok)) {
    cat("  NOTE: one or more regression checks failed -- see the table above.\n\n")
}
