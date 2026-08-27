#!/usr/bin/env Rscript
#
# estimating_enloc_test_density.R
#
# How many fastEnloc tests is each GWAS locus actually subject to?
#
# WHY THIS FILE EXISTS. A colocalization "hit" is only interpretable against the
# number of comparisons that produced it. The simulations in this repo put 100
# GTEx traits and 50 GWAS causal loci in a single 2 Mb region, so every simulated
# GWAS locus is compared against ~100 eQTL signals. This script measures the real
# figure, so that number can be judged rather than assumed.
#
# THE DATA. Barbeira et al. 2021 (the GTEx v8 GWAS companion) fastEnloc results:
#   <ENLOC_DIR>/<TRAIT>__PM__<TISSUE>.enloc.rst.gz    87 traits x 49 tissues
#   gwas_locus  molecular_qtl_trait  locus_gwas_pip  locus_rcp  lead_coloc_SNP  lead_snp_rcp
#
# WHY NOT THE .enloc.sig.out FILES. The other human fastEnloc results on disk
# (Hukku 2021, enloc/henloc/hukku_2021/*/*/*.enloc.sig.out) cannot answer this
# question at all: their Signal column is a bare eQTL signal ("ENSG00000224969:2")
# with no "(@)<gwas_locus>" suffix, so RCP is reported per eQTL signal against the
# GWAS as a whole and there is no GWAS-locus grouping to count within. Verified:
# 0 of the 49 files on disk contain "(@)". Only the pig results
# (enloc/penloc/<Tissue>/*.enloc.sig.out) carry GWAS locus ids; this script is
# human-only by design and pig is deliberately out of scope.
#
# UNITS, stated plainly because they are not the same as the pig files':
#   - a GWAS locus is Barbeira's `gwas_locus` ("region1450"), an LD-independent
#     region, and it is TRAIT-SPECIFIC: region1450 is a different observation for
#     Eosinophil counts than for Height.
#   - one "test" is one row = one (gwas_locus, gene) pair. `molecular_qtl_trait`
#     is a GENE, not a DAP-G signal cluster, so multiple independent eQTL signals
#     in the same gene are collapsed. Rows are unique within a file (verified:
#     8,187 rows / 8,187 distinct pairs), so "count rows" and "count distinct
#     genes" are the same number.
#   - one observation is one (trait, tissue, gwas_locus) cell. Tissues are NOT
#     summed; a locus contributes up to 49 separate observations.
#
# THE TWO FLOORS. `locus_rcp` is printed to 3 decimals, so the pig files' RCP >
# 1e-4 floor is not representable here -- 79% of rows print as exactly 0.000. So
# two counts are reported per cell:
#   n_tests_all     every row in the file (floor 0; the complete denominator)
#   n_tests_rcp001  rows with locus_rcp > 0.001 (the nearest representable analogue)
# Under the 0.001 floor ~14% of cells lose every test. Those cells are summarized
# BOTH ways -- zeros kept (same locus universe as the unfiltered floor, so the two
# rows are comparable) and zeros dropped -- because the choice moves the mean a
# long way and should be visible rather than baked in.
#
# THE COLOCALIZING / NEVER-COLOCALIZING SPLIT.
#   A (trait, gwas_locus) COLOCALIZES if locus_rcp > 0.5 in at least one tissue.
#   For the colocalizing pool, only the tissues that actually reached RCP > 0.5
#   contribute observations -- not all 49. For the never-colocalizing pool, every
#   tissue contributes, since none of them cleared the cut.
#
# TRAIT GROUPS. All 87 traits are pooled for the headline, and additionally split
# into the 14 Astle et al. 2016 blood-cell traits (heavily correlated with each
# other, sharing many of the same regions) and the 73 others, so it is visible
# whether the Astle family drives the pooled result.
#
# USAGE
#   Rscript estimating_enloc_test_density.R [output_dir]
# Default output_dir is the directory holding this script. Runtime ~1.5 min.

suppressPackageStartupMessages(library(data.table))

# The summary table is 27 rows x 11 columns; without this it wraps and the max
# column is torn off onto its own block.
options(width = 200)

# ---------------------------------------------------------------- configuration

ENLOC_DIR <- "~/eqtl_selection/data/enloc/henloc/barbeira_2021/results_eur"

RCP_FLOOR <- 0.001   # secondary floor; strict >, matching "RCP > 0.001"
COLOC_CUT <- 0.5     # strict >, matching "RCP > 0.5"

ASTLE_PREFIX <- "Astle_et_al_2016"

# Resolve the script's own directory the way build_dashboard_data.R does, so the
# outputs land next to the script no matter where Rscript was invoked from.
script_dir <- function() {
    ca <- commandArgs(trailingOnly = FALSE)
    hit <- grep("^--file=", ca, value = TRUE)
    if (length(hit) == 0) return(normalizePath("."))
    normalizePath(dirname(sub("^--file=", "", hit[1])))
}

args <- commandArgs(trailingOnly = TRUE)
HERE <- script_dir()
OUTDIR <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else HERE

PER_LOCUS_OUT <- file.path(OUTDIR, "enloc_test_density_per_locus.tsv.gz")
SUMMARY_OUT   <- file.path(OUTDIR, "enloc_test_density_summary.tsv")

# ------------------------------------------------------------------- read + fold

# One pass over the 4,263 files. Each file is reduced to one row per gwas_locus
# BEFORE being kept, so peak memory stays at ~1.3M rows rather than the ~35M rows
# on disk. That is also why no intermediate cache is written: the full pass costs
# about 90 seconds, which is cheaper than managing a stale cache.
read_enloc_density <- function(enloc_dir) {
    dir <- path.expand(enloc_dir)
    if (!dir.exists(dir)) stop("enloc directory not found: ", dir)

    # The anchored pattern matters: results_eur/ also holds barbeira_2021_rcp_10.tsv
    # and barbeira_2021_rcp_50.tsv, which are pre-thresholded summaries and must
    # not be mistaken for per-trait-per-tissue results.
    files <- list.files(dir, pattern = "[.]enloc[.]rst[.]gz$", full.names = TRUE)
    if (length(files) == 0) stop("no *.enloc.rst.gz found under: ", dir)
    message(sprintf("Reading %d enloc result files from %s", length(files), dir))

    stems <- sub("[.]enloc[.]rst[.]gz$", "", basename(files))
    if (any(!grepl("__PM__", stems, fixed = TRUE)))
        stop("filenames not in <TRAIT>__PM__<TISSUE> form: ",
             paste(head(stems[!grepl("__PM__", stems, fixed = TRUE)], 3), collapse = ", "))
    traits  <- sub("__PM__.*$", "", stems)
    tissues <- sub("^.*__PM__", "", stems)

    tick <- max(1L, length(files) %/% 20L)
    out <- vector("list", length(files))
    empty <- logical(length(files))
    for (i in seq_along(files)) {
        d <- fread(files[i], showProgress = FALSE)
        # 245 of the files are header-only -- all 49 tissues of 5 traits that
        # contributed no GWAS regions to the enloc run at all. They must be
        # skipped explicitly, not just left to fall out of rbindlist: data.table
        # evaluates the j expression once on an empty table to infer the result
        # structure, so max() on a zero-length vector fires "no non-missing
        # arguments to max; returning -Inf" twice per file (490 warnings) even
        # though the returned table is correctly empty.
        if (nrow(d) == 0L) { empty[i] <- TRUE; next }
        # Fold to one row per GWAS locus. n_tests_all counts every evaluated
        # (locus, gene) pair; max_locus_rcp is what the 0.5 classification uses.
        out[[i]] <- d[, .(n_tests_all    = .N,
                          n_tests_rcp001 = sum(locus_rcp > RCP_FLOOR),
                          max_locus_rcp  = max(locus_rcp),
                          locus_gwas_pip = max(locus_gwas_pip)),
                      by = gwas_locus]
        out[[i]][, `:=`(trait = traits[i], tissue = tissues[i])]
        if (i %% tick == 0) message(sprintf("  ... %d / %d files", i, length(files)))
    }

    cells <- rbindlist(out)
    setcolorder(cells, c("trait", "tissue", "gwas_locus", "n_tests_all",
                         "n_tests_rcp001", "max_locus_rcp", "locus_gwas_pip"))
    # Carried out for the report so the drop from 87 named traits to 82 with data
    # is stated rather than inferred.
    setattr(cells, "n_files", length(files))
    setattr(cells, "n_empty_files", sum(empty))
    setattr(cells, "empty_traits", sort(unique(traits[empty])))
    setattr(cells, "all_traits", sort(unique(traits)))
    cells[]
}

# ---------------------------------------------------------------- classification

# locus_colocalizes is a property of the (trait, gwas_locus) pair pooled over
# tissues; cell_colocalizes is a property of the individual (trait, tissue,
# gwas_locus) cell. The colocalizing pool uses the latter, which is what
# restricts it to "the tissue or tissues with RCP > 0.5".
add_coloc_flags <- function(cells) {
    cells[, cell_colocalizes := max_locus_rcp > COLOC_CUT]
    cells[, locus_colocalizes := any(cell_colocalizes), by = .(trait, gwas_locus)]
    cells[, pool := fifelse(!locus_colocalizes, "never_colocalizing",
                    fifelse(cell_colocalizes, "colocalizing", "colocalizing_other_tissue"))]
    cells[]
}

# ---------------------------------------------------------------------- summarize

describe <- function(x) {
    q <- unname(quantile(x, c(0.25, 0.5, 0.75), type = 7))
    data.table(n_obs = length(x), mean = mean(x), median = q[2],
               q25 = q[1], q75 = q[3], min = min(x), max = max(x))
}

# One summary row per (floor x zeros-policy x pool x trait group). `subsets` is
# built explicitly rather than by a cross-join so that the impossible cells --
# a zeros policy for the unfiltered floor, where no cell can be zero -- never
# appear and cannot be misread as a result.
summarize_density <- function(cells) {
    trait_groups <- list(
        all_traits = function(dt) dt,
        astle      = function(dt) dt[startsWith(trait, ASTLE_PREFIX)],
        non_astle  = function(dt) dt[!startsWith(trait, ASTLE_PREFIX)]
    )
    pools <- list(
        all_loci           = function(dt) dt,
        colocalizing       = function(dt) dt[pool == "colocalizing"],
        never_colocalizing = function(dt) dt[pool == "never_colocalizing"]
    )
    floors <- list(
        list(label = "all_rows",      zeros = "n/a",     col = "n_tests_all",    drop0 = FALSE),
        list(label = "rcp_gt_0.001",  zeros = "kept",    col = "n_tests_rcp001", drop0 = FALSE),
        list(label = "rcp_gt_0.001",  zeros = "dropped", col = "n_tests_rcp001", drop0 = TRUE)
    )

    rows <- list()
    for (fl in floors) for (pn in names(pools)) for (tg in names(trait_groups)) {
        dt <- pools[[pn]](trait_groups[[tg]](cells))
        if (fl$drop0) dt <- dt[get(fl$col) > 0L]
        if (nrow(dt) == 0) next
        rows[[length(rows) + 1L]] <- cbind(
            data.table(floor = fl$label, zeros = fl$zeros, pool = pn, traits = tg,
                       n_trait_loci = uniqueN(dt, by = c("trait", "gwas_locus"))),
            describe(dt[[fl$col]]))
    }
    rbindlist(rows)
}

# --------------------------------------------------------------------------- run

cells <- add_coloc_flags(read_enloc_density(ENLOC_DIR))
summary_dt <- summarize_density(cells)

fwrite(cells, PER_LOCUS_OUT, sep = "\t", quote = FALSE, compress = "gzip")
fwrite(summary_dt, SUMMARY_OUT, sep = "\t", quote = FALSE)

# ----------------------------------------------------------------- terminal report

hr <- function() cat(strrep("-", 96), "\n")

cat("\n")
hr()
cat("fastEnloc TEST DENSITY per GWAS locus -- Barbeira 2021, GTEx v8 EUR\n")
hr()
cat(sprintf("files                     : %d  (%d header-only, excluded)\n",
            attr(cells, "n_files"), attr(cells, "n_empty_files")))
cat(sprintf("traits                    : %d of %d named  (%d Astle, %d non-Astle)\n",
            uniqueN(cells$trait), length(attr(cells, "all_traits")),
            uniqueN(cells[startsWith(trait, ASTLE_PREFIX)]$trait),
            uniqueN(cells[!startsWith(trait, ASTLE_PREFIX)]$trait)))
if (length(attr(cells, "empty_traits")) > 0)
    cat(sprintf("  no enloc regions anywhere: %s\n",
                paste(attr(cells, "empty_traits"), collapse = ", ")))
cat(sprintf("tissues                   : %d\n", uniqueN(cells$tissue)))
cat(sprintf("distinct regions          : %d\n", uniqueN(cells$gwas_locus)))
cat(sprintf("(trait, locus) pairs      : %d\n", uniqueN(cells, by = c("trait", "gwas_locus"))))
cat(sprintf("(trait, tissue, locus)    : %d observations\n", nrow(cells)))
cat(sprintf("  of these, colocalizing  : %d cells from %d (trait, locus) pairs\n",
            nrow(cells[pool == "colocalizing"]),
            uniqueN(cells[pool == "colocalizing"], by = c("trait", "gwas_locus"))))
cat(sprintf("  never colocalizing      : %d cells from %d (trait, locus) pairs\n",
            nrow(cells[pool == "never_colocalizing"]),
            uniqueN(cells[pool == "never_colocalizing"], by = c("trait", "gwas_locus"))))
cat(sprintf("  colocalizing locus, but tissue below cut (excluded): %d cells\n",
            nrow(cells[pool == "colocalizing_other_tissue"])))
hr()

pr <- summary_dt[, .(floor, zeros, pool, traits, n_obs,
                     mean = round(mean, 2), median, q25, q75, max)]
print(pr, nrows = 100, class = FALSE, row.names = FALSE)
hr()
cat("wrote:\n  ", PER_LOCUS_OUT, "\n  ", SUMMARY_OUT, "\n\n")
