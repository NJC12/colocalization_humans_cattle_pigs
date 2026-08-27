#!/usr/bin/env Rscript
# ===========================================================================
# build_dashboard_data.R
#
# Extracts everything the side-by-side simulation dashboard needs into one
# JavaScript file (`dashboard_data.js`), so the dashboard is a single static
# page that opens straight from disk -- file:// blocks fetch(), which is why
# the payload is a `var` assignment and not a .json read at load time.
#
# The loaders, the filename grammar and the outcome definitions all come from
# figures_and_tables/figure2_revision2.ipynb; its helper cells are extracted
# and sourced verbatim (see NB_CELLS) rather than copied, so the dashboard and
# the report PDFs cannot drift apart.
#
# What is exported, and why that and not the plots themselves:
#
#   traits   one row per (GWAS trait, replicate, eQTL panel) -- the
#            enloc_best_rows() reduction, which is what every fastEnloc number
#            in the report is built from. It does not depend on the RCP cutoff
#            or on the GWAS size, so exporting it lets the browser recompute
#            power, false-positive rate, the MAF-bin breakdown and the
#            precision-recall sweep at ANY cutoff instead of at the two the
#            PDFs hardcode.
#   af       allele-frequency histograms and MAF-by-selection-coefficient
#            box statistics, pre-binned. These come from ~35k variants per
#            replicate, so they are aggregated here; the two MAF floors the
#            report uses (0 and 1%) are both computed, and the histogram bin
#            width is 0.01 so the 1% floor falls exactly on a bin edge.
#
# Usage:  Rscript figures_and_tables/dashboard/build_dashboard_data.R [outfile]
# ===========================================================================

suppressPackageStartupMessages({
    library(jsonlite)
    library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
HERE <- tryCatch({
    a <- commandArgs(trailingOnly = FALSE)
    normalizePath(dirname(sub('^--file=', '', a[grep('^--file=', a)])))
}, error = function(e) getwd())
REPO    <- normalizePath(file.path(HERE, '..', '..'))
NB      <- file.path(REPO, 'figures_and_tables', 'figure2_revision2.ipynb')
OUTFILE <- if (length(args) >= 1) args[1] else file.path(HERE, 'dashboard_data.js')

# ---------------------------------------------------------------------------
# Source the notebook's helper cells. Cells 2-14 are setup, the filename
# grammar, select_runs(), the readers, load_sims(), the GWAS-significance
# calculation and the outcome definitions; cell 16 onward builds PDF pages and
# is not needed here.
# ---------------------------------------------------------------------------
NB_CELLS <- c(2, 4, 6, 8, 10, 12, 14)

extract_cells <- function(nb_path, cells) {
    nb <- fromJSON(nb_path, simplifyVector = FALSE)
    src <- vapply(cells, function(i) {
        cell <- nb$cells[[i + 1L]]          # jsonlite is 1-based, nbformat 0-based
        if (!identical(cell$cell_type, 'code'))
            stop(sprintf('notebook cell %d is %s, not code', i, cell$cell_type))
        paste(unlist(cell$source), collapse = '')
    }, character(1))
    paste(src, collapse = '\n')
}

tmp_helpers <- tempfile(fileext = '.R')
writeLines(extract_cells(NB, NB_CELLS), tmp_helpers)
suppressWarnings(suppressMessages(source(tmp_helpers)))
cat(sprintf('sourced notebook cells %s from %s\n',
            paste(NB_CELLS, collapse = ','), basename(NB)))

# ---------------------------------------------------------------------------
# find_file() override: tolerate the duplicated *_vars_* files.
#
# The two psamp arms relabelled cmaf0 (10x and 30x) carry BOTH namings of every
# variant table -- `..._maf_0.01.tsv` and `..._maf_0.tsv` -- because the
# relabel copied rather than moved. All 88 pairs are byte-identical, so the
# ambiguity is in the names only, not the data. The notebook's find_file()
# stops on any multiple match, which would drop 22 of the 140 replicates.
#
# Identical content is resolved to one file (deterministically, the first in
# sorted order) and recorded; differing content still stops, because then the
# ambiguity is real.
# ---------------------------------------------------------------------------
.dup_resolved <- new.env(parent = emptyenv())

find_file <- function(dir, pattern) {
    hits <- sort(list.files(dir, pattern = pattern))
    if (length(hits) == 0) {
        stop(sprintf('no file matching /%s/ in %s\n  directory contains: %s',
                     pattern, dir, paste(list.files(dir), collapse = ', ')),
             call. = FALSE)
    }
    if (length(hits) > 1) {
        paths <- file.path(dir, hits)
        sums <- unname(tools::md5sum(paths))
        if (length(unique(sums)) > 1) {
            stop(sprintf('%d DIFFERING files match /%s/ in %s: %s',
                         length(hits), pattern, dir, paste(hits, collapse = ', ')),
                 call. = FALSE)
        }
        .dup_resolved[[dir]] <- c(.dup_resolved[[dir]], hits[1])
    }
    file.path(dir, hits[1])
}

# ---------------------------------------------------------------------------
# Binning for the allele-frequency panels.
#
# Width 0.01 rather than the report's 30-bins-over-the-data-range: the MAF > 1%
# view then drops exactly the first bin, so one exported histogram serves both
# floors and the floor stays a live control in the browser.
# ---------------------------------------------------------------------------
AF_BIN_W  <- 0.01
AF_BREAKS <- seq(0, 0.5, by = AF_BIN_W)
AF_NBINS  <- length(AF_BREAKS) - 1L
AF_FLOORS <- c(0, 0.01)

GTEX_CATS <- c('gtex', 'gtex_smaller')

hex <- function(nm) {
    v <- grDevices::col2rgb(nm)
    sprintf('#%02X%02X%02X', v[1], v[2], v[3])
}

# Box statistics matching geom_boxplot(): hinges are the quartiles and the
# whiskers reach the furthest observation within 1.5 IQR of the hinge.
box_stats <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NULL)
    q <- unname(quantile(x, c(0.25, 0.5, 0.75), type = 7))
    iqr <- q[3] - q[1]
    lo <- min(x[x >= q[1] - 1.5 * iqr]); hi <- max(x[x <= q[3] + 1.5 * iqr])
    list(n = length(x),
         lower = q[1], middle = q[2], upper = q[3],
         ymin = lo, ymax = hi,
         n_out_lo = sum(x < lo), n_out_hi = sum(x > hi),
         min = min(x), max = max(x))
}

# ---------------------------------------------------------------------------
# best_rows_one_run(): the enloc_best_rows() reduction for a single replicate.
#
# Written out rather than calling enloc_best_rows() so it can run per replicate
# (the notebook's version consumes a whole load_sims() bundle). The arrange /
# group_by / slice is character-for-character the same, which is what keeps the
# dashboard's numbers identical to the report's.
# ---------------------------------------------------------------------------
best_rows_one_run <- function(e, run_id) {
    e %>%
        mutate(run_id = run_id) %>%
        arrange(desc(is.na(correct)), desc(correct), desc(is.na(RCP)), desc(RCP)) %>%
        group_by(gwas, run_id) %>%
        slice(1) %>%
        ungroup() %>%
        mutate(RCP = ifelse(is.na(RCP), 0, RCP)) %>%
        filter(!is.na(maf_bin), gwas_p < 1)
}

# ---------------------------------------------------------------------------
# scan_one_run(): everything one replicate directory contributes.
# ---------------------------------------------------------------------------
scan_one_run <- function(run) {
    out <- list(run_id = run$run_id, sim = run$sim, ok = TRUE, err = NULL)

    res <- try({
        # --- variants (GWAS panel): allele-frequency panels + significance ----
        v <- read_vars(run, 'gwas') %>%
            mutate(selco_bin = cut(selco, breaks = SELCO_BREAKS,
                                   labels = SELCO_LABELS, right = FALSE))

        # vexp/(1 + vexp): the GWAS non-centrality parameter is n times this, so
        # the browser can move the study size without re-reading anything.
        causal <- v %>%
            filter(beta != 0) %>%
            transmute(trait,
                      vexp = 2 * maf * (1 - maf) * beta^2,
                      ncp_per_n = vexp / (1 + vexp)) %>%
            distinct()

        out$af <- lapply(setNames(as.character(AF_FLOORS), AF_FLOORS), function(fl) {
            vv <- v %>% filter(maf >= as.numeric(fl))
            h <- function(d) as.integer(table(cut(d$maf, breaks = AF_BREAKS,
                                                  include.lowest = TRUE)))
            list(n         = nrow(vv),
                 all       = h(vv),
                 neutral   = h(vv %>% filter(selco == 0)),
                 selected  = h(vv %>% filter(selco != 0)),
                 # MAF values per selection-coefficient bin, pooled across the
                 # replicates of a simulation before the quartiles are taken.
                 by_selco  = split(vv$maf, vv$selco_bin))
        })

        # --- fastEnloc outcomes, one reduction per eQTL panel size ------------
        out$traits <- list()
        for (panel in GTEX_CATS) {
            if (!has_panel(run, 'enloc', panel)) next
            b <- best_rows_one_run(read_enloc(run, panel), run$run_id) %>%
                left_join(causal %>% rename(gwas = trait), by = 'gwas')
            out$traits[[panel]] <- data.frame(
                rep       = run$sim,
                maf_bin   = as.integer(b$maf_bin),
                rcp       = b$RCP,
                # correct is NA when the trait has no fastEnloc signal at all
                # ("underpowered fine-mapping"); that distinction has to survive
                # the round trip, so it stays a three-valued field.
                correct   = ifelse(is.na(b$correct), NA_integer_, as.integer(b$correct)),
                ncp_per_n = b$ncp_per_n,
                stringsAsFactors = FALSE)
        }
        TRUE
    }, silent = TRUE)

    if (inherits(res, 'try-error')) {
        out$ok <- FALSE
        out$err <- as.character(res)
    }
    out
}

# ---------------------------------------------------------------------------
# Walk every loadable replicate.
# ---------------------------------------------------------------------------
manifest <- suppressWarnings(build_manifest(quiet = TRUE))
runs <- manifest %>% filter(ok) %>% arrange(arm, sim_category, replicate)
cat(sprintf('%d loadable replicate(s) in %d arm(s)\n',
            nrow(runs), length(unique(runs$arm))))

ncores <- max(1L, min(detectCores() - 1L, 8L))
t0 <- Sys.time()
scanned <- mclapply(seq_len(nrow(runs)), function(i) scan_one_run(runs[i, ]),
                    mc.cores = ncores)
cat(sprintf('read %d replicate(s) in %.1f s on %d core(s)\n',
            length(scanned), as.numeric(difftime(Sys.time(), t0, units = 'secs')), ncores))

failed <- Filter(function(x) !isTRUE(x$ok), scanned)
for (f in failed) cat(sprintf('  FAILED %s: %s', f$run_id, f$err))

names(scanned) <- vapply(scanned, function(x) x$run_id, character(1))

# ---------------------------------------------------------------------------
# Group the replicates into simulations. One simulation = one simulation
# category within one arm: the arm fixes every parameter, the category fixes
# the demography and the trait-variant model, and that pair is exactly what the
# category reports compare.
# ---------------------------------------------------------------------------
one <- function(x) { u <- unique(x[!is.na(x)]); if (length(u)) u[1] else NA }

sims <- list()
groups <- runs %>% group_by(arm, sim_category) %>% group_split()

for (g in groups) {
    gid <- paste0(g$arm[1], '::', g$sim_category[1])
    ok_rows <- vapply(g$run_id, function(r) isTRUE(scanned[[r]]$ok), logical(1))
    gg <- g[ok_rows, , drop = FALSE]
    if (!nrow(gg)) next

    parts <- scanned[gg$run_id]

    # --- fastEnloc traits, pooled over replicates -------------------------
    traits <- list()
    for (panel in GTEX_CATS) {
        d <- bind_rows(lapply(parts, function(p) p$traits[[panel]]))
        if (!nrow(d)) next
        # I() everywhere a JSON ARRAY is required: auto_unbox = TRUE turns a
        # length-1 vector into a bare scalar, which would break the browser's
        # array handling for any single-replicate, single-trait arm.
        traits[[panel]] <- list(
            n         = nrow(d),
            rep       = I(d$rep),
            maf_bin   = I(d$maf_bin),
            rcp       = I(round(d$rcp, 5)),
            # JSON has no NA; the three states become 1 / 0 / null.
            correct   = I(lapply(d$correct, function(x) if (is.na(x)) NULL else x)),
            ncp_per_n = I(signif(ifelse(is.na(d$ncp_per_n), NA_real_, d$ncp_per_n), 6)))
    }

    # --- allele frequencies, pooled over replicates -----------------------
    af <- list()
    for (fl in as.character(AF_FLOORS)) {
        pieces <- lapply(parts, function(p) p$af[[fl]])
        add <- function(field) Reduce(`+`, lapply(pieces, function(p) p[[field]]))
        pooled <- lapply(setNames(SELCO_LABELS, SELCO_LABELS), function(lbl)
            unlist(lapply(pieces, function(p) p$by_selco[[lbl]]), use.names = FALSE))
        af[[fl]] <- list(
            n        = sum(vapply(pieces, function(p) p$n, numeric(1))),
            all      = I(add('all')),
            neutral  = I(add('neutral')),
            selected = I(add('selected')),
            box      = lapply(pooled, box_stats))
    }

    cat_letter <- g$sim_category[1]
    power_sampled <- identical(as.character(one(g$causal_sampling)), 'power')

    sims[[length(sims) + 1]] <- list(
        id           = gid,
        arm          = g$arm[1],
        category     = cat_letter,
        category_label = unname(CATEGORY_LABEL[cat_letter]),
        species      = one(g$species),
        color        = hex(color_key[[paste(cat_letter, 'coloc')]] %||%
                           color_key[[one(g$species)]] %||% 'gray35'),
        # Under power-weighted sampling causal_min_maf is NOT a floor on the
        # causal pool -- it only names the stage-2 directory and gates the
        # flanking GTEx-only loci -- so the causative cutoff is reported as the
        # sampling scheme instead of as a frequency. Same rule as
        # report_param_table() in the notebook.
        causal_sampling = if (power_sampled) 'power' else 'uniform',
        causal_cut   = if (power_sampled) 'power' else as.character(one(g$causal_min_maf)),
        causal_min_maf = one(g$causal_min_maf),
        sampling_gwas_n = one(g$sampling_gwas_n),
        fm_min_maf   = one(g$fm_min_maf),
        min_maf      = one(g$min_maf),
        fm_r2        = one(g$fm_r2),
        gwas_mult    = one(g$gwas_mult),
        gtex_mult    = one(g$gtex_mult),
        region       = one(g$region),
        replicates   = nrow(gg),
        reps         = I(gg$sim),
        traits       = traits,
        af           = af)
}

cat(sprintf('%d simulation(s) (arm x category)\n', length(sims)))

payload <- list(
    meta = list(
        generated      = format(Sys.time(), '%Y-%m-%d %H:%M:%S'),
        sim_root       = SIM_ROOT,
        maf_bins       = I(levels(cut(0.1, breaks = MAF_BREAKS))),
        selco_bins     = I(SELCO_LABELS),
        af_bin_width   = AF_BIN_W,
        af_breaks      = I(AF_BREAKS),
        af_floors      = I(AF_FLOORS),
        gtex_cats      = I(GTEX_CATS),
        panel_n        = as.list(PANEL_N[GTEX_CATS]),
        gwas_sig_chisq = GWAS_SIG_CHISQ,
        gwas_sig_p     = GWAS_SIG_P,
        gwas_n_default = GWAS_N_DEFAULT,
        n_failed       = length(failed)),
    category_labels = as.list(CATEGORY_LABEL),
    colors = list(
        outcome = list(
            fp  = hex('black'),
            uc  = hex('gray40'),
            ufm = hex('gray80')),
        species = list(human = hex(color_key[['human coloc']]),
                       cattle = hex(color_key[['cattle coloc']])),
        category = setNames(lapply(names(CATEGORY_LABEL),
                                   function(k) hex(color_key[[paste(k, 'coloc')]])),
                            names(CATEGORY_LABEL))),
    sims = sims)

json <- toJSON(payload, auto_unbox = TRUE, digits = NA, null = 'null', na = 'null')
writeLines(c('// Generated by figures_and_tables/dashboard/build_dashboard_data.R',
             '// Do not edit by hand.',
             paste0('window.SIM_DATA = ', json, ';')), OUTFILE)

cat(sprintf('wrote %s (%.1f MB)\n', OUTFILE, file.info(OUTFILE)$size / 1e6))
