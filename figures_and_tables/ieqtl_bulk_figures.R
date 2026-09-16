#!/usr/bin/env Rscript
# =============================================================================
# ieqtl_bulk_figures.R
#
# WHY THIS FILE EXISTS
# --------------------
# Draws the figure for ieqtl_bulk_detectability.R. That script answers three
# questions; this one shows them:
#
#   Q1  Are human and pig ieQTLs detectable in bulk tissue at similar rates?
#   Q2  Does the human-vs-pig similarity change with TSS distance?
#   Q3  Does it change with MAF?
#
# "Detectable in bulk" means pval_g * tests_emt < 0.05 -- the eigenMT-corrected
# genotype main effect, which is the bulk-tissue effect at average cell
# composition. See the header of ieqtl_bulk_detectability.R for why.
#
# INPUT   ieqtl_bulk_rates.tsv  (written by ieqtl_bulk_detectability.R)
# OUTPUT  ieqtl_bulk_figures.png
#
# CONVENTIONS
#   * Colours are the project's: human = firebrick, pig = royalblue. Same family
#     as red_herring_histograms.R COL and figure2_revision2.ipynb color_key.
#   * Every panel uses the SAME fixed 0-100% y axis. These are all the same
#     quantity measured on different subsets, so free scales would make
#     genuinely different rates look alike. That is a deliberate choice.
#   * Error bars are Wilson 95% binomial intervals. The matched-pair panels have
#     very few human genes per bin (as low as 1), so their intervals are wide on
#     purpose -- do not read a point estimate there without its interval.
#
# USAGE
#   Rscript ieqtl_bulk_figures.R [output_dir]
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
})

script_dir <- function() {
    ca  <- commandArgs(trailingOnly = FALSE)
    hit <- grep("^--file=", ca, value = TRUE)
    if (length(hit) == 0) return(normalizePath("."))
    normalizePath(dirname(sub("^--file=", "", hit[1])))
}

args   <- commandArgs(trailingOnly = TRUE)
HERE   <- script_dir()
OUTDIR <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else HERE

IN_RATES <- file.path(HERE, "ieqtl_bulk_rates.tsv")
OUT_PNG  <- file.path(OUTDIR, "ieqtl_bulk_figures.png")

if (!file.exists(IN_RATES)) {
    stop("ieqtl_bulk_rates.tsv not found -- run ieqtl_bulk_detectability.R first: ",
         IN_RATES, call. = FALSE)
}

r <- fread(IN_RATES, showProgress = FALSE)

COL <- c(human = "firebrick3", pig = "royalblue3")

COMPARISON_LABEL <- c(
    pooled              = "All tissue x cell-type pairs",
    lung_epithelium     = "Lung epithelium (matched)",
    brain_neurons       = "Brain neurons (matched)",
    pooled_rankmatched  = "All pairs, rank-matched",
    lung_epithelium_rankmatched = "Lung epithelium, rank-matched",
    brain_neurons_rankmatched   = "Brain neurons, rank-matched"
)

PRIMARY <- c("pooled", "lung_epithelium", "brain_neurons")
RANKED  <- c("pooled_rankmatched", "lung_epithelium_rankmatched",
             "brain_neurons_rankmatched")

r[, species := factor(species, levels = c("human", "pig"))]

# n's for the subtitles
lab_n <- function(d) {
    s <- d[, .(n = sum(n_ieqtl)), by = species]
    paste(sprintf("%s n = %s", s$species, format(s$n, big.mark = ",", trim = TRUE)),
          collapse = "   ")
}

THEME <- theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text       = element_text(size = 8.5),
          plot.title       = element_text(face = "bold", size = 11),
          plot.subtitle    = element_text(size = 8.5, colour = "grey30"),
          legend.position  = "bottom",
          legend.title     = element_blank())

pct_axis <- scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(100 * x, "%"),
                               expand = expansion(mult = c(0.02, 0.06)))

# ---- panel 1: overall rates -------------------------------------------------

d1 <- r[stratifier == "overall" & comparison %in% c(PRIMARY, RANKED)]
d1[, set := factor(ifelse(comparison %in% PRIMARY, "published FDR < 0.05",
                          "rank-matched"),
                   levels = c("published FDR < 0.05", "rank-matched"))]
d1[, base := sub("_rankmatched$", "", comparison)]
d1[, base := factor(COMPARISON_LABEL[base],
                    levels = COMPARISON_LABEL[PRIMARY])]

p1 <- ggplot(d1, aes(species, rate, colour = species,
                     shape = set, group = interaction(species, set))) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.10,
                  position = position_dodge(width = 0.5), linewidth = 0.4) +
    geom_point(position = position_dodge(width = 0.5), size = 2.2) +
    facet_wrap(~ base, nrow = 1) +
    scale_colour_manual(values = COL) +
    scale_shape_manual(values = c(16, 1)) +
    pct_axis +
    labs(title = "1. ieQTLs that are also detectable in bulk tissue",
         subtitle = paste0("pval_g x tests_emt < 0.05, protein-coding, MAF >= 0.10.  ",
                           lab_n(r[stratifier == "overall" & comparison == "pooled"])),
         x = NULL, y = "bulk-detectable") +
    # Species is already on this panel's x axis, so its colour guide is
    # redundant here -- and because this is the only panel that also maps
    # shape, its colour guide is not identical to the others', so patchwork
    # would collect it as a second copy rather than deduplicating.
    guides(colour = "none") +
    THEME

# ---- panels 2 and 3: stratified ---------------------------------------------

strat_panel <- function(stratifier_name, comparisons, title, subtitle, xlab,
                        show_facet = TRUE) {
    d <- r[stratifier == stratifier_name & comparison %in% comparisons]
    if (nrow(d) == 0) return(NULL)
    d[, facet := factor(COMPARISON_LABEL[comparison],
                        levels = COMPARISON_LABEL[comparisons])]
    d <- d[order(facet, species, stratum_idx)]
    d[, stratum := factor(stratum, levels = unique(stratum[order(stratum_idx)]))]
    p <- ggplot(d, aes(stratum, rate, colour = species, group = species)) +
        geom_line(linewidth = 0.5) +
        geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, linewidth = 0.35) +
        geom_point(size = 1.8) +
        scale_colour_manual(values = COL) +
        pct_axis +
        labs(title = title, subtitle = subtitle, x = xlab, y = "bulk-detectable") +
        THEME +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5))
    # free_x because the matched pairs use coarser MAF bins than the pooled set;
    # a shared x would print every pooled bin as an empty slot in the matched
    # facets, which reads as missing data rather than a different binning.
    if (show_facet) p <- p + facet_wrap(~ facet, nrow = 1, scales = "free_x")
    p
}

p2 <- strat_panel(
    "tss_distance", PRIMARY,
    "2. By distance to the TSS",
    paste("Both species fall off with distance; the pooled species x log10(distance)",
          "interaction is positive (p = 1.4e-4), i.e. the gap narrows as distance grows."),
    "|distance to TSS|")

# The human 0.05-0.10 MAF bin has no pig counterpart (pig was mapped at a 0.10
# floor). It is carried as its own stratum so the lowest-frequency human
# behaviour is visible, and drawn as an open point so it is never mistaken for
# part of the cross-species comparison.
d3 <- r[stratifier == "maf" & comparison %in% PRIMARY]
lowmaf <- r[comparison == "human_lowmaf"]
if (nrow(lowmaf) > 0) {
    lowmaf <- copy(lowmaf)[, `:=`(comparison = "pooled", stratum = "[0.05,0.1]",
                                  stratum_idx = 0L)]
    d3 <- rbind(d3, lowmaf, fill = TRUE)
}
d3[, facet := factor(COMPARISON_LABEL[comparison], levels = COMPARISON_LABEL[PRIMARY])]
d3 <- d3[order(facet, species, stratum_idx)]
d3[, stratum := factor(stratum, levels = unique(stratum[order(stratum_idx)]))]
d3[, harmonized := stratum_idx > 0L]

p3 <- ggplot(d3, aes(stratum, rate, colour = species, group = interaction(species, harmonized))) +
    geom_line(linewidth = 0.5) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, linewidth = 0.35) +
    geom_point(aes(shape = harmonized), size = 1.8) +
    facet_wrap(~ facet, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = COL) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    pct_axis +
    labs(title = "3. By minor allele frequency",
         subtitle = paste("Open point is the human-only 0.05-0.10 bin: pig was mapped at a 0.10",
                          "floor, so it has no pig counterpart and is excluded from every",
                          "cross-species number."),
         x = "MAF", y = "bulk-detectable") +
    THEME +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5))

# ---- panel 4: rank-matched sensitivity --------------------------------------

p4a <- strat_panel(
    "tss_distance", "pooled_rankmatched",
    "4. Rank-matched sensitivity: TSS distance",
    paste("Equal depth into each species' interaction ranking, so the species",
          "difference is not merely that pig calls more hits."),
    "|distance to TSS|", show_facet = FALSE)
p4b <- strat_panel(
    "maf", "pooled_rankmatched",
    "5. Rank-matched sensitivity: MAF",
    "Same gene set as panel 4, stratified by frequency instead of distance.",
    "MAF", show_facet = FALSE)

p4 <- p4a + p4b + plot_layout(widths = c(1, 1))

# One guide collection, at the outermost level only. Collecting inside p4 as
# well produced a second, clipped copy of the species legend.
fig <- p1 / p2 / p3 / p4 +
    plot_layout(heights = c(1, 1, 1, 1), guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal")

png(OUT_PNG, width = 3000, height = 3400, res = 200)
print(fig)
invisible(dev.off())

cat(sprintf("wrote %s (%s)\n", OUT_PNG,
            format(structure(file.size(OUT_PNG), class = "object_size"), units = "auto")))
