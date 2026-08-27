#!/usr/bin/env Rscript
#
# red_herring_histograms.R
#
# The four-panel red-herring figure, repeated three times across the page:
#
#   instance 1  "All"           every colocalizing / non-colocalizing eQTL  (unchanged)
#   instance 2  "RH weaker"     pairs where the red herring eQTL has a WEAKER |effect size|
#                               than its colocalizing partner
#   instance 3  "RH stronger"   pairs where it has a STRONGER |effect size|
#
# and within each instance:
#   panel 1  |TSS distance| of the COLOCALIZING eQTL
#   panel 2  |TSS distance| of the NON-colocalizing ("red herring") eQTLs of the same
#            gene-tissue pairs
#   panel 3  (2) minus (1)
#   panel 4  fraction of (1) and (2) within 500 bp of a promoter / of an enhancer
#
# HOW THE THREE INSTANCES RELATE, stated because they are not nested the way they look.
# "The red herring has a weaker effect than the colocalizing one" is a statement about a
# PAIR, so instances 2 and 3 are built pairwise: every (red herring, colocalizing partner)
# pair of a gene-tissue is assigned to one of them, and panels 1 and 2 then show the
# distinct sets participating in that stratum's pairs. A colocalizing set with both a
# weaker and a stronger partner therefore appears in BOTH instances -- correctly, since it
# forms pairs of both kinds. Instance 1 is left exactly as the standalone figure was:
# panel 3 there measures each red herring against the MEDIAN |TSS distance| of its
# gene-tissue's colocalizing set(s), one value per set, rather than pairwise. So instance 1
# is not the row-wise union of 2 and 3, and its panel-3 n is smaller. Keeping it unchanged
# was the explicit request; the alternative would have been to make all three pairwise.
#
# EFFECT SIZE is `lead_abs_slope` from red_herring_eqtls.tsv.gz: |slope| of the set's
# top-PIP member, joined there from the GTEx v8 / PigGTEx nominal significant-pair tables.
# A slope exists only where the variant-gene pair cleared significance, so a pair enters
# instance 2 or 3 only when BOTH sides have one; the rest stay in instance 1 and the count
# dropped is reported. Comparing |slope| between two eQTLs is only meaningful because both
# describe the same gene in the same tissue, i.e. the same expression phenotype in the same
# units -- these are never compared across genes or tissues.
#
# OTHER CHOICES, unchanged from the standalone version:
#   - "the eQTL" is the eqtl_wmedian row (credible set at its PIP-weighted median).
#   - colocalizing vs not is `same_as_coloc` (exact evidence, or within 500 bp).
#   - pig is restricted to cs_truncated == FALSE, the subset comparable to human CS95.
#   - panels 1-3 bin on log10(distance), not by transforming the axis afterwards.
#   - panel 4 computes abs(pro_dist) <= 500 from the continuous distance columns; the
#     precomputed boolean flags stop at 250 bp.
#
# USAGE  Rscript red_herring_histograms.R [output_dir]

suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(patchwork)
})

EQTL_TYPE <- "eqtl_wmedian"   # or "eqtl_median"
ELEM_BP   <- 500

script_dir <- function() {
    ca <- commandArgs(trailingOnly = FALSE)
    hit <- grep("^--file=", ca, value = TRUE)
    if (length(hit) == 0) return(normalizePath("."))
    normalizePath(dirname(sub("^--file=", "", hit[1])))
}
args   <- commandArgs(trailingOnly = TRUE)
HERE   <- script_dir()
OUTDIR <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else HERE
IN     <- file.path(HERE, "red_herring_eqtls.tsv.gz")
OUT    <- file.path(OUTDIR, "red_herring_histograms.png")

INST <- c("All", "RH weaker", "RH stronger")
COL  <- c(human = "firebrick3", pig = "royalblue3")
COL4 <- c("human Colocalizing" = "firebrick4", "human Red herring" = "firebrick1",
          "pig Colocalizing"   = "navyblue",   "pig Red herring"   = "royalblue1")

# ------------------------------------------------------------------------ data

if (!file.exists(IN)) stop("input not found: ", IN)
d <- fread(IN, showProgress = FALSE)
if (!"lead_abs_slope" %in% names(d))
    stop("lead_abs_slope missing -- rerun red_herring_eqtls.R, which now joins effect sizes")

q <- d[type == EQTL_TYPE & !is.na(abs_tss_dist)]
q <- q[species == "human" | cs_truncated == FALSE]
if (nrow(q) == 0) stop("no eQTL rows after filtering -- is EQTL_TYPE right?")

SETCOLS <- c("species", "tissue", "gene", "cs_id", "abs_tss_dist", "pro_dist", "enh_dist")
coloc_all <- q[same_as_coloc == TRUE,  ..SETCOLS]
other_all <- q[same_as_coloc == FALSE, ..SETCOLS]

# --- instance 1: exactly the standalone figure -------------------------------
ref1 <- coloc_all[, .(coloc_abs = as.numeric(stats::median(abs_tss_dist))),
                  by = .(species, tissue, gene)]
d1 <- merge(other_all[, .(species, tissue, gene, abs_tss_dist)], ref1,
            by = c("species", "tissue", "gene"))
d1[, delta := abs_tss_dist - coloc_abs]

# --- instances 2 and 3: pairwise on effect size ------------------------------
oS <- q[same_as_coloc == FALSE & !is.na(lead_abs_slope),
        .(species, tissue, gene, o_id = cs_id, o_abs = abs_tss_dist,
          o_pro = pro_dist, o_enh = enh_dist, o_slope = lead_abs_slope)]
cS <- q[same_as_coloc == TRUE  & !is.na(lead_abs_slope),
        .(species, tissue, gene, c_id = cs_id, c_abs = abs_tss_dist,
          c_pro = pro_dist, c_enh = enh_dist, c_slope = lead_abs_slope)]
pr <- merge(oS, cS, by = c("species", "tissue", "gene"), allow.cartesian = TRUE)
n_tie <- pr[o_slope == c_slope, .N]
pr <- pr[o_slope != c_slope]
pr[, instance := fifelse(o_slope < c_slope, "RH weaker", "RH stronger")]
pr[, delta := o_abs - c_abs]

# --- assemble the long frames ------------------------------------------------
p1_dat <- rbindlist(list(
    cbind(instance = "All", coloc_all[, .(species, abs_tss_dist, pro_dist, enh_dist)]),
    unique(pr[, .(instance, species, cs_id = c_id, abs_tss_dist = c_abs,
                  pro_dist = c_pro, enh_dist = c_enh)])[, .(instance, species, abs_tss_dist,
                                                            pro_dist, enh_dist)]))
p2_dat <- rbindlist(list(
    cbind(instance = "All", other_all[, .(species, abs_tss_dist, pro_dist, enh_dist)]),
    unique(pr[, .(instance, species, cs_id = o_id, abs_tss_dist = o_abs,
                  pro_dist = o_pro, enh_dist = o_enh)])[, .(instance, species, abs_tss_dist,
                                                            pro_dist, enh_dist)]))
p3_dat <- rbindlist(list(cbind(instance = "All", d1[, .(species, delta)]),
                         pr[, .(instance, species, delta)]))

for (x in list(p1_dat, p2_dat, p3_dat))
    x[, `:=`(instance = factor(instance, levels = INST),
             facet    = factor(paste(instance, species), levels =
                               as.vector(t(outer(INST, c("human", "pig"), paste)))))]
p1_dat[, lx := log10(abs_tss_dist)]
p2_dat[, lx := log10(abs_tss_dist)]
slog <- function(x) sign(x) * log10(abs(x) + 1)
p3_dat[, sx := slog(delta)]

# ---------------------------------------------------------------------- panels

BP <- c(1, 100, 1e4, 1e6)
bp_lab <- function(x) ifelse(abs(x) >= 1e6, sprintf("%g Mb", x / 1e6),
                      ifelse(abs(x) >= 1e3, sprintf("%g kb", x / 1e3), sprintf("%g bp", x)))

dens_panel <- function(dt, xcol, title, sub, xlab, brk_lab, brk_pos) {
    m <- dt[, .(v = stats::median(get(xcol))), by = .(facet, species)]
    ggplot(dt, aes(x = .data[[xcol]], y = after_stat(density), fill = species)) +
        geom_histogram(binwidth = 0.15, colour = "white", linewidth = 0.1) +
        geom_vline(data = m, aes(xintercept = v), linetype = "dashed",
                   colour = "grey20", linewidth = 0.35) +
        facet_wrap(~ facet, nrow = 1, scales = "free_y") +
        scale_fill_manual(values = COL, guide = "none") +
        scale_x_continuous(breaks = brk_pos, labels = brk_lab) +
        labs(title = title, subtitle = sub, x = xlab, y = "density") +
        theme_bw(base_size = 10) +
        theme(panel.grid.minor = element_blank(),
              strip.background = element_rect(fill = "grey92", colour = NA),
              strip.text = element_text(size = 8.5),
              plot.title = element_text(face = "bold", size = 11),
              plot.subtitle = element_text(size = 8.5, colour = "grey30"))
}

lab_n <- function(dt) paste(dt[, .N, by = facet][order(facet)][
    , sprintf("%s %s", facet, format(N, big.mark = ","))], collapse = "  |  ")

p1 <- dens_panel(p1_dat, "lx", "1. Colocalizing eQTL: distance to the TSS",
        lab_n(p1_dat), "|distance to TSS|", bp_lab(BP), log10(BP))
p2 <- dens_panel(p2_dat, "lx", "2. Red herring eQTLs (non-colocalizing, same gene-tissues)",
        lab_n(p2_dat), "|distance to TSS|", bp_lab(BP), log10(BP))
DBP <- c(-1e6, -1e3, 0, 1e3, 1e6)
p3 <- dens_panel(p3_dat, "sx", "3. Difference: red herring minus colocalizing",
        paste0("positive = the red herring sits FARTHER from the TSS.  ", lab_n(p3_dat)),
        "|TSS dist| red herring  -  |TSS dist| colocalizing",
        ifelse(DBP == 0, "0", bp_lab(DBP)), slog(DBP)) +
      geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.4)

# --- panel 4 ------------------------------------------------------------------
el <- rbind(cbind(category = "Colocalizing", p1_dat[, .(instance, species, facet, pro_dist, enh_dist)]),
            cbind(category = "Red herring",  p2_dat[, .(instance, species, facet, pro_dist, enh_dist)]))
el <- melt(el, id.vars = c("category", "instance", "species", "facet"),
           measure.vars = c("pro_dist", "enh_dist"),
           variable.name = "element", value.name = "dist")
el[, element := factor(ifelse(element == "pro_dist", "Prom", "Enh"), levels = c("Prom", "Enh"))]

wilson <- function(k, n, z = 1.96) {
    ph <- k / n; dd <- 1 + z^2 / n
    c(lo = (ph + z^2 / (2 * n) - z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / dd,
      hi = (ph + z^2 / (2 * n) + z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / dd)
}
e4 <- el[, { k <- sum(abs(dist) <= ELEM_BP); n <- .N; ci <- wilson(k, n)
             .(k = k, n = n, frac = k / n, lo = ci[["lo"]], hi = ci[["hi"]]) },
         by = .(facet, instance, species, element, category)]
e4[, fillkey := paste(species, category)]

p4 <- ggplot(e4, aes(x = element, y = frac, fill = fillkey, group = category)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68,
             colour = "grey25", linewidth = 0.15) +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.75),
                  width = 0.15, linewidth = 0.3, colour = "grey20") +
    geom_text(aes(y = hi, label = sprintf("%.0f%%", 100 * frac)),
              position = position_dodge(width = 0.75), vjust = -0.5, size = 2.5) +
    # Fixed y: the species differ ~4x on promoters, and a free scale would draw pig's
    # 10% at the same height as human's 37%.
    facet_wrap(~ facet, nrow = 1, scales = "fixed") +
    scale_fill_manual(values = COL4, breaks = names(COL4),
                      labels = tolower(sub(" ", " · ", names(COL4))),
                      guide = guide_legend(title = NULL, nrow = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.20))) +
    labs(title = sprintf("4. Within %d bp of a promoter (Prom) or enhancer (Enh)", ELEM_BP),
         subtitle = "bars are 95% Wilson intervals", x = NULL,
         y = sprintf("fraction within %d bp", ELEM_BP)) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(size = 8.5),
          legend.position = "top", legend.margin = margin(b = -4),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8.5, colour = "grey30"))

png(OUT, width = 3500, height = 2900, res = 200)
print(p1 / p2 / p3 / p4 + plot_layout(heights = c(1, 1, 1, 1.02)))
invisible(dev.off())

# ------------------------------------------------------------------- report

hr <- function() cat(strrep("-", 104), "\n")
fmt <- function(x) format(round(x), big.mark = ",")
cat("\n"); hr()
cat(sprintf("TSS-distance histograms x 3 effect-size strata  (eQTL summarized as %s)\n", EQTL_TYPE))
hr()
cat(sprintf("pairs usable for the effect-size split: %s (both sides have a slope)\n",
            format(nrow(pr), big.mark = ",")))
cat(sprintf("  red herring sets with no slope: %s of %s | colocalizing: %s of %s | exact ties dropped: %d\n",
            format(other_all[, .N] - uniqueN(oS, by = c("species","tissue","gene","o_id")), big.mark = ","),
            format(other_all[, .N], big.mark = ","),
            format(coloc_all[, .N] - uniqueN(cS, by = c("species","tissue","gene","c_id")), big.mark = ","),
            format(coloc_all[, .N], big.mark = ","), n_tie))
cat("\npanel 1 / 2 medians of |TSS distance| (bp):\n")
s12 <- rbind(cbind(panel = "1 colocalizing", p1_dat[, .(n = .N, median = stats::median(abs_tss_dist)), by = .(instance, species)]),
             cbind(panel = "2 red herring",  p2_dat[, .(n = .N, median = stats::median(abs_tss_dist)), by = .(instance, species)]))
s12[, median := fmt(median)]
print(s12[order(panel, instance, species)], class = FALSE, row.names = FALSE)
cat("\npanel 3, difference (positive = red herring farther from the TSS):\n")
s3 <- p3_dat[, .(n = .N, median = fmt(stats::median(delta)),
                 frac_positive = round(mean(delta > 0), 3),
                 p = signif(stats::wilcox.test(delta)$p.value, 3)), by = .(instance, species)]
print(s3[order(instance, species)], class = FALSE, row.names = FALSE)
cat(sprintf("\npanel 4, fraction within %d bp of an element:\n", ELEM_BP))
p4t <- dcast(e4, instance + species + element ~ category, value.var = "frac")
setnames(p4t, c("Colocalizing", "Red herring"), c("coloc", "red_herring"))
p4t[, `:=`(coloc = round(coloc, 3), red_herring = round(red_herring, 3))]
p4t[, ratio := round(coloc / red_herring, 2)]
print(p4t[order(instance, species, element)], class = FALSE, row.names = FALSE)
cat("  (ratio > 1 = the colocalizing eQTL is more often near the element)\n")
hr()
cat("wrote:\n  ", OUT, "\n\n")
