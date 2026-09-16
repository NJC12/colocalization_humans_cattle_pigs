#!/usr/bin/env Rscript
# =============================================================================
# ieqtl_effect_size_tss.R
#
# WHY THIS FILE EXISTS
# --------------------
# ieqtl_bulk_detectability.R showed that the RATE at which an ieQTL is also
# detectable in bulk tissue collapses with distance to the TSS, in both species
# (human 82.7% -> 5.0%, pig 28.0% -> 5.1%). This script asks the obvious
# follow-up: is that because the cell-type-specific effect itself gets weaker
# with distance?
#
#   For each species, how does |b_gi| relate to |TSS distance|, among
#   significant sites only?
#
# It is not. |b_gi| is close to flat with distance in both species, while the
# bulk effect |b_g| at the very same sites falls off a cliff. Distal ieQTLs are
# not weaker cell-type modifiers -- they are just no longer visible in bulk.
#
# WHAT b_gi IS
#   The model is expression ~ g + i + g:i + covariates. b_gi is the interaction
#   coefficient: how much the per-allele genotype effect changes per unit of the
#   (inverse-normal transformed) cell-type score. It is the cell-type-specific
#   effect, and the quantity pval_gi tests. b_g, by contrast, is the effect at
#   i = 0, i.e. at average cell composition -- the bulk effect.
#
# THE SELECTION PROBLEM, WHICH LIMITS HOW FAR THIS CAN BE PUSHED
# --------------------------------------------------------------
# These are top-association-per-gene files restricted to significant sites, so
# every |b_gi| here is winner's-cursed: median |b_gi|/SE is 5.9 (human) and 5.0
# (pig), which is deep in the selected regime. Two consequences:
#
#   * The LEVEL of |b_gi| is inflated, and inflated more at smaller sample size.
#     Pig's median |b_gi| (0.230) exceeding human's (0.185) is therefore not
#     evidence that pig ieQTLs are stronger. Do not compare levels ACROSS
#     species; only trends WITHIN a species.
#   * The spread is compressed against the detection floor, so a true
#     relationship would look flatter here than it is.
#
# The trend across distance bins is still readable because the floor itself does
# not move with distance: median b_gi_se is 0.029-0.032 across every human bin
# and 0.045-0.046 across every pig bin. The script re-checks that at runtime.
#
# ONE CONFOUNDER WORTH KNOWING
# ---------------------------
# tests_emt (the eigenMT burden, which sets each gene's significance threshold)
# correlates weakly with distance: human rho = -0.109, pig rho = +0.035. In
# human, genes whose top hit is distant face a slightly LOWER burden, hence an
# easier threshold, hence a slightly smaller detectable |b_gi| -- which pushes in
# the same direction as the observed human decline. The human slope should be
# read as a mild upper bound on the true decline, not a point estimate.
#
# INPUT   ieqtl_bulk_detectability.tsv.gz  (from ieqtl_bulk_detectability.R)
# OUTPUT  ieqtl_effect_size_tss.tsv        per-bin medians and slopes
#         ieqtl_effect_size_tss.png        two-panel figure
#
# USAGE
#   Rscript ieqtl_effect_size_tss.R [output_dir]
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
})

options(width = 190)
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

IN_GENES <- file.path(HERE, "ieqtl_bulk_detectability.tsv.gz")
OUT_TSV  <- file.path(OUTDIR, "ieqtl_effect_size_tss.tsv")
OUT_PNG  <- file.path(OUTDIR, "ieqtl_effect_size_tss.png")

if (!file.exists(IN_GENES))
    stop("run ieqtl_bulk_detectability.R first -- missing ", IN_GENES, call. = FALSE)

MAF_FLOOR  <- 0.10      # harmonized frequency range (pig was mapped at 0.10)
TSS_LEVELS <- c("<10kb", "10-50kb", "50-100kb", "100-250kb", "250-500kb", "500kb-1Mb")
COL        <- c(human = "firebrick3", pig = "royalblue3")

d <- fread(cmd = sprintf("gzcat %s", shQuote(IN_GENES)), showProgress = FALSE)
d <- d[maf >= MAF_FLOOR]
d[, `:=`(abs_bgi = abs(b_gi), abs_bg = abs(b_g),
         species = factor(species, levels = c("human", "pig")),
         tss_bin = factor(tss_bin, levels = TSS_LEVELS))]

stopifnot(nrow(d) > 0, !anyNA(d$tss_bin), all(d$is_ieqtl))

# ---- per-bin summary --------------------------------------------------------

bins <- d[, .(n          = .N,
              med_bgi    = median(abs_bgi),
              q1_bgi     = quantile(abs_bgi, .25),
              q3_bgi     = quantile(abs_bgi, .75),
              med_bg     = median(abs_bg),
              q1_bg      = quantile(abs_bg, .25),
              q3_bg      = quantile(abs_bg, .75),
              med_se_gi  = median(b_gi_se),
              med_maf    = median(maf)),
          by = .(species, tss_bin)][order(species, tss_bin)]
bins[, `:=`(rel_bgi = med_bgi / med_bgi[1], rel_bg = med_bg / med_bg[1]), by = species]

# ---- slopes -----------------------------------------------------------------
# log-log, so the slope is an elasticity: proportional change in |effect| per
# tenfold increase in distance. Adjusted models carry MAF and sample size
# because both set the detection floor, and MAF itself declines with distance
# among human hits (0.323 -> 0.230), which masks part of the trend.

slope_of <- function(x, form) {
    m  <- lm(form, data = x)
    cf <- coef(summary(m)); ci <- confint(m)
    k  <- grep("log10\\(abs_tss_dist\\)", rownames(cf))[1]
    data.table(term = deparse(form[[2]]), model = deparse(form[[3]]),
               slope = cf[k, 1], ci_lo = ci[k, 1], ci_hi = ci[k, 2],
               p_value = cf[k, 4], n = nrow(x))
}

forms <- list(
    log(abs_bgi) ~ log10(abs_tss_dist),
    log(abs_bgi) ~ log10(abs_tss_dist) + maf + log10(n_samples),
    log(abs_bg)  ~ log10(abs_tss_dist),
    log(abs_bg)  ~ log10(abs_tss_dist) + maf + log10(n_samples)
)
slopes <- rbindlist(lapply(levels(d$species), function(s)
    rbindlist(lapply(forms, function(f) cbind(species = s, slope_of(d[species == s], f))))))

# Spearman, distribution-free, as a check that the log-log fit is not doing the work
spear <- rbindlist(lapply(levels(d$species), function(s) {
    x  <- d[species == s]
    ct <- suppressWarnings(cor.test(x$abs_bgi, x$abs_tss_dist, method = "spearman"))
    data.table(species = s, rho = unname(ct$estimate), p_value = ct$p.value, n = nrow(x))
}))

# Detection floor: if b_gi_se trended with distance, a flat |b_gi| could be an
# artifact of the threshold moving rather than a real result.
floor_chk <- rbindlist(lapply(levels(d$species), function(s) {
    m <- lm(log(b_gi_se) ~ log10(abs_tss_dist), data = d[species == s])
    cf <- coef(summary(m))
    data.table(species = s, se_slope = cf[2, 1], p_value = cf[2, 4])
}))

fwrite(bins, OUT_TSV, sep = "\t", quote = FALSE, na = "NA")

# ---- figure -----------------------------------------------------------------

THEME <- theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          plot.title       = element_text(face = "bold", size = 11),
          plot.subtitle    = element_text(size = 8.5, colour = "grey30"),
          legend.position  = "bottom", legend.title = element_blank(),
          axis.text.x      = element_text(angle = 45, hjust = 1, size = 7.5))

panel <- function(med, q1, q3, title, subtitle, ylab) {
    ggplot(bins, aes(tss_bin, .data[[med]], colour = species, group = species)) +
        geom_ribbon(aes(ymin = .data[[q1]], ymax = .data[[q3]], fill = species),
                    alpha = 0.13, colour = NA) +
        geom_line(linewidth = 0.6) +
        geom_point(size = 1.9) +
        scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
        labs(title = title, subtitle = subtitle, x = "|distance to TSS|", y = ylab) +
        THEME
}

p1 <- panel("med_bgi", "q1_bgi", "q3_bgi",
            "1. Cell-type-specific effect |b_gi|",
            "Near-flat in both species. Band is the IQR. Levels are winner's-cursed and not comparable across species.",
            "median |b_gi|")
# Compute the steepness ratio rather than hardcoding it, so the subtitle cannot
# drift from the numbers. Uses the adjusted slopes, which are the ones the
# report leads with.
adj   <- slopes[grepl("maf", model)]
ratio <- adj[term == "log(abs_bg)", .(species, bg = slope)][
             adj[term == "log(abs_bgi)", .(species, bgi = slope)], on = "species"][
             , .(species, r = bg / bgi)]
ratio_txt <- paste(sprintf("%.0fx in %s", ratio$r, ratio$species), collapse = ", ")

p2 <- panel("med_bg", "q1_bg", "q3_bg",
            "2. Bulk effect |b_g|, same sites",
            sprintf("Collapses with distance -- steeper than |b_gi| by %s (adjusted log-log slopes).",
                    ratio_txt),
            "median |b_g|")

fig <- p1 + p2 + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

png(OUT_PNG, width = 2400, height = 1150, res = 200)
print(fig)
invisible(dev.off())

# ---- report -----------------------------------------------------------------

cat("\n"); hr()
cat("|b_gi| AND |b_g| BY |TSS DISTANCE|   significant protein-coding ieGenes, MAF >= 0.10\n"); hr()
print(bins[, .(species, tss_bin, n,
               `med |b_gi|` = sprintf("%.3f", med_bgi),
               `rel`        = sprintf("%3.0f%%", 100 * rel_bgi),
               `med |b_g|`  = sprintf("%.3f", med_bg),
               `rel `       = sprintf("%3.0f%%", 100 * rel_bg),
               `med SE_gi`  = sprintf("%.3f", med_se_gi),
               `med MAF`    = sprintf("%.3f", med_maf))])

cat("\n"); hr(); cat("SPEARMAN  |b_gi| vs |TSS distance|\n"); hr()
print(spear[, .(species, rho = sprintf("%+.4f", rho),
                p = format.pval(p_value, digits = 3), n)])

cat("\n"); hr()
cat("LOG-LOG SLOPES  (proportional change in |effect| per 10x distance)\n"); hr()
print(slopes[, .(species, effect = sub("^log\\((.*)\\)$", "\\1", term),
                 adjusted = ifelse(grepl("maf", model), "yes", "no"),
                 slope = sprintf("%+.4f", slope),
                 ci    = sprintf("[%+.4f, %+.4f]", ci_lo, ci_hi),
                 p     = format.pval(p_value, digits = 3))])

cat("\n"); hr(); cat("DETECTION-FLOOR CHECK  log(b_gi_se) ~ log10(distance)\n"); hr()
print(floor_chk[, .(species, se_slope = sprintf("%+.4f", se_slope),
                    p = format.pval(p_value, digits = 3))])
cat("\n  A flat SE means the significance floor does not move with distance, so\n",
    " the |b_gi| trend is not an artifact of the threshold shifting.\n", sep = "")

cat("\n"); hr(); cat("WROTE\n"); hr()
for (p in c(OUT_TSV, OUT_PNG))
    cat(sprintf("  %-40s %s\n", basename(p),
                format(structure(file.size(p), class = "object_size"), units = "auto")))
cat("\n")
