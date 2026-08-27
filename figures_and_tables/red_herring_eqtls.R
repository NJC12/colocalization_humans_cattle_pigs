#!/usr/bin/env Rscript
#
# red_herring_eqtls.R
#
# Where does a colocalizing variant sit, relative to the fine-mapped eQTLs of the
# same gene and tissue?
#
# WHY THIS FILE EXISTS. A gene-tissue-trait colocalization says that an eQTL and a
# GWAS signal share a causal variant. It does not say that the shared variant is a
# typical eQTL for that gene. If colocalizing variants sit systematically farther
# from (or closer to) the TSS than the ordinary fine-mapped eQTL signals of the same
# gene-tissue pair, then the colocalizing signal is its own class of variant, and
# reading it as "the eQTL that explains the trait" is a red herring. This script
# builds the table that lets that be measured: one row per colocalization, plus one
# row per fine-mapped eQTL credible set of the same gene-tissue pair, each carrying
# a distance to that gene's TSS.
#
# THE DATA.
#   human coloc  Barbeira et al. 2021 GTEx v8 EUR, pre-thresholded to RCP > 0.5:
#                barbeira_2021/results_eur/barbeira_2021_rcp_50.tsv
#                28,196 rows; 58 traits x 49 tissues; 3,477 genes; 16,174 gene-tissues.
#   human eQTL   GTEx v8 DAP-G 95% credible sets (GTEx_v8_finemapping_DAPG.CS95.txt.gz),
#                the same fine-mapping Barbeira's enloc run consumed.
#   pig coloc    PigGTEx fastEnloc, pre-thresholded: penloc/all_tis_rcp_gt_50.tsv
#                22,817 rows; 110 traits x 21 tissues; 4,212 genes; 6,036 gene-tissues.
#   pig eQTL     the PIP_qtl column of the pig .enloc.snp.out files -- i.e. the DAP-G
#                fine-mapping fastEnloc itself ran on.
#
# WHY NOT PigGTEx SuSiE FOR THE PIG eQTLs. all_tissues.susieinf.txt with cs != -1
# yields credible sets for only 1,191 of the 6,036 pig coloc gene-tissue pairs (20%):
# Muscle 35%, Adipose 18%, Spleen 2%, and seven coloc tissues at 0%. That is not a
# tissue-name mismatch (all 21 names match); SuSiE-inf simply produces few
# purity-passing sets at PigGTEx sample sizes. Taking the eQTL signals straight out of
# the enloc .snp.out files instead gives ~100% coverage AND uses exactly the signals
# the colocalization was computed from, which mirrors the human side.
#
# THE TSS CONVENTION, which is the reason the two row types are comparable at all.
# GTEx v8 and PigGTEx both define tss_distance = variant_pos - TSS_pos, with NO strand
# flip; TSS_pos is txStart+1 for + genes and txEnd for - genes. Verified on minus-strand
# genes in both species (human ENSG00000227232/WASH7P: 666028 - 29553 = 636475, matching
# the GTEx eGenes row; pig ENSSSCG00000005383: 240512625 - 240984823 = -472198, matching
# pig_enloc_tss_pos.tsv). So this script does not join anyone's precomputed tss_dist --
# it recovers one TSS position per gene by inverting that identity, then applies
# tss_dist = pos - tss_pos identically to every row of both types in both species.
# CONSEQUENCE: tss_dist is signed by GENOME COORDINATE, NOT by strand. The sign says which
# side of the TSS the variant sits on in chromosome coordinates; it does NOT say
# upstream/downstream. For a minus-strand gene the biological reading is inverted --
# positive tss_dist is upstream (5'), while on a plus-strand gene positive is downstream.
# Verified directly against pig strand from pig_encode_genes.txt. Use the sign only to
# locate the variant relative to the TSS; `abs_tss_dist` is carried as its own column and
# is the quantity every distance statistic in this script and its report uses.
#
# THE SAME-eQTL RULE. A gene-tissue's colocalizing signal is also one of its fine-mapped
# eQTL signals, and fine-mapping routinely splits one signal across neighbouring DAP-G
# clusters. Left alone, the colocalizing eQTL therefore sits inside the very background it
# is being compared against, which pulls the contrast toward "no difference". So each set
# carries `dist_to_coloc` -- the distance from its nearest MEMBER variant to the nearest
# colocalizing position of its own gene-tissue -- and `same_as_coloc`, true when the set is
# the coloc signal by exact evidence (member-of-set for human, matching signal id for pig)
# OR lies within COLOC_SAME_BP (default 500) of it. Filter to `same_as_coloc == FALSE` for
# the genuinely additional eQTLs; the report prints the headline both ways. `dist_to_coloc`
# is member-based and so identical across a set's median and wmedian rows;
# `dist_to_coloc_pos` is the same measurement taken at each row's own summary position and
# is provided because the two disagree whenever the synthetic median drifts.
#
# EXPECT PIG DISTANCES ~10x HUMAN ONES, and do not read that as a bug in this script.
# It is in the published data: median |tss_distance| over the shipped fine-mapped eQTL
# tables is 241,664 bp for pig (pig_eqtl_tss_dist.tsv.gz, n=460,845) against 23,710 bp
# for human (finemapped_eqtl_tss_dist.tsv.gz, n=9,818,957). PigGTEx fine-mapping
# localizes far less sharply within its cis window than GTEx v8 does.
#
# PROMOTER / ENHANCER COLUMNS. Every position also gets a signed distance to the nearest
# promoter and enhancer, from the cached element BEDs in snps_elements_enrichments/
# (hpro/henh/ppro/penh.tsv). `pro_dist` reproduces `bedtools closest -D ref` exactly -- the
# convention was reverse-engineered from those tables' own pro_dist/enh_dist columns and is
# 100.0000% exact on all 350,370 rows available to check against, so no bedtools install,
# no ~/tmp/bedtools scratch dir, and no temp files are needed. `in_pro`/`in_enh` are the
# enrichment analysis's default and headline criterion (get_enr(cutoff = 0): strict
# overlap); near_*_50 and near_*_100 are its published sensitivity pair; near_*_250
# corresponds to a bed extension used only in the superseded enh_pro_enrichment.ipynb.
#
# WHY THERE ARE ALSO frac_mem_* COLUMNS. On an eqtl_* row, `pos` is a credible-set MEDIAN:
# a synthetic point that need not be any real variant and can land in a gap between
# elements. TSS distance tolerates that because it is smooth in position; element overlap
# does not, because it is a step function. So each set additionally carries the fraction of
# its actual member variants inside an element, at each cutoff. PREFER frac_mem_* FOR ANY
# REGULATORY CLAIM ABOUT A CREDIBLE SET; pro_dist on such a row describes the summary
# point only. The terminal report prints how far the two diverge. These fractions are NA on
# enloc rows, where the position is a single real variant and in_pro/in_enh already answer
# the question, and identical between a set's median and wmedian rows.
#
# THE ELEMENT DEFINITIONS ARE NOT MATCHED ACROSS SPECIES. Human promoters/enhancers are
# EpiMap 18-state TssA / EnhA1+EnhA2 thresholded at >= 41 of 833 epigenomes (~top 5%); pig
# are Pan et al. 2021 E1 / E6 present in >= 1 of 14 tissues, with no threshold at all.
# Genome coverage happens to land close (promoters 1.08% human vs 1.15% pig, enhancers
# 6.54% vs 4.79%), so the columns are broadly comparable in aggregate, but a human-vs-pig
# difference in one element class is a difference in annotation depth at least as much as
# in biology. Compare within a species by preference.
#
# WHAT A ROW IS.
#   type = "enloc"        one colocalization. Position is the colocalizing SNP:
#                         lead_coloc_SNP for human, top-SCP SNP from .snp.out for pig.
#   type = "eqtl_median"  one fine-mapped eQTL credible set of the same gene-tissue,
#                         summarized at the plain median of its member positions.
#   type = "eqtl_wmedian" the same set, summarized at the PIP-weighted median.
# eQTLs are a property of a gene-tissue pair, not of a trait, so the eqtl rows are
# emitted ONCE per (species, tissue, gene) and carry trait = NA. Pair an enloc row to
# its eQTL background by joining on (species, tissue, gene).
# sub("_.*", "", type) collapses this back to the plain enloc/eqtl distinction.
#
# COVERAGE, stated because it is not complete and should not look complete.
#   - Only 10,963 of 16,174 human coloc gene-tissue pairs (68%) have any CS95 credible
#     set, because CS95 requires signal-level PIP > 0.95. The other 5,211 pairs keep
#     their enloc row and contribute no eqtl rows. The report prints this count.
#   - gene -> TSS covers 3,468/3,477 human coloc genes and 4,175/4,212 pig. The rest
#     get tss_dist = NA and are named in the report.
#   - A pig .snp.out lists only the SNPs fastEnloc evaluated (the eQTL signal
#     intersected with the GWAS region), so a signal's credible set can be truncated
#     before reaching 95% cumulative PIP. Those rows are flagged cs_truncated, not
#     silently passed off as complete sets. This matters more than it sounds: 67,153
#     of the 104,534 pig sets (64%) are truncated, and they are weak, not merely
#     clipped -- median signal-level PIP 0.073, against 0.982 for the rest. FILTER TO
#     cs_truncated == FALSE (37,381 sets) FOR THE SUBSET COMPARABLE TO HUMAN CS95,
#     which is itself a signal-level PIP > 0.95 cut. The unfiltered pig rows are kept
#     in the file because the flag makes the choice reversible; human rows are all
#     cs_truncated = FALSE by construction.
#
# GOTCHAS THAT SILENTLY CORRUPT THE JOIN.
#   - CS95 truncates exactly two tissue names: Skin_Not_Sun_Exposed and
#     Skin_Sun_Exposed. The other 47 match GTEx/Barbeira. Asserted below.
#   - CS95 signal ids look like ENSG00000117620.14:1. Strip ONLY the version. A naive
#     sub("\\..*", "") eats the ":1" cluster index too and merges every signal of a
#     gene into one.
#   - The pig .enloc.{sig,snp}.out files are NOT tab-delimited: the header is
#     tab-separated but data rows are space-padded. Split on whitespace.
#   - 11,937 of 22,817 pig coloc rows carry a "LocNNNNNN" GWAS-locus id with no
#     coordinate. That does not matter here, because the position comes from the SCP
#     column rather than the (@) suffix -- but it IS why the cached penloc.tsv from
#     figures3and4.ipynb has only 9,646 rows. Do not cross-check against that file.
#   - tar -xO emits member contents with no filename markers, so the pig trait name is
#     recovered by counting header lines against the tar -t listing. Verified one
#     header per member (Milk: 268 members, 268 "Signal" lines) and asserted per tissue.
#
# RUNTIME ~15-25 min, dominated by the pig pass: 21 tissues of .enloc.snp.out, of which
# Muscle alone is 4.7 GB and 15.2M matching rows. Those 15.2M rows collapse to 545k
# distinct (signal, SNP) pairs, so the dedup happens in awk and only the collapsed
# result is ever handed to R. No intermediate cache is written -- managing a stale one
# across two source trees costs more than the rerun.
#
# USAGE
#   Rscript red_herring_eqtls.R [output_dir]
# Default output_dir is the directory holding this script.

suppressPackageStartupMessages(library(data.table))

options(width = 200)

# ---------------------------------------------------------------- configuration

HUMAN_COLOC   <- "~/eqtl_selection/data/enloc/henloc/barbeira_2021/results_eur/barbeira_2021_rcp_50.tsv"
HUMAN_CS95    <- "~/eqtl_selection/data/gtex/hgtex/finemapped/GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.CS95.txt.gz"
HUMAN_TSS_A   <- "~/eqtl_selection/data/gtex/hgtex/finemapped/GTEx_v8_finemapping_DAPG/finemapped_eqtl_tss_dist.tsv.gz"
HUMAN_TSS_B   <- "~/data/gtex_v8/all_eqtl_tss_dist.tsv.gz"

PIG_COLOC     <- "~/eqtl_selection/data/enloc/penloc/all_tis_rcp_gt_50.tsv"
PIG_ENLOC_DIR <- "~/eqtl_selection/data/enloc/penloc"
PIG_PERM_DIR  <- "~/eqtl_selection/data/gtex/pgtex/eqtls_permutation/PigGTEx_v0.permutations_eQTL"
PIG_UCSC      <- "~/data/tss_ucsc/pig_encode_genes.txt"

CS_COVERAGE   <- 0.95   # credible-set coverage for the pig PIP_qtl sets, matching CS95

# Promoter / enhancer BEDs, and the cutoff ladder from the enrichment analysis. 0 is the
# default and headline value there (get_enr(cutoff = 0) -- strict overlap); 50 and 100 are
# its published sensitivity pair (enrichment_figure.ipynb cell 46). 250 appears only as a
# bed extension in the superseded enh_pro_enrichment.ipynb and is carried here for the
# ladder, not because the current pipeline uses it.
ELEM_DIR     <- "~/comparative_colocalization/data/intermediate_data/snps_elements_enrichments"
ELEM_CUTOFFS <- c(0L, 50L, 100L, 250L)

# The enrichment script's own null model: element bp over a hard-coded genome size.
GENOME_BP <- c(human = 3e9, pig = 2.7e9)

# An eQTL credible set whose nearest member variant lies within this many bp of a
# colocalizing position for the same gene-tissue is taken to BE the colocalizing eQTL
# rather than an additional, independent one. Fine-mapping splits one signal across
# neighbouring DAP-G clusters often enough that without this the coloc eQTL is counted
# again inside its own comparison background, which biases the coloc-vs-background
# contrast toward "no difference".
COLOC_SAME_BP <- 500

# eQTL effect sizes. Not in the fine-mapping output at all -- DAP-G and the .snp.out files
# carry PIPs only -- so `slope` is joined from the nominal significant-pair tables. That
# means an effect size exists only where the variant-gene pair cleared significance
# (~90% of credible-set members), and sets whose lead member did not are left NA.
HUMAN_EQTL_TAR  <- "~/data/gtex_v8/GTEx_Analysis_v8_eQTL.tar"
PIG_SIGNIF_DIR  <- "~/eqtl_selection/data/gtex/pgtex/eqtls/PigGTEx_v0.significant_eQTL"

# CS95 abbreviates these two; every other tissue name is already the GTEx form.
SKIN_REMAP <- c(Skin_Not_Sun_Exposed = "Skin_Not_Sun_Exposed_Suprapubic",
                Skin_Sun_Exposed     = "Skin_Sun_Exposed_Lower_leg")

# Resolve the script's own directory the way estimating_enloc_test_density.R does, so
# the output lands next to the script no matter where Rscript was invoked from.
script_dir <- function() {
    ca <- commandArgs(trailingOnly = FALSE)
    hit <- grep("^--file=", ca, value = TRUE)
    if (length(hit) == 0) return(normalizePath("."))
    normalizePath(dirname(sub("^--file=", "", hit[1])))
}

args   <- commandArgs(trailingOnly = TRUE)
HERE   <- script_dir()
OUTDIR <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else HERE
OUT    <- file.path(OUTDIR, "red_herring_eqtls.tsv.gz")

# ---------------------------------------------------------------------- helpers

pexp <- function(p) path.expand(p)

need_file <- function(p, what) {
    f <- pexp(p)
    if (!file.exists(f)) stop(what, " not found: ", f)
    f
}

need_dir <- function(p, what) {
    d <- pexp(p)
    if (!dir.exists(d)) stop(what, " not found: ", d)
    d
}

# awk programs are written to temp files and invoked with -f. Inlining them into a
# shell string would need three levels of quoting and is where this sort of script
# usually acquires a silent bug.
awk_file <- function(prog) {
    f <- tempfile(fileext = ".awk")
    writeLines(prog, f)
    f
}

key_file <- function(x) {
    f <- tempfile(fileext = ".txt")
    writeLines(x, f)
    f
}

# The two summaries of a credible set. wmedian sorts by POSITION (not by PIP) and walks
# the normalized cumulative weight to 0.5, which is the weighted median proper -- not
# "the position of the median-PIP variant".
wmedian <- function(pos, w) {
    o <- order(pos)
    pos <- pos[o]
    w <- w[o]
    s <- sum(w)
    if (!is.finite(s) || s <= 0) return(as.numeric(stats::median(pos)))
    pos[which(cumsum(w) / s >= 0.5)[1L]]
}

hr <- function() cat(strrep("-", 104), "\n")

# ------------------------------------------------------------ eQTL effect sizes

# Keyed on (tissue, gene, chr, pos) rather than on the variant id string: the two sources
# agree on position but a variant id also encodes ref/alt order, and one flipped pair would
# silently drop a row instead of erroring.
AWK_SLOPE <- 'BEGIN { while ((getline ln < KF) > 0) want[ln] = 1 }
$1 == "variant_id" || $1 == "phenotype_id" { next }
{
    if (MODE == "human") { v = $1; g = $2; sl = $8 }
    else                 { g = $1; v = $2; sl = $8 }
    sub(/\\.[0-9]+$/, "", g)
    n = split(v, p, "_")
    c = (MODE == "human") ? p[1] : "chr" p[1]
    k = TIS "\\t" g "\\t" c "\\t" p[2]
    if (!(k in want) || (k in seen)) next
    seen[k] = 1
    # Strip the sign as TEXT. Negating numerically routes the value through CONVFMT,
    # which defaults to %.6g and silently rounds. Positive slopes stay exact (the field
    # is never converted) while negative ones lose digits, so only half the data
    # degrades -- exactly the kind of asymmetry that never gets noticed.
    # (No apostrophes here: this block lives inside a single-quoted R string.)
    print k "\\t" (substr(sl, 1, 1) == "-" ? substr(sl, 2) : sl)
}'

# `leads` is (tissue, gene, chr, pos) for the top-PIP member of every set -- only those rows
# are wanted, so the very large nominal tables are filtered in awk and R never sees them.
read_slopes <- function(species, leads) {
    kf <- key_file(unique(leads[, paste(tissue, gene, chr, pos, sep = "\t")]))
    af <- awk_file(AWK_SLOPE)
    on.exit(unlink(c(kf, af)), add = TRUE)

    if (species == "human") {
        tarf <- need_file(HUMAN_EQTL_TAR, "GTEx v8 eQTL tar")
        tmpd <- file.path(tempdir(), "gtex_signif")
        dir.create(tmpd, showWarnings = FALSE, recursive = TRUE)
        # One pass over the 1.5 GB archive; extracting per tissue would rescan it 49 times.
        if (length(list.files(tmpd, pattern = "signif", recursive = TRUE)) == 0)
            system(sprintf("tar -xf %s -C %s --include='*signif_variant_gene_pairs*'",
                           shQuote(tarf), shQuote(tmpd)))
        files <- list.files(tmpd, pattern = "signif_variant_gene_pairs", recursive = TRUE,
                            full.names = TRUE)
        tis <- sub("[.]v8[.].*$", "", basename(files))
    } else {
        d <- need_dir(PIG_SIGNIF_DIR, "PigGTEx significant eQTL directory")
        files <- list.files(d, pattern = "[.]cis_qtl_pairs[.]significant[.]txt[.]gz$",
                            full.names = TRUE)
        tis <- sub("[.]cis_qtl_pairs.*$", "", basename(files))
    }
    if (length(files) == 0) stop("no significant-pair files found for ", species)

    keep <- tis %in% unique(leads$tissue)
    files <- files[keep]; tis <- tis[keep]
    out <- vector("list", length(files))
    for (i in seq_along(files)) {
        cmd <- sprintf("gzip -cd %s | awk -F'\\t' -v KF=%s -v TIS=%s -v MODE=%s -f %s",
                       shQuote(files[i]), shQuote(kf), shQuote(tis[i]), shQuote(species),
                       shQuote(af))
        r <- tryCatch(fread(cmd = cmd, header = FALSE, showProgress = FALSE,
                            col.names = c("tissue", "gene", "chr", "pos", "abs_slope")),
                      error = function(e) NULL)
        if (!is.null(r) && nrow(r) > 0) out[[i]] <- r
    }
    rbindlist(out)
}

# Lead member of a set: highest PIP, ties broken by position so the pick is deterministic.
# Returns one row per set carrying that member's position and PIP.
set_leads <- function(members, by, pipcol) {
    m <- members[order(-get(pipcol), pos)]
    # setdiff: `gene` is part of `by` for human but not for pig, and naming a by-column in
    # .SDcols emits it twice, which then breaks the downstream merge on a duplicate name.
    m[, .SD[1L], by = by, .SDcols = setdiff(c("gene", "chr", "pos", pipcol), by)]
}

# --------------------------------------------------- distance to the colocalizing locus

# Distance from a credible set to the nearest colocalizing position of its own gene-tissue,
# measured over the set's MEMBER variants rather than its median. The median is synthetic --
# it need not be any real variant -- and using it here would let a set that literally
# contains the coloc SNP score as far away. `link` are the gene-tissue columns shared with
# the coloc table; `by` identifies a set.
min_dist_to_coloc <- function(members, coloc_pos, by, link) {
    # unique(): `by` and `link` overlap (gene/tissue appear in both), and selecting a
    # duplicated name is a hard error in data.table.
    mm <- members[, unique(c(by, link, "pos")), with = FALSE]
    j <- merge(mm, coloc_pos, by = link, allow.cartesian = TRUE)
    j[, d := abs(pos - coloc_pos)]
    j[, .(dist_to_coloc = min(d)), by = by]
}

# ------------------------------------------------ promoter / enhancer element distance

# The cached element BEDs are 3-column, chr-prefixed, and already sorted AND merged (zero
# out-of-order and zero overlapping intervals in all four, verified). Merged + sorted is
# what makes the exact nearest-element distance two rolling joins instead of a bedtools
# shell-out, so the property is asserted rather than assumed.
read_elements <- function(species) {
    d <- need_dir(ELEM_DIR, "element BED directory")
    pre <- c(human = "h", pig = "p")[[species]]
    out <- lapply(c(pro = "pro", enh = "enh"), function(k) {
        f <- need_file(file.path(d, paste0(pre, k, ".tsv")), paste(species, k, "BED"))
        el <- fread(f, header = TRUE, sep = "\t", showProgress = FALSE,
                    select = c("chr", "start", "end"))
        setnames(el, c("chr", "s", "e"))
        el[, `:=`(chr = as.character(chr), s = as.numeric(s), e = as.numeric(e))]
        setorder(el, chr, s)
        bad <- el[, sum(s < shift(e, fill = -Inf) & chr == shift(chr, fill = "")), ]
        if (bad > 0)
            stop(sprintf("%s %s BED is not merged (%d overlapping intervals); the distance
                          join assumes disjoint sorted intervals", species, k, bad))
        setkey(el, chr, s)
        el[]
    })
    out
}

# Signed distance to the nearest element, reproducing `bedtools closest -D ref` exactly.
# Reverse-engineered from the cached pro_dist/enh_dist columns rather than from the docs,
# then checked against them: 100.0000% exact on all 350,370 available rows (henloc 65,605,
# penloc 9,646, heqtls 99,934, both element types), zero mismatches.
#
# For a 1-based `pos` and a BED element [s, e) -- which covers 1-based s+1 .. e:
#     overlap        (s <  pos <= e)  ->  0
#     element left        (e < pos)   ->  e - pos       (negative)
#     element right      (s >= pos)   ->  s + 1 - pos   (positive)
# and the nearer of the two candidates wins on absolute value.
#
# THE `s < pos` IS STRICT AND LOAD-BEARING. Written `s <= pos` it treats a variant sitting
# one base before an element as inside it; that one character produced 14 mismatches
# against the cached tables, all of the form "bedtools 1, mine 0".
elem_dist <- function(chr, pos, el) {
    q <- data.table(chr = as.character(chr), pos = as.numeric(pos))
    iL <- el[q, on = .(chr, s <  pos), mult = "last",  which = TRUE]
    iR <- el[q, on = .(chr, s >= pos), mult = "first", which = TRUE]
    dL <- ifelse(is.na(iL), NA_real_, ifelse(el$e[iL] >= q$pos, 0, el$e[iL] - q$pos))
    dR <- ifelse(is.na(iR), NA_real_, el$s[iR] + 1 - q$pos)
    ifelse(is.na(dL), dR, ifelse(is.na(dR), dL, ifelse(abs(dL) <= abs(dR), dL, dR)))
}

# Per credible set, the fraction of its MEMBER variants within each cutoff. These cannot be
# recovered from the finished file the way the summary-position flags can (the member
# positions are not in it), which is why all four cutoffs are precomputed here.
member_element_fracs <- function(members, el, by) {
    m <- members[, c(by, "chr", "pos"), with = FALSE]
    m[, `:=`(.pd = elem_dist(chr, pos, el$pro), .ed = elem_dist(chr, pos, el$enh))]
    j <- c(sprintf("frac_mem_pro_%d = mean(abs(.pd) <= %d)", ELEM_CUTOFFS, ELEM_CUTOFFS),
           sprintf("frac_mem_enh_%d = mean(abs(.ed) <= %d)", ELEM_CUTOFFS, ELEM_CUTOFFS))
    m[, eval(parse(text = paste0(".(", paste(j, collapse = ", "), ")"))), by = by]
}

# ------------------------------------------------------------ TSS reconstruction

# Both projects publish tss_distance = variant_pos - TSS_pos, so one row per gene is
# enough to recover TSS_pos exactly. `want` keeps these tables to the genes actually
# needed, which is what makes a 9.8M-row and a 1.6 GB file cheap to mine.
AWK_TSS_INVERT <- 'BEGIN { while ((getline ln < GF) > 0) want[ln] = 1 }
$1 == "phenotype_id" { next }
{
    g = $(GCOL)
    sub(/\\.[0-9]+$/, "", g)
    if (!(g in want) || (g in s)) next
    split($(VCOL), v, "_")
    s[g] = v[2] - $(DCOL)
}
END { for (g in s) print g "\\t" s[g] }'

tss_from_distance <- function(cmd_in, genes, gcol, vcol, dcol) {
    gf <- key_file(genes)
    af <- awk_file(AWK_TSS_INVERT)
    cmd <- sprintf("%s | awk -F'\\t' -v GF=%s -v GCOL=%d -v VCOL=%d -v DCOL=%d -f %s",
                   cmd_in, shQuote(gf), gcol, vcol, dcol, shQuote(af))
    out <- fread(cmd = cmd, header = FALSE, col.names = c("gene", "tss_pos"),
                 colClasses = c("character", "numeric"), showProgress = FALSE)
    unlink(c(gf, af))
    out
}

# UCSC genePred for pig: TSS is txStart+1 on +, txEnd on -, taken over all transcripts
# of the gene. name2 (field 13) is a VERSIONED ENSSSCG.
AWK_UCSC_PIG <- '$1 ~ /^#/ { next }
{
    g = $13
    sub(/\\.[0-9]+$/, "", g)
    t = ($4 == "+") ? $5 + 1 : $6
    if (!(g in s)) { s[g] = t; st[g] = $4 }
    else if (st[g] == "+" && t < s[g]) s[g] = t
    else if (st[g] == "-" && t > s[g]) s[g] = t
}
END { for (g in s) print g "\\t" s[g] }'

build_human_tss <- function(genes) {
    a <- tss_from_distance(sprintf("gzip -cd %s", shQuote(need_file(HUMAN_TSS_A, "human TSS key A"))),
                           genes, gcol = 3, vcol = 2, dcol = 4)
    missing <- setdiff(genes, a$gene)
    if (length(missing) > 0 && file.exists(pexp(HUMAN_TSS_B))) {
        message(sprintf("  %d human genes absent from the fine-mapped key; scanning all significant pairs",
                        length(missing)))
        b <- tss_from_distance(sprintf("gzip -cd %s", shQuote(pexp(HUMAN_TSS_B))),
                               missing, gcol = 3, vcol = 2, dcol = 4)
        a <- rbindlist(list(a, b))
    }
    unique(a, by = "gene")
}

build_pig_tss <- function(genes) {
    perm <- need_dir(PIG_PERM_DIR, "PigGTEx permutation directory")
    files <- list.files(perm, pattern = "[.]cis_qtl_fdr0[.]05[.]txt[.]gz$", full.names = TRUE)
    if (length(files) == 0) stop("no *.cis_qtl_fdr0.05.txt.gz under: ", perm)
    a <- tss_from_distance(sprintf("gzip -cd %s", paste(shQuote(files), collapse = " ")),
                           genes, gcol = 1, vcol = 7, dcol = 8)

    af <- awk_file(AWK_UCSC_PIG)
    ucsc <- fread(cmd = sprintf("awk -F'\\t' -f %s %s", shQuote(af),
                                shQuote(need_file(PIG_UCSC, "pig UCSC gene models"))),
                  header = FALSE, col.names = c("gene", "tss_pos"),
                  colClasses = c("character", "numeric"), showProgress = FALSE)
    unlink(af)

    # The UCSC fill is only legitimate if the two sources agree where they overlap.
    # They were checked at 17,431/17,431 exact during development; keep the check so a
    # future PigGTEx or UCSC refresh fails loudly instead of mixing two conventions.
    chk <- merge(a, ucsc, by = "gene", suffixes = c("_perm", "_ucsc"))
    if (nrow(chk) > 0) {
        disagree <- chk[tss_pos_perm != tss_pos_ucsc]
        if (nrow(disagree) > 0)
            stop(sprintf(paste0("pig TSS sources disagree on %d of %d shared genes (e.g. %s: ",
                                "PigGTEx %.0f vs UCSC %.0f). The UCSC fill is only valid while ",
                                "these agree exactly -- investigate before trusting tss_dist."),
                         nrow(disagree), nrow(chk), disagree$gene[1],
                         disagree$tss_pos_perm[1], disagree$tss_pos_ucsc[1]))
        message(sprintf("  pig TSS cross-check: %d shared genes, all exact", nrow(chk)))
    }

    rbindlist(list(a, ucsc[gene %in% setdiff(genes, a$gene)]))[, .(tss_pos = tss_pos[1]), by = gene]
}

# ------------------------------------------------------------------ human coloc

read_human_coloc <- function() {
    f <- need_file(HUMAN_COLOC, "human coloc table")
    d <- fread(f, header = FALSE, sep = "\t", showProgress = FALSE,
               col.names = c("gwas_locus", "molecular_qtl_trait", "locus_gwas_pip",
                             "locus_rcp", "lead_coloc_SNP", "lead_snp_rcp", "file"))

    stem <- sub("[.]enloc[.]rst[.]gz$", "", d$file)
    bad <- !grepl("__PM__", stem, fixed = TRUE)
    if (any(bad))
        stop("source filenames not in <TRAIT>__PM__<TISSUE> form: ",
             paste(head(unique(stem[bad]), 3), collapse = ", "))

    # This file is already a strict RCP > 0.5 cut (min 0.501 as shipped). Re-filtering
    # would be a no-op; asserting catches the day someone points HUMAN_COLOC at the
    # rcp_10 sibling that lives in the same directory.
    if (min(d$locus_rcp) <= 0.5)
        stop(sprintf("HUMAN_COLOC is not thresholded at RCP > 0.5 (min %.4f) -- wrong file?",
                     min(d$locus_rcp)))

    v <- tstrsplit(d$lead_coloc_SNP, "_", fixed = TRUE)
    d[, `:=`(trait  = sub("__PM__.*$", "", stem),
             tissue = sub("^.*__PM__", "", stem),
             gene   = sub("[.][0-9]+$", "", molecular_qtl_trait),
             chr    = v[[1]],
             pos    = as.numeric(v[[2]]),
             rcp    = locus_rcp,
             pip    = lead_snp_rcp,
             cs_id  = gwas_locus)]
    d[, .(species = "human", trait, tissue, gene, chr, pos, rcp, pip, cs_id)]
}

# CS95 is VCF-shaped: field 6 is a "|"-separated list of <gene>.<ver>:<cluster>@<tissue>=<PIP>.
# The explode is done in awk and filtered to the coloc gene-tissue pairs on the way past,
# so R never sees the 7.6M-record full expansion.
AWK_CS95 <- 'BEGIN { while ((getline ln < KF) > 0) want[ln] = 1 }
{
    n = split($6, r, "|")
    for (i = 1; i <= n; i++) {
        if (r[i] == "") continue
        split(r[i], a, "@")
        split(a[2], b, "=")
        traw = b[1]
        t = traw
        if (t == "Skin_Not_Sun_Exposed") t = "Skin_Not_Sun_Exposed_Suprapubic"
        else if (t == "Skin_Sun_Exposed") t = "Skin_Sun_Exposed_Lower_leg"
        split(a[1], s, ":")
        g = s[1]
        sub(/\\.[0-9]+$/, "", g)
        if (!((g "\\t" t) in want)) continue
        print g "\\t" s[2] "\\t" t "\\t" traw "\\t" "chr" $1 "\\t" $2 "\\t" b[2]
    }
}'

read_human_credible_sets <- function(coloc, el) {
    f  <- need_file(HUMAN_CS95, "GTEx DAP-G CS95 file")
    kf <- key_file(unique(paste(coloc$gene, coloc$tissue, sep = "\t")))
    af <- awk_file(AWK_CS95)
    message("Exploding CS95 credible sets (1.4M rows) ...")
    long <- fread(cmd = sprintf("gzip -cd %s | awk -F'\\t' -v KF=%s -f %s",
                                shQuote(f), shQuote(kf), shQuote(af)),
                  header = FALSE, showProgress = FALSE,
                  col.names = c("gene", "cluster", "tissue", "tissue_raw", "chr", "pos", "pip"),
                  colClasses = c("character", "character", "character", "character",
                                 "character", "numeric", "numeric"))
    unlink(c(kf, af))
    if (nrow(long) == 0) stop("CS95 explode produced no rows -- check the annotation field layout")

    # Assert the Skin remap is the only name surgery needed, and that it is needed.
    rewritten <- unique(long[tissue_raw != tissue]$tissue_raw)
    if (!setequal(rewritten, names(SKIN_REMAP)))
        stop("CS95 tissue remap changed unexpectedly. Rewritten: ",
             paste(rewritten, collapse = ", "), "; expected: ",
             paste(names(SKIN_REMAP), collapse = ", "))
    stray <- setdiff(unique(long$tissue), unique(coloc$tissue))
    if (length(stray) > 0)
        stop("CS95 tissue names not present in the coloc table after remap: ",
             paste(head(stray, 5), collapse = ", "))

    # Mark members that ARE a colocalizing position for their gene-tissue, before the
    # sets are collapsed and the member positions are gone.
    coloc_keys <- unique(paste(coloc$gene, coloc$tissue, coloc$pos))
    long[, is_member := paste(gene, tissue, pos) %chin% coloc_keys]

    sets <- long[, .(chr          = chr[1],
                     pos_median   = round(as.numeric(stats::median(pos))),
                     pos_wmedian  = wmedian(pos, pip),
                     cs_size      = .N,
                     cs_pip       = sum(pip),
                     is_coloc_set = any(is_member),
                     cs_truncated = FALSE),
                 by = .(gene, tissue, cluster)]
    sets[, cs_id := paste(gene, cluster, sep = ":")]
    sets[, species := "human"]

    # Member-level element overlap, computed while the member positions still exist.
    mf <- member_element_fracs(long, el, by = c("gene", "tissue", "cluster"))
    sets <- merge(sets, mf, by = c("gene", "tissue", "cluster"), all.x = TRUE)

    # Distance to the nearest colocalizing position of the same gene-tissue. Coloc
    # positions are deduped first: several traits routinely colocalize at the same SNP.
    cp <- unique(coloc[, .(gene, tissue, coloc_pos = pos)])
    dc <- min_dist_to_coloc(long, cp, by = c("gene", "tissue", "cluster"),
                            link = c("gene", "tissue"))
    sets <- merge(sets, dc, by = c("gene", "tissue", "cluster"), all.x = TRUE)

    ld <- set_leads(long, by = c("gene", "tissue", "cluster"), pipcol = "pip")
    setnames(ld, c("pos", "pip"), c("lead_pos", "lead_pip"))
    sets <- merge(sets, ld[, .(gene, tissue, cluster, lead_pos, lead_pip)],
                  by = c("gene", "tissue", "cluster"), all.x = TRUE)
    sets[]
}

# -------------------------------------------------------------------- pig coloc

read_pig_coloc <- function() {
    f <- need_file(PIG_COLOC, "pig coloc table")
    # Not tab-delimited: data rows are space-padded with a tab only before the filename.
    # Squeezing all whitespace runs to single tabs is safe because no field contains a space.
    d <- fread(cmd = sprintf("tr -s ' \\t' '\\t' < %s", shQuote(f)),
               header = FALSE, sep = "\t", showProgress = FALSE,
               col.names = c("Signal", "Num_SNP", "CPIP_qtl", "CPIP_gwas_marginal",
                             "CPIP_gwas_qtl_prior", "RCP", "file"))
    if (min(d$RCP) < 0.5)
        stop(sprintf("PIG_COLOC is not thresholded at RCP > 0.5 (min %.4f) -- wrong file?",
                     min(d$RCP)))

    # Signal is gene:chr:start:end:cluster(@)<gwas_locus>. Split on the literal "(@)"
    # with index/substr: in awk and in regex generally, "(@)" is a GROUP matching "@",
    # so a regex split leaves a stray "(" on the left half.
    at <- regexpr("(@)", d$Signal, fixed = TRUE)
    d[, eqtl_signal := ifelse(at > 0, substr(Signal, 1, at - 1), Signal)]
    d[, `:=`(gene   = sub(":.*$", "", eqtl_signal),
             tissue = sub("/.*$", "", file),
             trait  = sub("[.]enloc[.]sig[.]out$", "", sub("^.*/[^.]+[.]", "", file)))]
    d[]
}

# One pass per tissue over the .enloc.snp.out stream, producing two small files:
#   OUTA  distinct (tissue, eqtl_signal, SNP, PIP_qtl) -- the eQTL fine-mapping
#   OUTB  per (trait, full signal) the SNP with the highest SCP -- the coloc position.
#         Keyed on trait AND signal, not signal alone: the same eQTL signal is tested
#         against many traits, and a coloc row exists only for the pairs that cleared
#         RCP > 0.5. Matching on signal alone emitted 245 positions for Milk's 2 rows.
# The dedup for OUTA happens here rather than in R because the raw match volume is
# ~28x the deduped volume (Muscle: 15,254,088 rows -> 545,039 pairs).
AWK_PIG_SNP <- 'BEGIN {
    while ((getline ln < GF) > 0) keep[ln] = 1
    while ((getline ln < SF) > 0) colocsig[ln] = 1
    nm = 0
    while ((getline ln < ML) > 0) { nm++; member[nm] = ln }
    idx = 0
}
$1 == "Signal" {
    idx++
    tr = member[idx]
    sub(/.*\\//, "", tr)
    sub(/^[^.]*\\./, "", tr)
    sub(/\\.enloc\\.snp\\.out$/, "", tr)
    next
}
{
    sf = $1
    p = index(sf, "(@)")
    es = (p > 0) ? substr(sf, 1, p - 1) : sf
    split(es, b, ":")
    k = tr "\\t" sf
    if (k in colocsig) {
        if (!(k in best) || $6 + 0 > best[k]) { best[k] = $6 + 0; bestsnp[k] = $2 }
    }
    if (b[1] in keep) {
        k2 = es SUBSEP $2
        if (!(k2 in seen)) { seen[k2] = 1; print TIS "\\t" es "\\t" $2 "\\t" $3 > OUTA }
    }
}
END {
    for (k in best) {
        split(k, q, "\\t")
        print TIS "\\t" q[1] "\\t" q[2] "\\t" bestsnp[k] "\\t" best[k] > OUTB
    }
    if (idx != nm) {
        printf "member/header mismatch in %s: %d headers streamed, %d members listed\\n", TIS, idx, nm > "/dev/stderr"
        exit 3
    }
}'

# macOS bsdtar. --include selects members inside the archive; the empty-listing check
# below is what fails loudly if a different tar ignores it.
# `colockeys` are "<trait>\t<full Signal>" strings, not bare signal ids -- see AWK_PIG_SNP.
stream_pig_tissue <- function(tissue, genes, colockeys, af) {
    base  <- need_dir(PIG_ENLOC_DIR, "pig enloc directory")
    dir_t <- file.path(base, tissue)
    tar_t <- file.path(base, paste0(tissue, ".fastENLOC.tar.gz"))

    if (dir.exists(dir_t)) {
        list_cmd   <- sprintf("ls -1 %s/*.enloc.snp.out", shQuote(dir_t))
        stream_cmd <- sprintf("cat %s/*.enloc.snp.out", shQuote(dir_t))
    } else if (file.exists(tar_t)) {
        list_cmd   <- sprintf("tar -tzf %s --include='*.enloc.snp.out'", shQuote(tar_t))
        stream_cmd <- sprintf("tar -xOzf %s --include='*.enloc.snp.out'", shQuote(tar_t))
    } else {
        stop("no extracted directory or tarball for pig tissue: ", tissue)
    }

    ml <- tempfile(fileext = ".txt")
    if (system(sprintf("%s > %s", list_cmd, shQuote(ml))) != 0)
        stop("could not list .enloc.snp.out members for tissue: ", tissue)
    n_members <- length(readLines(ml))
    if (n_members == 0) stop("no .enloc.snp.out members found for tissue: ", tissue)

    gf <- key_file(genes)
    sf <- key_file(colockeys)
    outa <- tempfile(fileext = ".tsv")
    outb <- tempfile(fileext = ".tsv")
    file.create(outa, outb)

    status <- system(sprintf("%s | awk -v TIS=%s -v GF=%s -v SF=%s -v ML=%s -v OUTA=%s -v OUTB=%s -f %s",
                             stream_cmd, shQuote(tissue), shQuote(gf), shQuote(sf),
                             shQuote(ml), shQuote(outa), shQuote(outb), shQuote(af)))
    if (status != 0)
        stop("pig .snp.out pass failed for tissue ", tissue, " (awk status ", status, ")")

    rd <- function(p, nm) {
        if (file.info(p)$size == 0) return(NULL)
        fread(p, header = FALSE, sep = "\t", col.names = nm, showProgress = FALSE)
    }
    out <- list(snps  = rd(outa, c("tissue", "eqtl_signal", "snp", "pip_qtl")),
                coloc = rd(outb, c("tissue", "trait", "Signal", "snp", "scp")),
                n_members = n_members)
    unlink(c(ml, gf, sf, outa, outb))
    out
}

read_pig <- function(coloc) {
    af <- awk_file(AWK_PIG_SNP)
    tissues <- sort(unique(coloc$tissue))
    message(sprintf("Streaming pig .enloc.snp.out for %d tissues (the slow step) ...", length(tissues)))
    snps <- vector("list", length(tissues))
    hits <- vector("list", length(tissues))
    for (i in seq_along(tissues)) {
        ti <- tissues[i]
        sub <- coloc[tissue == ti]
        r <- stream_pig_tissue(ti, unique(sub$gene),
                               unique(paste(sub$trait, sub$Signal, sep = "\t")), af)
        snps[[i]] <- r$snps
        hits[[i]] <- r$coloc
        message(sprintf("  [%2d/%2d] %-20s %d members, %s eQTL (signal, SNP) pairs, %s coloc positions",
                        i, length(tissues), ti, r$n_members,
                        format(if (is.null(r$snps)) 0L else nrow(r$snps), big.mark = ","),
                        format(if (is.null(r$coloc)) 0L else nrow(r$coloc), big.mark = ",")))
    }
    unlink(af)
    list(snps = rbindlist(snps), coloc = rbindlist(hits))
}

# The pig analogue of CS95: sort a signal's SNPs by PIP_qtl and take members until the
# cumulative PIP reaches CS_COVERAGE. Where a signal's PIP_qtl never sums that high --
# because .snp.out only lists the SNPs fastEnloc evaluated -- take every member and say so.
build_pig_sets <- function(snps, el, coloc_pos) {
    v <- tstrsplit(snps$snp, "_", fixed = TRUE)
    snps[, `:=`(chr = paste0("chr", v[[1]]), pos = as.numeric(v[[2]]))]
    setorder(snps, tissue, eqtl_signal, -pip_qtl)
    snps[, `:=`(rank = seq_len(.N),
                cum  = cumsum(pip_qtl),
                tot  = sum(pip_qtl)), by = .(tissue, eqtl_signal)]
    snps[, cutoff := { w <- which(cum >= CS_COVERAGE)[1L]; if (is.na(w)) .N else w },
         by = .(tissue, eqtl_signal)]

    sets <- snps[rank <= cutoff,
                 .(chr          = chr[1],
                   pos_median   = round(as.numeric(stats::median(pos))),
                   pos_wmedian  = wmedian(pos, pip_qtl),
                   cs_size      = .N,
                   cs_pip       = sum(pip_qtl),
                   cs_truncated = tot[1] < CS_COVERAGE),
                 by = .(tissue, eqtl_signal)]
    sets[, `:=`(species = "pig",
                gene    = sub(":.*$", "", eqtl_signal),
                cs_id   = eqtl_signal)]

    # Same member-level pass, over exactly the members that made the 95% cut, so the
    # fractions and cs_size describe the same set of variants.
    mem <- snps[rank <= cutoff]
    mf <- member_element_fracs(mem, el, by = c("tissue", "eqtl_signal"))
    sets <- merge(sets, mf, by = c("tissue", "eqtl_signal"), all.x = TRUE)

    mem[, gene := sub(":.*$", "", eqtl_signal)]
    dc <- min_dist_to_coloc(mem, coloc_pos, by = c("tissue", "eqtl_signal"),
                            link = c("tissue", "gene"))
    sets <- merge(sets, dc, by = c("tissue", "eqtl_signal"), all.x = TRUE)

    ld <- set_leads(mem, by = c("tissue", "eqtl_signal"), pipcol = "pip_qtl")
    setnames(ld, c("pos", "pip_qtl"), c("lead_pos", "lead_pip"))
    sets <- merge(sets, ld[, .(tissue, eqtl_signal, lead_pos, lead_pip)],
                  by = c("tissue", "eqtl_signal"), all.x = TRUE)
    sets[]
}

# ------------------------------------------------------------------- is_coloc_set

# Most colocalizing positions are literal members of one of their gene-tissue's credible
# sets, and `any(is_member)` has already caught those. A colocalizing SNP can fall
# outside every 95% set, though (CS95 keeps only the top 95% of each signal's mass), so
# those gene-tissue pairs fall back to the set whose median sits closest. Marked the same
# way either route, because in both cases it is "the set this colocalization came from".
mark_nearest_coloc_set <- function(sets, coloc) {
    pairs_with_sets <- unique(sets[, .(gene, tissue)])
    already <- sets[is_coloc_set == TRUE, unique(paste(gene, tissue))]
    need <- merge(coloc[, .(gene, tissue, cpos = pos)], pairs_with_sets, by = c("gene", "tissue"))
    need <- need[!paste(gene, tissue) %chin% already]
    if (nrow(need) == 0) return(0L)

    cand <- merge(sets[, .(gene, tissue, cs_id, pos_median)], need,
                  by = c("gene", "tissue"), allow.cartesian = TRUE)
    cand[, d := abs(pos_median - cpos)]
    pick <- cand[cand[, .I[which.min(d)], by = .(gene, tissue)]$V1]
    sets[pick, on = .(gene, tissue, cs_id), is_coloc_set := TRUE]
    nrow(pick)
}

# --------------------------------------------------------------------------- run

message("Reading colocalizations ...")
hcoloc <- read_human_coloc()
pcoloc <- read_pig_coloc()
message(sprintf("  human %s coloc rows, pig %s coloc rows",
                format(nrow(hcoloc), big.mark = ","), format(nrow(pcoloc), big.mark = ",")))

message("Reading promoter / enhancer element BEDs ...")
ELEM <- list(human = read_elements("human"), pig = read_elements("pig"))

hsets <- read_human_credible_sets(hcoloc, ELEM$human)
n_human_fallback <- mark_nearest_coloc_set(hsets, hcoloc)

pig <- read_pig(pcoloc)

# Pig enloc rows: the coloc table joined to the top-SCP SNP of its own signal. Built
# BEFORE the credible sets, because their distance-to-coloc needs these positions.
pig_enloc <- merge(pcoloc, pig$coloc, by = c("tissue", "trait", "Signal"), all.x = TRUE)
pv <- tstrsplit(pig_enloc$snp, "_", fixed = TRUE)
pig_enloc[, `:=`(species = "pig",
                 chr     = ifelse(is.na(snp), NA_character_, paste0("chr", pv[[1]])),
                 pos     = as.numeric(pv[[2]]),
                 rcp     = RCP,
                 pip     = scp,
                 cs_id   = Signal)]
n_pig_nopos <- sum(is.na(pig_enloc$pos))

psets <- build_pig_sets(pig$snps, ELEM$pig,
                        unique(pig_enloc[!is.na(pos), .(tissue, gene, coloc_pos = pos)]))
psets[, is_coloc_set := paste(tissue, eqtl_signal) %chin%
          unique(paste(pcoloc$tissue, pcoloc$eqtl_signal))]

# Effect sizes for the lead member of every set. Joined after both species' sets exist so
# the very large nominal tables are read once each, filtered to just the leads.
message("Joining eQTL effect sizes to set lead variants ...")
attach_slopes <- function(sets, species) {
    sl <- read_slopes(species, sets[, .(tissue, gene, chr, pos = lead_pos)])
    if (nrow(sl) == 0) { sets[, lead_abs_slope := NA_real_]; return(invisible(sets)) }
    setnames(sl, "pos", "lead_pos")
    sets[sl, on = .(tissue, gene, chr, lead_pos), lead_abs_slope := i.abs_slope]
    if (!"lead_abs_slope" %in% names(sets)) sets[, lead_abs_slope := NA_real_]
    invisible(sets)
}
attach_slopes(hsets, "human")
attach_slopes(psets, "pig")
message(sprintf("  lead effect size present: human %s / %s sets, pig %s / %s",
                format(hsets[!is.na(lead_abs_slope), .N], big.mark = ","),
                format(nrow(hsets), big.mark = ","),
                format(psets[!is.na(lead_abs_slope), .N], big.mark = ","),
                format(nrow(psets), big.mark = ",")))

message("Reconstructing TSS positions ...")
htss <- build_human_tss(unique(hcoloc$gene))
ptss <- build_pig_tss(unique(pcoloc$gene))
message(sprintf("  human TSS for %s / %s coloc genes; pig %s / %s",
                format(nrow(htss), big.mark = ","), format(uniqueN(hcoloc$gene), big.mark = ","),
                format(nrow(ptss), big.mark = ","), format(uniqueN(pcoloc$gene), big.mark = ",")))

# ------------------------------------------------------------------ assemble rows

WIDE  <- setdiff(ELEM_CUTOFFS, 0L)
# Summary-position columns. These stay recoverable from pro_dist/enh_dist afterwards --
# they are precomputed for convenience, not necessity.
ELEM_COLS <- c("pro_dist", "enh_dist", "in_pro", "in_enh",
               paste0("near_pro_", WIDE), paste0("near_enh_", WIDE))
# Member-level columns. NOT recoverable from the finished file, because the member
# positions are not in it.
FRAC_COLS <- c(paste0("frac_mem_pro_", ELEM_CUTOFFS), paste0("frac_mem_enh_", ELEM_CUTOFFS))

COLS <- c("species", "trait", "tissue", "gene", "type", "tss_dist", "abs_tss_dist",
          "chr", "pos", "rcp", "pip", "cs_id", "cs_size", "is_coloc_set", "cs_truncated",
          "dist_to_coloc", "dist_to_coloc_pos", "same_as_coloc",
          "lead_pos", "lead_pip", "lead_abs_slope",
          ELEM_COLS, FRAC_COLS)

as_enloc <- function(d) {
    out <- d[, .(species, trait, tissue, gene, type = "enloc",
                 chr, pos, rcp = as.numeric(rcp), pip = as.numeric(pip), cs_id,
                 cs_size = NA_integer_, is_coloc_set = NA, cs_truncated = NA)]
    # An enloc row is one real variant, so in_pro/in_enh already say everything a member
    # fraction could; leaving these NA keeps "fraction of members" meaning only that.
    out[, (FRAC_COLS) := NA_real_]
    out[, `:=`(dist_to_coloc = NA_real_, same_as_coloc = NA,
               lead_pos = NA_real_, lead_pip = NA_real_, lead_abs_slope = NA_real_)]
    out[]
}

as_eqtl <- function(d, stat) {
    p <- if (stat == "median") d$pos_median else d$pos_wmedian
    out <- data.table(species = d$species, trait = NA_character_, tissue = d$tissue,
                      gene = d$gene, type = paste0("eqtl_", stat), chr = d$chr,
                      pos = as.numeric(p), rcp = NA_real_, pip = as.numeric(d$cs_pip),
                      cs_id = d$cs_id, cs_size = as.integer(d$cs_size),
                      is_coloc_set = d$is_coloc_set, cs_truncated = d$cs_truncated)
    # Identical across a set's median and wmedian rows -- both summarize the same members.
    out[, (FRAC_COLS) := d[, FRAC_COLS, with = FALSE]]
    # Set-level, so identical across a set's median and wmedian rows.
    out[, dist_to_coloc := as.numeric(d$dist_to_coloc)]
    out[, same_as_coloc := d$is_coloc_set | (!is.na(dist_to_coloc) & dist_to_coloc <= COLOC_SAME_BP)]
    out[, `:=`(lead_pos       = as.numeric(d$lead_pos),
               lead_pip       = as.numeric(d$lead_pip),
               lead_abs_slope = as.numeric(d$lead_abs_slope))]
    out[]
}

# Distances are computed once per distinct (species, chr, pos) -- 292,389 rows collapse to
# 170,412 positions -- then joined back.
annotate_elements <- function(final) {
    u <- rbindlist(lapply(names(ELEM), function(sp) {
        q <- unique(final[species == sp & !is.na(pos), .(species, chr, pos)])
        q[, `:=`(pro_dist = elem_dist(chr, pos, ELEM[[sp]]$pro),
                 enh_dist = elem_dist(chr, pos, ELEM[[sp]]$enh))]
        q
    }))
    final[u, on = .(species, chr, pos),
          `:=`(pro_dist = i.pro_dist, enh_dist = i.enh_dist)]
    final[, `:=`(in_pro = abs(pro_dist) <= 0L, in_enh = abs(enh_dist) <= 0L)]
    for (k in WIDE) {
        final[, (paste0("near_pro_", k)) := abs(pro_dist) <= k]
        final[, (paste0("near_enh_", k)) := abs(enh_dist) <= k]
    }
    # Every chromosome in the table is present in every element BED, so a missing distance
    # would mean the join silently failed rather than that the data is genuinely absent.
    if (final[!is.na(pos) & (is.na(pro_dist) | is.na(enh_dist)), .N] > 0)
        stop("element distance is NA for rows that have a position -- chromosome naming mismatch?")
    final
}

final <- rbindlist(list(
    as_enloc(hcoloc), as_enloc(pig_enloc),
    as_eqtl(hsets, "median"), as_eqtl(hsets, "wmedian"),
    as_eqtl(psets, "median"), as_eqtl(psets, "wmedian")))

tss <- rbindlist(list(cbind(species = "human", htss), cbind(species = "pig", ptss)))
final[tss, on = .(species, gene), tss_pos := i.tss_pos]
final[, tss_dist := pos - tss_pos]
final[, tss_pos := NULL]

# tss_dist is signed by GENOME COORDINATE, not by strand (see the header). abs_tss_dist is
# carried explicitly so that "how far is the TSS" never depends on remembering that.
final[, abs_tss_dist := abs(tss_dist)]

# The same-eQTL test at this row's own position, alongside the member-based dist_to_coloc
# that drives the flag. They differ because a summary position is synthetic.
cpos <- unique(final[type == "enloc" & !is.na(pos), .(species, tissue, gene, coloc_pos = pos)])
final[, .rid := .I]
rowd <- merge(final[type != "enloc" & !is.na(pos), .(.rid, species, tissue, gene, pos)],
              cpos, by = c("species", "tissue", "gene"), allow.cartesian = TRUE)
rowd <- rowd[, .(dist_to_coloc_pos = min(abs(pos - coloc_pos))), by = .rid]
final[, dist_to_coloc_pos := NA_real_]
final[rowd, on = ".rid", dist_to_coloc_pos := i.dist_to_coloc_pos]
final[, .rid := NULL]

final <- annotate_elements(final)
setcolorder(final, COLS)
setorder(final, species, tissue, gene, type, cs_id)

fwrite(final, OUT, sep = "\t", quote = FALSE, compress = "gzip", na = "NA")

# ---------------------------------------------------------------- terminal report

n_na_tss  <- final[is.na(tss_dist), .N]
na_genes  <- final[is.na(tss_dist), uniqueN(gene), by = species]
h_pairs   <- uniqueN(hcoloc, by = c("gene", "tissue"))
h_covered <- uniqueN(hsets, by = c("gene", "tissue"))
p_pairs   <- uniqueN(pcoloc, by = c("gene", "tissue"))
p_covered <- uniqueN(psets[, .(gene, tissue)])

cat("\n")
hr()
cat("COLOCALIZING VARIANTS vs FINE-MAPPED eQTLs OF THE SAME GENE-TISSUE\n")
hr()
cat(sprintf("human  coloc rows %s | traits %d | tissues %d | genes %s | gene-tissue pairs %s\n",
            format(nrow(hcoloc), big.mark = ","), uniqueN(hcoloc$trait), uniqueN(hcoloc$tissue),
            format(uniqueN(hcoloc$gene), big.mark = ","), format(h_pairs, big.mark = ",")))
cat(sprintf("       CS95 credible sets %s over %s pairs (%.1f%%); %s pairs have none\n",
            format(nrow(hsets), big.mark = ","), format(h_covered, big.mark = ","),
            100 * h_covered / h_pairs, format(h_pairs - h_covered, big.mark = ",")))
cat(sprintf("       coloc set identified by membership for %s sets, by nearest median for %s pairs\n",
            format(hsets[is_coloc_set == TRUE, .N] - n_human_fallback, big.mark = ","),
            format(n_human_fallback, big.mark = ",")))
cat(sprintf("pig    coloc rows %s | traits %d | tissues %d | genes %s | gene-tissue pairs %s\n",
            format(nrow(pcoloc), big.mark = ","), uniqueN(pcoloc$trait), uniqueN(pcoloc$tissue),
            format(uniqueN(pcoloc$gene), big.mark = ","), format(p_pairs, big.mark = ",")))
cat(sprintf("       PIP_qtl credible sets %s over %s pairs (%.1f%%); %s truncated below %.2f coverage\n",
            format(nrow(psets), big.mark = ","), format(p_covered, big.mark = ","),
            100 * p_covered / p_pairs, format(psets[cs_truncated == TRUE, .N], big.mark = ","),
            CS_COVERAGE))
if (n_pig_nopos > 0)
    cat(sprintf("       %s coloc rows had no .snp.out row for their signal (position NA)\n",
                format(n_pig_nopos, big.mark = ",")))
# Truncated pig sets are weak signals, not merely clipped ones, so show the split that
# tells the two apart rather than leaving the 64% figure to be read as a defect.
ptr <- psets[, .(sets = .N, median_signal_pip = round(stats::median(cs_pip), 3),
                 median_set_size = stats::median(cs_size)), by = cs_truncated]
setorder(ptr, cs_truncated)
cat("       pig sets by truncation (filter to FALSE for the subset comparable to human CS95):\n")
cat(paste0("         ", capture.output(print(ptr, class = FALSE, row.names = FALSE)), collapse = "\n"), "\n")
cat(sprintf("tss_dist NA in %s of %s rows (%s)\n",
            format(n_na_tss, big.mark = ","), format(nrow(final), big.mark = ","),
            paste(sprintf("%s: %d genes", na_genes$species, na_genes$V1), collapse = ", ")))
hr()

# The headline the file exists to support. abs() because tss_dist is signed by genome
# coordinate, not by strand (see the header block).
# ADDITIONAL eQTLs only: sets within COLOC_SAME_BP of a colocalizing position for their own
# gene-tissue are the colocalizing eQTL itself, so counting them in the background compares
# the coloc signal against a background that contains it.
tss_tab <- function(dt) {
    t <- dt[!is.na(tss_dist), {
        q <- unname(quantile(abs_tss_dist, c(0.25, 0.5, 0.75), type = 7))
        .(n_rows = .N, median_abs = q[2], q25 = q[1], q75 = q[3],
          median_signed = as.numeric(stats::median(tss_dist)))
    }, by = .(species, type)]
    setorder(t, species, type)
    t
}
cat("ALL rows (eQTL background still contains the colocalizing eQTL):\n")
print(tss_tab(final), class = FALSE, row.names = FALSE)
cat(sprintf("\nADDITIONAL eQTLs only (sets within %d bp of the coloc position removed):\n",
            COLOC_SAME_BP))
addl <- final[type == "enloc" | !same_as_coloc]
print(tss_tab(addl), class = FALSE, row.names = FALSE)
hr()

# Promoter / enhancer overlap, mirroring get_enr(cutoff = 0): the fraction of positions
# inside an element, over the same null (element bp / genome bp) the enrichment script uses.
nulls <- rbindlist(lapply(names(ELEM), function(sp) data.table(
    species = sp,
    null_pro = ELEM[[sp]]$pro[, sum(e - s)] / GENOME_BP[[sp]],
    null_enh = ELEM[[sp]]$enh[, sum(e - s)] / GENOME_BP[[sp]])))

etab <- addl[!is.na(pos), .(n = .N,
                            pro_frac = mean(in_pro), enh_frac = mean(in_enh)),
             by = .(species, type)]
etab <- merge(etab, nulls, by = "species")
etab[, `:=`(pro_enrich = round(pro_frac / null_pro, 2),
            enh_enrich = round(enh_frac / null_enh, 2),
            pro_frac = round(pro_frac, 4), enh_frac = round(enh_frac, 4))]
setorder(etab, species, type)
cat("PROMOTER / ENHANCER OVERLAP at the summary position (enrichment cutoff 0,\n")
cat("ADDITIONAL eQTLs only)\n")
print(etab[, .(species, type, n, pro_frac, pro_enrich, enh_frac, enh_enrich)],
      class = FALSE, row.names = FALSE)

cat(sprintf("\nSETS ABSORBED AS \"the colocalizing eQTL\" (member within %d bp)\n", COLOC_SAME_BP))
ab <- final[type == "eqtl_median", .(
        sets = .N,
        by_exact_evidence = sum(is_coloc_set),
        by_distance_only  = sum(!is_coloc_set & dist_to_coloc <= COLOC_SAME_BP),
        absorbed = sum(same_as_coloc),
        additional = sum(!same_as_coloc)), by = species]
print(ab, class = FALSE, row.names = FALSE)
cat(sprintf("  median dist_to_coloc among absorbed sets: human %s bp, pig %s bp\n",
            format(final[species=="human" & type=="eqtl_median" & same_as_coloc,
                         stats::median(dist_to_coloc)], big.mark = ","),
            format(final[species=="pig" & type=="eqtl_median" & same_as_coloc,
                         stats::median(dist_to_coloc)], big.mark = ",")))
cat("  gene-tissue pairs left with NO additional eQTL set: ")
cat(sprintf("human %s, pig %s\n",
    format(uniqueN(final[species=="human" & type=="eqtl_median"], by=c("tissue","gene")) -
           uniqueN(final[species=="human" & type=="eqtl_median" & !same_as_coloc], by=c("tissue","gene")), big.mark=","),
    format(uniqueN(final[species=="pig" & type=="eqtl_median"], by=c("tissue","gene")) -
           uniqueN(final[species=="pig" & type=="eqtl_median" & !same_as_coloc], by=c("tissue","gene")), big.mark=",")))

# The member fractions exist because a summary position is synthetic. Show how far the two
# actually diverge rather than asserting that it matters.
mm <- final[type == "eqtl_median" & !is.na(pos)]
cat(sprintf("\nsummary position vs its own members (eqtl_median rows, n = %s):\n",
            format(nrow(mm), big.mark = ",")))
cat(sprintf("  median position in a promoter but NO member is: %s (%.1f%%)\n",
            format(mm[in_pro & frac_mem_pro_0 == 0, .N], big.mark = ","),
            100 * mm[, mean(in_pro & frac_mem_pro_0 == 0)]))
cat(sprintf("  a member in a promoter but the median position is not: %s (%.1f%%)\n",
            format(mm[!in_pro & frac_mem_pro_0 > 0, .N], big.mark = ","),
            100 * mm[, mean(!in_pro & frac_mem_pro_0 > 0)]))
cat(sprintf("  mean fraction of members in a promoter: %.4f | in an enhancer: %.4f\n",
            mm[, mean(frac_mem_pro_0)], mm[, mean(frac_mem_enh_0)]))
hr()
cat("wrote:\n  ", OUT, "\n\n")
