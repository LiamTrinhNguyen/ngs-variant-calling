#!/usr/bin/env Rscript
# =============================================================================
# STEP 9b · VARIANT STATISTICS (R) – tests that go beyond counting
# =============================================================================
# WHAT: reads one sample's 06_reports/<sample>.all_calls.tsv (written by
#       filter_and_report() with Polars) and writes two tables:
#         06_reports/<sample>.r_stats.tsv     one row per statistic:
#                                             Metric | Value | Interpretation
#         06_reports/<sample>.r_windows.tsv   trusted variants per genome window,
#                                             with a Poisson test for hotspots
# WHY R: these are classic statistical tests (binomial CI, Wilcoxon,
#        Spearman, Kolmogorov–Smirnov, Poisson + Benjamini–Hochberg) that are
#        one line each in base R. Python (Polars) does the data wrangling,
#        R does the statistics, and dashboard.py shows both on one page.
#
# Only base R is used (no packages), so `r-base` in ngs-env is all you need.
#
# Usage:
#   Rscript R/variant_stats.R <all_calls.tsv> <reference.fasta.fai> <out_dir> <sample> [window_bp]
# Example:
#   Rscript R/variant_stats.R results/06_reports/SRR1553607.all_calls.tsv \
#       results/01_reference/NC_002549.1.fasta.fai results/06_reports SRR1553607 1000
#
# Called by: run_r_stats() in run_interactive_pipeline.py (interactive script),
#            rule r_stats (Snakemake) and process R_STATS (Nextflow).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript variant_stats.R <all_calls.tsv> <reference.fai> <out_dir> <sample> [window_bp]")
}
calls_file <- args[1]
fai_file   <- args[2]
out_dir    <- args[3]
sample     <- args[4]
window_bp  <- if (length(args) >= 5) as.integer(args[5]) else 1000L

# -----------------------------------------------------------------------------
# Read the inputs
#   colClasses for CHROM/REF/ALT: WHY "character": stops R turning a REF of
#   "T" into TRUE, or a chromosome name like "1" into a number.
# -----------------------------------------------------------------------------
calls <- read.delim(calls_file, colClasses = c(CHROM = "character", REF = "character",
                                               ALT = "character", FILTER = "character"),
                    na.strings = c("", "NA", "null"))
fai <- read.delim(fai_file, header = FALSE, colClasses = c("character", "numeric"),
                  col.names = c("name", "length", "offset", "linebases", "linewidth"))
genome_len <- sum(fai$length)

pass    <- calls[calls$FILTER == "PASS", ]
flagged <- calls[calls$FILTER != "PASS", ]

# Collect the results as rows of a long table (Metric, Value, Interpretation).
# WHY long format: the dashboard and Power BI can show any number of
#      statistics without changing their column layout.
rows <- list()
add <- function(metric, value, interpretation) {
  rows[[length(rows) + 1]] <<- data.frame(Metric = metric, Value = value,
                                          Interpretation = interpretation)
}
# fmt: 3 significant digits, no padding. pfmt: "p = 0.012" or "p < 1e-16".
fmt  <- function(x, d = 3) if (is.na(x)) NA_character_ else trimws(formatC(x, digits = d, format = "g"))
pfmt <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 1e-16) "p < 1e-16" else paste("p =", format(signif(p, 2)))
}

# -----------------------------------------------------------------------------
# 1 · Counts and density
# -----------------------------------------------------------------------------
add("Calls (all)", nrow(calls), "Every raw call, PASS or flagged.")
add("Calls (PASS)", nrow(pass), "Calls that passed QUAL, depth and AF filters.")
add("PASS variants per kb", fmt(1000 * nrow(pass) / genome_len),
    sprintf("Over a %s bp genome.", format(genome_len, big.mark = ",")))

# -----------------------------------------------------------------------------
# 2 · Transition/transversion ratio (Ts/Tv) with a 95% confidence interval
#   WHY: Ts/Tv is a standard sanity check. Real mutations are biased towards
#        transitions (A<->G, C<->T), so Ts/Tv is usually well above 0.5;
#        random sequencing errors give about 0.5 (there are twice as many
#        possible transversions). The CI comes from an exact binomial test on
#        the share of transitions, converted to a ratio p / (1 - p).
# -----------------------------------------------------------------------------
snv <- pass[nchar(pass$REF) == 1 & nchar(pass$ALT) == 1, ]
is_ts <- paste(snv$REF, snv$ALT) %in% c("A G", "G A", "C T", "T C")
ts <- sum(is_ts); tv <- sum(!is_ts)
add("SNVs / indels (PASS)", sprintf("%d / %d", nrow(snv), nrow(pass) - nrow(snv)),
    "Single-base changes vs insertions/deletions.")
if (ts + tv > 0 && tv > 0) {
  ci <- binom.test(ts, ts + tv)$conf.int
  add("Ts/Tv ratio (95% CI)",
      sprintf("%s (%s–%s)", fmt(ts / tv), fmt(ci[1] / (1 - ci[1])),
              if (ci[2] < 1) fmt(ci[2] / (1 - ci[2])) else "Inf"),
      if (ci[1] / (1 - ci[1]) > 0.5) "Above 0.5 (random errors): looks like real mutations."
      else "Not clearly above 0.5: could include sequencing errors.")
} else {
  add("Ts/Tv ratio (95% CI)", NA, "Not enough SNVs to estimate.")
}

# -----------------------------------------------------------------------------
# 3 · Do PASS and flagged calls really differ? (Wilcoxon rank-sum tests)
#   WHY: if the filters work, PASS calls should have clearly higher depth and
#        QUAL than flagged ones. Wilcoxon compares medians without assuming
#        normal data (depth and QUAL are heavily skewed).
# -----------------------------------------------------------------------------
compare <- function(col, label) {
  a <- pass[[col]]; b <- flagged[[col]]
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (length(a) >= 2 && length(b) >= 2) {
    p <- suppressWarnings(wilcox.test(a, b)$p.value)   # warnings = ties, expected for depth
    add(sprintf("Median %s: PASS vs flagged", label),
        sprintf("%s vs %s (%s)", fmt(median(a)), fmt(median(b)), pfmt(p)),
        if (p < 0.05) "Significantly different (Wilcoxon): the filter separates the two groups."
        else "Not significantly different: the filter may not be informative here.")
  } else {
    add(sprintf("Median %s: PASS vs flagged", label), NA, "Too few calls in one group to test.")
  }
}
compare("DP", "depth")
compare("QUAL", "QUAL")

# -----------------------------------------------------------------------------
# 4 · Is caller confidence driven by depth? (Spearman correlation)
#   WHY: a strong positive rho means QUAL mostly reflects how many reads there
#        are, so low-depth regions are where calls become unreliable.
# -----------------------------------------------------------------------------
ok <- !is.na(calls$QUAL) & !is.na(calls$DP)
if (sum(ok) >= 3) {
  ct <- suppressWarnings(cor.test(calls$QUAL[ok], calls$DP[ok], method = "spearman"))
  add("Spearman rho, QUAL vs depth (all calls)",
      sprintf("%s (%s)", fmt(unname(ct$estimate)), pfmt(ct$p.value)),
      if (ct$estimate > 0.5) "Strong: confidence rises with depth; low-coverage regions limit calling."
      else "Weak or moderate: QUAL is not just a function of depth.")
}

# -----------------------------------------------------------------------------
# 5 · Allele frequency: is the sample mixed?
#   WHY: in a clean haploid sample nearly every call has AF close to 1.
#        Many calls between 0.2 and 0.8 suggest a mixed infection,
#        contamination or reads from a related genome.
# -----------------------------------------------------------------------------
af <- calls$AF[!is.na(calls$AF)]
if (length(af)) {
  q <- quantile(af, c(0.25, 0.5, 0.75))
  mixed <- mean(af >= 0.2 & af < 0.8)
  add("AF median (IQR), all calls", sprintf("%s (%s–%s)", fmt(q[2]), fmt(q[1]), fmt(q[3])),
      "Share of reads showing the variant; ~1.0 is expected for a haploid genome.")
  add("Calls with AF 0.2–0.8", sprintf("%.1f%%", 100 * mixed),
      if (mixed > 0.1) "More than 10%: check for a mixed or contaminated sample."
      else "Low: no sign of a mixed sample.")
}

# -----------------------------------------------------------------------------
# 6 · Are the variants evenly spread along the genome? (Kolmogorov–Smirnov)
#   WHY: mutations accumulate unevenly between strains, but clustering can also
#        point to a problem region (repeats, mis-mapping, poor coverage).
#        KS compares the positions with a uniform spread over the genome.
#        Contigs are joined end to end, in .fai order.
# -----------------------------------------------------------------------------
offset <- setNames(c(0, cumsum(fai$length))[seq_len(nrow(fai))], fai$name)
gpos <- if (nrow(pass)) offset[pass$CHROM] + pass$POS else numeric(0)
if (length(gpos) >= 5) {
  ks <- suppressWarnings(ks.test(gpos, "punif", 0, genome_len))
  add("Even spread along genome (KS test)",
      sprintf("D = %s (%s)", fmt(unname(ks$statistic)), pfmt(ks$p.value)),
      if (ks$p.value < 0.05) "Not even: variants cluster in some regions (see the hotspot windows)."
      else "Consistent with an even spread.")
}

# -----------------------------------------------------------------------------
# 7 · Hotspot windows (Poisson test + Benjamini–Hochberg correction)
#   WHAT: counts PASS variants in fixed windows (default 1 kb) and asks, for
#         each window: is this count higher than expected if variants were
#         spread evenly? Expected = total variants × window length / genome.
#   WHY BH: many windows are tested at once; the correction keeps false
#           "hotspots" rare (q < 0.05 = less than 5% expected to be false).
# -----------------------------------------------------------------------------
starts <- seq(0, genome_len - 1, by = window_bp)
win <- data.frame(Window_start = starts + 1,
                  Window_end   = pmin(starts + window_bp, genome_len))
win$Variants <- vapply(seq_len(nrow(win)), function(i)
  sum(gpos >= win$Window_start[i] & gpos <= win$Window_end[i]), integer(1))
win$Expected <- round(length(gpos) * (win$Window_end - win$Window_start + 1) / genome_len, 2)
# P(X >= observed) under Poisson(expected)
win$P_value  <- signif(ppois(win$Variants - 1, win$Expected, lower.tail = FALSE), 3)
win$Q_value  <- signif(p.adjust(win$P_value, method = "BH"), 3)
win$Hotspot  <- win$Q_value < 0.05 & win$Variants > win$Expected
hot <- sum(win$Hotspot)
add(sprintf("Hotspot windows (%s bp, BH q < 0.05)", format(window_bp, big.mark = ",")),
    sprintf("%d of %d", hot, nrow(win)),
    if (hot) "Windows with more variants than chance: check genes or coverage there."
    else "No window has more variants than expected by chance.")

# -----------------------------------------------------------------------------
# Write the two tables (tab-separated, like the other 06_reports files)
# -----------------------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stats <- do.call(rbind, rows)
stats_file   <- file.path(out_dir, paste0(sample, ".r_stats.tsv"))
windows_file <- file.path(out_dir, paste0(sample, ".r_windows.tsv"))
write.table(stats, stats_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(win, windows_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("R statistics written -> %s\nHotspot windows written -> %s\n", stats_file, windows_file))
