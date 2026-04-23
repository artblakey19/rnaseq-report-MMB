# Per-sample QC metric table.
# Joins samples.tsv metadata, MultiQC general stats, and counts-derived
# assigned_reads. Metrics-only: no pass/fail columns, no thresholds.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "message")
sink(log_con, type = "output")

counts_path  <- snakemake@input[["counts"]]
multiqc_dir  <- snakemake@params[["multiqc_dir"]]
samples_path <- snakemake@input[["samples"]]
out_path     <- snakemake@output[["summary"]]

stopifnot(file.exists(counts_path))
stopifnot(file.exists(samples_path))

REQUIRED_COLS <- c(
  "sample", "condition", "replicate", "batch",
  "total_reads", "mapping_rate", "uniquely_mapped_pct",
  "rrna_pct", "duplication", "dupradar_slope",
  "strandedness", "exonic_pct",
  "assigned_reads"
)

samples <- read_tsv(samples_path, show_col_types = FALSE)
stopifnot("sample" %in% names(samples))
for (col in c("condition", "replicate", "batch")) {
  if (!col %in% names(samples)) samples[[col]] <- NA
}

counts <- read_tsv(counts_path, show_col_types = FALSE)
sample_cols <- setdiff(names(counts), c("gene_id", "gene_name"))

missing_in_counts <- setdiff(samples$sample, sample_cols)
if (length(missing_in_counts) > 0) {
  stop(sprintf(
    "qc_summary: samples in samples.tsv missing from counts matrix: %s",
    paste(missing_in_counts, collapse = ", ")
  ))
}

assigned_tbl <- tibble(
  sample = samples$sample,
  assigned_reads = vapply(
    samples$sample,
    function(s) sum(counts[[s]], na.rm = TRUE),
    numeric(1)
  )
)

# Pull a metric column from a general_stats-style data frame by trying a set
# of alias names. Matches exact or suffix (for module-prefixed MultiQC keys
# such as "STAR_mqc-generalstats-star-uniquely_mapped_percent").
metric_aliases <- list(
  total_reads         = c("total_reads", "total_sequences",
                          "m_reads", "m_seqs", "m_reads_mapped"),
  # Generic "overall % aligned" — from Salmon/STAR aggregate columns. Includes
  # multi-mappers, so it can look healthy even when unique alignments are low.
  mapping_rate        = c("mapping_rate", "pct_mapped", "percent_mapped"),
  # STAR unique alignment rate — the DE-relevant fraction. Big gap between
  # this and mapping_rate flags repeat/pseudo-gene multi-mapping issues.
  uniquely_mapped_pct = c("uniquely_mapped_percent", "uniquely_mapped_pct",
                          "pct_uniquely_mapped"),
  rrna_pct            = c("rrna_pct", "rrna", "percent_rrna", "pct_rrna"),
  duplication         = c("duplication", "percent_duplication",
                          "pct_duplication", "percent_duplicates")
)

normalise_name <- function(x) gsub("[^a-z0-9]", "_", tolower(x))

pick_metric <- function(df, aliases) {
  norm_names <- normalise_name(names(df))
  for (a in aliases) {
    key <- normalise_name(a)
    hit <- which(norm_names == key)
    if (length(hit) > 0) return(df[[hit[1]]])
    hit <- grep(paste0("(^|_)", key, "$"), norm_names)
    if (length(hit) > 0) return(df[[hit[1]]])
  }
  rep(NA_real_, nrow(df))
}

load_multiqc_general_stats <- function(dir) {
  if (is.null(dir) || !nzchar(dir) || !dir.exists(dir)) return(NULL)
  candidates <- c(
    file.path(dir, "multiqc_general_stats.txt"),
    file.path(dir, "multiqc_general_stats.tsv"),
    file.path(dir, "general_stats.tsv"),
    file.path(dir, "general_stats.txt")
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) return(NULL)
  df <- read_tsv(found[1], show_col_types = FALSE)
  if (ncol(df) < 2) return(NULL)
  names(df)[1] <- "sample"
  df$sample <- as.character(df$sample)
  df
}

# Read a per-tool MultiQC TSV (samples-in-col-1). Returns NULL if absent.
load_mqc_module <- function(dir, basename) {
  if (is.null(dir) || !nzchar(dir) || !dir.exists(dir)) return(NULL)
  path <- file.path(dir, basename)
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(read_tsv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df) || ncol(df) < 2) return(NULL)
  names(df)[1] <- "sample"
  df$sample <- as.character(df$sample)
  df
}

# Strandedness from RSeQC infer_experiment. MultiQC emits per-orientation
# fractions (pe_* for paired-end, se_* for single-end); we merge them and
# classify by the dominant category. Also return sense/antisense fractions
# for the report.
parse_infer_experiment <- function(dir) {
  df <- load_mqc_module(dir, "multiqc_rseqc_infer_experiment.txt")
  if (is.null(df)) return(NULL)
  merge_cols <- function(df, suffix) {
    pe <- paste0("pe_", suffix)
    se <- paste0("se_", suffix)
    if (pe %in% names(df) && se %in% names(df)) {
      return(ifelse(is.na(df[[pe]]) | df[[pe]] == 0, df[[se]], df[[pe]]))
    }
    if (pe %in% names(df)) return(df[[pe]])
    if (se %in% names(df)) return(df[[se]])
    rep(NA_real_, nrow(df))
  }
  sense <- suppressWarnings(as.numeric(merge_cols(df, "sense")))
  anti  <- suppressWarnings(as.numeric(merge_cols(df, "antisense")))
  und   <- suppressWarnings(as.numeric(merge_cols(df, "undetermined")))
  # Scale to percent if MultiQC reported fractions (max <=1).
  mx <- suppressWarnings(max(c(sense, anti, und), na.rm = TRUE))
  if (is.finite(mx) && mx <= 1.01) {
    sense <- sense * 100; anti <- anti * 100; und <- und * 100
  }
  call <- mapply(function(s, a, u) {
    vals <- c(forward = s, reverse = a, unstranded = u)
    if (all(is.na(vals))) return(NA_character_)
    names(which.max(vals))
  }, sense, anti, und, USE.NAMES = FALSE)
  tibble(
    sample = df$sample,
    strandedness = call,
    strand_sense_pct = sense,
    strand_antisense_pct = anti
  )
}

# Exonic percentage from RSeQC read_distribution. The *_tag_pct columns are
# the useful per-region numbers; exonic = CDS + 5'UTR + 3'UTR.
parse_read_distribution <- function(dir) {
  df <- load_mqc_module(dir, "multiqc_rseqc_read_distribution.txt")
  if (is.null(df)) return(NULL)
  pct_col <- function(name) {
    hit <- grep(paste0("^", name, "$"), names(df), ignore.case = TRUE, value = TRUE)
    if (length(hit) == 0) return(rep(NA_real_, nrow(df)))
    suppressWarnings(as.numeric(df[[hit[1]]]))
  }
  cds  <- pct_col("cds_exons_tag_pct")
  utr5 <- pct_col("5_utr_exons_tag_pct")
  utr3 <- pct_col("3_utr_exons_tag_pct")
  exonic <- rowSums(cbind(cds, utr5, utr3), na.rm = TRUE)
  exonic[is.na(cds) & is.na(utr5) & is.na(utr3)] <- NA_real_
  tibble(sample = df$sample, exonic_pct = exonic)
}

# dupRadar slope. Low slope => PCR-driven duplication (real complexity loss);
# high slope (~1) => duplication tracks expression (benign for high-expr genes).
parse_dupradar <- function(dir) {
  df <- load_mqc_module(dir, "multiqc_dupradar.txt")
  if (is.null(df)) return(NULL)
  slope_col <- grep("slope", names(df), ignore.case = TRUE, value = TRUE)
  if (length(slope_col) == 0) return(NULL)
  tibble(
    sample = df$sample,
    dupradar_slope = suppressWarnings(as.numeric(df[[slope_col[1]]]))
  )
}

mqc <- load_multiqc_general_stats(multiqc_dir)

if (!is.null(mqc)) {
  mqc_tbl <- tibble(
    sample              = mqc$sample,
    total_reads         = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$total_reads))),
    mapping_rate        = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$mapping_rate))),
    uniquely_mapped_pct = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$uniquely_mapped_pct))),
    rrna_pct            = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$rrna_pct))),
    duplication         = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$duplication)))
  )
  message(sprintf("qc_summary: loaded MultiQC general stats (%d samples)", nrow(mqc_tbl)))
} else {
  message("qc_summary: no MultiQC general_stats file found; QC metrics set to NA")
  mqc_tbl <- tibble(
    sample              = samples$sample,
    total_reads         = NA_real_,
    mapping_rate        = NA_real_,
    uniquely_mapped_pct = NA_real_,
    rrna_pct            = NA_real_,
    duplication         = NA_real_
  )
}

# Per-tool MultiQC modules (optional — absence does not fail the rule).
strand_tbl   <- parse_infer_experiment(multiqc_dir)
readdist_tbl <- parse_read_distribution(multiqc_dir)
dupradar_tbl <- parse_dupradar(multiqc_dir)
for (tbl_name in c("strand_tbl", "readdist_tbl", "dupradar_tbl")) {
  if (!is.null(get(tbl_name))) {
    message(sprintf("qc_summary: loaded %s (%d rows)",
                    tbl_name, nrow(get(tbl_name))))
  }
}

summary_tbl <- samples %>%
  select(sample, condition, replicate, batch) %>%
  left_join(mqc_tbl, by = "sample")
for (extra in list(strand_tbl, readdist_tbl, dupradar_tbl)) {
  if (!is.null(extra)) {
    summary_tbl <- left_join(summary_tbl, extra, by = "sample")
  }
}
summary_tbl <- left_join(summary_tbl, assigned_tbl, by = "sample")

for (col in REQUIRED_COLS) {
  if (!col %in% names(summary_tbl)) summary_tbl[[col]] <- NA
}
summary_tbl <- summary_tbl[, REQUIRED_COLS]

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write_tsv(summary_tbl, out_path)

message(sprintf("qc_summary: wrote %d rows x %d cols -> %s",
                nrow(summary_tbl), ncol(summary_tbl), out_path))
