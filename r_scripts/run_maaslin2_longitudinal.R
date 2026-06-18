#!/usr/bin/env Rscript
# ============================================================================
# 16S Analyzer — MaAsLin2 longitudinal differential abundance
#
# Uses random effects for subject (repeated measures) and time as fixed effect.
# Standard interface: --counts, --metadata, --subject_col, --time_col,
#                     --group_col (optional), --output
# Output: TSV with feature, coef, log2fc, pvalue, qvalue
# ============================================================================

suppressPackageStartupMessages({
  library(Maaslin2)
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  make_option("--counts",      type="character", help="Count matrix TSV"),
  make_option("--metadata",    type="character", help="Metadata TSV"),
  make_option("--subject_col", type="character", help="Subject/individual ID column"),
  make_option("--time_col",    type="character", help="Time variable column"),
  make_option("--group_col",   type="character", default=NULL, help="Optional group fixed effect"),
  make_option("--output",      type="character", help="Output TSV path"),
  make_option("--threads",     type="integer", default=1, help="Number of parallel cores")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Read data
counts <- read.csv(opt$counts, sep="\t", row.names=1, check.names=FALSE)
meta   <- read.csv(opt$metadata, sep="\t", check.names=FALSE)

# Align samples
common <- intersect(colnames(counts), meta$SampleID)
if (length(common) < 3) {
  cat(toJSON(list(status = "error", message = "Too few matching samples"), auto_unbox=TRUE))
  cat("\n")
  quit(status=1)
}

counts <- counts[, common, drop=FALSE]
meta   <- meta[match(common, meta$SampleID), ]
rownames(meta) <- meta$SampleID

# Ensure time column is numeric
meta[[opt$time_col]] <- as.numeric(meta[[opt$time_col]])

# Build fixed and random effects
fixed_fx <- c(opt$time_col)
if (!is.null(opt$group_col) && nchar(opt$group_col) > 0) {
  fixed_fx <- c(fixed_fx, opt$group_col)
}
random_fx <- c(opt$subject_col)

# MaAsLin2 expects samples as rows, features as columns
counts_t <- as.data.frame(t(counts))

# Run MaAsLin2 with random effects
output_dir <- tempdir()
maaslin_dir <- file.path(output_dir, "maaslin2_longitudinal")

n_cores <- max(1L, opt$threads)

tryCatch({
  fit <- Maaslin2(
    input_data      = counts_t,
    input_metadata  = meta,
    output          = maaslin_dir,
    fixed_effects   = fixed_fx,
    random_effects  = random_fx,
    normalization   = "TSS",
    transform       = "LOG",
    analysis_method = "LM",
    min_abundance   = 0.0,
    min_prevalence  = 0.1,
    plot_heatmap    = FALSE,
    plot_scatter    = FALSE,
    cores           = n_cores
  )
}, error = function(e) {
  cat(toJSON(list(status = "error", message = paste("MaAsLin2 failed:", e$message)),
             auto_unbox=TRUE))
  cat("\n")
  quit(status=1)
})

# Read results — filter to time variable coefficient
res <- fit$results
res_time <- res[res$metadata == opt$time_col, ]

if (nrow(res_time) == 0) {
  # Fall back to all results if time filter yields nothing
  res_time <- res
}

results_df <- data.frame(
  feature = res_time$feature,
  coef    = res_time$coef,
  log2fc  = res_time$coef / log(2),
  pvalue  = res_time$pval,
  qvalue  = res_time$qval,
  stringsAsFactors = FALSE
)

results_df <- results_df[order(results_df$qvalue), ]

write.table(results_df, file=opt$output, sep="\t", row.names=FALSE, quote=FALSE)

n_sig <- sum(results_df$qvalue < 0.05, na.rm=TRUE)
cat(toJSON(list(
  status = "success",
  n_features = nrow(results_df),
  n_significant = n_sig
), auto_unbox=TRUE))
cat("\n")
