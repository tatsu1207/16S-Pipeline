#!/usr/bin/env Rscript
# ============================================================================
# MicrobiomeDash — SILVA 138.1 Taxonomy Assignment
#
# Uses DADA2's assignTaxonomy + addSpecies with SILVA reference databases.
# With --species_train (long reads), species are instead predicted by the
# naive Bayesian classifier using the SILVA species-level training set; a
# unique exact match from --silva_species takes precedence when given.
# Outputs: taxonomy.tsv (ASV_ID, Kingdom, Phylum, ..., Species)
# ============================================================================

suppressPackageStartupMessages({
  library(dada2)
  library(optparse)
  library(jsonlite)
})

# ── CLI Arguments ─────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--rep_seqs",      type="character", help="Path to representative sequences FASTA"),
  make_option("--output",        type="character", help="Output taxonomy TSV path"),
  make_option("--silva_train",   type="character", help="SILVA training set path"),
  make_option("--silva_species", type="character", help="SILVA species assignment path"),
  make_option("--threads",       type="integer",   default=1, help="Number of threads"),
  make_option("--species_train", type="character", default=NULL,
              help="SILVA species-level training set (Bayesian species assignment)"),
  make_option("--skip_species",  action="store_true", default=FALSE,
              help="Skip species-level assignment (faster)")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$rep_seqs) || is.null(opt$output) || is.null(opt$silva_train)) {
  stop("--rep_seqs, --output, and --silva_train are required")
}
if (!opt$skip_species && is.null(opt$silva_species) && is.null(opt$species_train)) {
  stop("--silva_species or --species_train is required unless --skip_species is set")
}

# ── Main ──────────────────────────────────────────────────────────────────────

tryCatch({

  # Read representative sequences
  cat("Reading representative sequences...\n")
  lines <- readLines(opt$rep_seqs)
  header_idx <- grep("^>", lines)
  asv_ids <- sub("^>", "", lines[header_idx])
  seqs <- lines[header_idx + 1]
  names(seqs) <- asv_ids

  cat("Loaded", length(seqs), "ASV sequences\n")

  use_mt <- ifelse(opt$threads > 1, opt$threads, FALSE)

  if (!opt$skip_species && !is.null(opt$species_train)) {
    # Kingdom-to-species in one Bayesian pass (bootstrap-supported species)
    cat("Assigning taxonomy with species-level training set (this may take a while)...\n")
    taxa <- assignTaxonomy(seqs, opt$species_train,
                           multithread=use_mt, tryRC=FALSE)
    # A unique 100% match (genus-consistent) overrides the Bayesian species
    if (!is.null(opt$silva_species)) {
      cat("Checking exact species matches...\n")
      exact <- addSpecies(taxa[, 1:6, drop=FALSE], opt$silva_species)
      hit <- !is.na(exact[, "Species"])
      changed <- hit & (is.na(taxa[, "Species"]) | taxa[, "Species"] != exact[, "Species"])
      taxa[hit, "Species"] <- exact[hit, "Species"]
      cat("Exact matches:", sum(hit), "ASVs (", sum(changed), "changed )\n")
    }
  } else {
    # Assign taxonomy to genus level
    cat("Assigning taxonomy (this may take a while)...\n")
    taxa <- assignTaxonomy(seqs, opt$silva_train,
                           multithread=use_mt, tryRC=FALSE)
  }

  # Add species-level assignment by exact matching (single-threaded, slow)
  if (!opt$skip_species && is.null(opt$species_train)) {
    cat("Adding species assignments...\n")
    taxa <- addSpecies(taxa, opt$silva_species)
  } else if (opt$skip_species) {
    cat("Skipping species assignment (--skip_species)\n")
  }

  # Build output data frame
  tax_df <- data.frame(
    ASV_ID  = asv_ids,
    Kingdom = taxa[, "Kingdom"],
    Phylum  = taxa[, "Phylum"],
    Class   = taxa[, "Class"],
    Order   = taxa[, "Order"],
    Family  = taxa[, "Family"],
    Genus   = taxa[, "Genus"],
    Species = if ("Species" %in% colnames(taxa)) taxa[, "Species"] else NA,
    stringsAsFactors = FALSE
  )

  # Write output
  dir.create(dirname(opt$output), recursive=TRUE, showWarnings=FALSE)
  write.table(tax_df, opt$output, sep="\t", row.names=FALSE, quote=FALSE)
  cat("Wrote:", opt$output, "\n")

  # Success
  status <- toJSON(list(
    status    = "success",
    asv_count = length(asv_ids)
  ), auto_unbox=TRUE)
  cat("\n", status, "\n", sep="")
  quit(status=0)

}, error = function(e) {
  status <- toJSON(list(
    status  = "error",
    message = conditionMessage(e)
  ), auto_unbox=TRUE)
  cat("\n", status, "\n", sep="")
  quit(status=1)
})
