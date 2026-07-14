#!/usr/bin/env Rscript

# DEXSeq interaction analysis for differential intron usage.
#
# This runner uses DSF's filtered focal counts and selected denominators:
#   focal = intron count
#   other = max(0, selected denominator - intron count)
#
# DEXSeq is then run in inclusion/exclusion mode via alternativeCountData.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--counts"),       type="character", help="Filtered intron count matrix (introns x samples)"),
  make_option(c("--denominators"), type="character", help="Filtered raw denominator matrix (introns x samples)"),
  make_option(c("--annotations"),  type="character", help="Filtered intron annotation file"),
  make_option(c("--samples"),      type="character", help="Sample metadata (sample_id, group, [batch])"),
  make_option(c("--output"),       type="character", help="Output prefix"),
  make_option(c("--group_col"),    type="character", default="group", help="Group column name [default: group]"),
  make_option(c("--batch_col"),    type="character", default=NULL, help="Optional batch column name"),
  make_option(c("--contrast"),     type="character", default=NULL, help="GroupA,GroupB (logFC = GroupA/GroupB)"),
  make_option(c("--fdr_threshold"),type="double", default=0.05, help="FDR threshold [default: 0.05]"),
  make_option(c("--min_logFC"),    type="double", default=0, help="Minimum |log2FC| [default: 0]")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

if (is.null(args$counts) || is.null(args$denominators) || is.null(args$annotations) ||
    is.null(args$samples) || is.null(args$output) || is.null(args$contrast)) {
  stop("Missing required arguments. Use --help for usage.")
}

if (!requireNamespace("DEXSeq", quietly=TRUE)) {
  stop(
    "The DEXSeq package is required for --intx_engine DEXSeq. ",
    "Install it with: BiocManager::install(\"DEXSeq\")",
    call.=FALSE
  )
}
suppressPackageStartupMessages(library(DEXSeq))

as_integer_counts <- function(x, label) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  if (any(is.na(x))) {
    stop(sprintf("%s contains NA values", label))
  }
  if (any(x < 0)) {
    stop(sprintf("%s contains negative values", label))
  }
  if (any(abs(x - round(x)) > 1e-6)) {
    stop(sprintf("%s contains non-integer values; DEXSeq requires count-like integer inputs", label))
  }
  rounded <- round(x)
  storage.mode(rounded) <- "integer"
  rounded
}

find_logfc_column <- function(results_df, group1, group2) {
  exact <- paste0("log2fold_", group1, "_", group2)
  if (exact %in% colnames(results_df)) {
    return(exact)
  }

  logfc_cols <- grep("^log2fold_", colnames(results_df), value=TRUE)
  candidates <- logfc_cols[grepl(group1, logfc_cols, fixed=TRUE) & grepl(group2, logfc_cols, fixed=TRUE)]
  if (length(candidates) == 1) {
    return(candidates)
  }
  if (length(logfc_cols) == 1) {
    return(logfc_cols)
  }

  stop(sprintf(
    "Could not identify a unique DEXSeq log2 fold-change column for %s vs %s. Available log2fold columns: %s",
    group1, group2, paste(logfc_cols, collapse=", ")
  ))
}

cat("=== DEXSeq Interaction Analysis ===\n\n")
cat("Model: ~ sample + exon + condition:exon [+ batch:exon]\n")
cat("Test:  LRT on condition:exon interaction\n\n")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cat(sprintf("Loading count matrix:       %s\n", args$counts))
counts <- read.table(args$counts, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
cat(sprintf("  %d introns x %d samples\n", nrow(counts), ncol(counts)))

cat(sprintf("Loading denominator matrix: %s\n", args$denominators))
denoms <- read.table(args$denominators, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

cat(sprintf("Loading annotations:        %s\n", args$annotations))
annotations <- read.table(args$annotations, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

cat(sprintf("Loading sample metadata:    %s\n", args$samples))
samples <- read.table(args$samples, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                      comment.char="#", fill=TRUE)
samples <- samples[match(colnames(counts), samples$sample_id), ]
if (any(is.na(samples$sample_id))) {
  stop("Sample metadata missing for some samples in count matrix")
}

if (!all(colnames(counts) == colnames(denoms))) {
  stop("Sample order mismatch between counts and denominators")
}
if (!all(rownames(counts) == rownames(denoms))) {
  stop("Intron order mismatch between counts and denominators")
}

# ---------------------------------------------------------------------------
# Parse contrast and subset samples
# ---------------------------------------------------------------------------
contrast_parts <- strsplit(args$contrast, ",")[[1]]
if (length(contrast_parts) != 2) {
  stop("Contrast must be GroupA,GroupB")
}
group1_raw <- trimws(contrast_parts[1])
group2_raw <- trimws(contrast_parts[2])
group1 <- make.names(group1_raw)
group2 <- make.names(group2_raw)

groups_raw <- as.character(samples[[args$group_col]])
if (!group1_raw %in% groups_raw) stop(sprintf("Group '%s' not found in metadata", group1_raw))
if (!group2_raw %in% groups_raw) stop(sprintf("Group '%s' not found in metadata", group2_raw))

use <- groups_raw %in% c(group1_raw, group2_raw)
counts_sub <- as.matrix(counts[, use, drop=FALSE])
denoms_sub <- as.matrix(denoms[, use, drop=FALSE])
samples_sub <- samples[use, ]
groups_sub_raw <- groups_raw[use]

cat(sprintf("\nContrast: %s vs %s  (%d samples)\n", group1_raw, group2_raw, ncol(counts_sub)))

# ---------------------------------------------------------------------------
# Focal and alternative count matrices
# ---------------------------------------------------------------------------
other_sub <- denoms_sub - counts_sub
other_sub[other_sub < 0] <- 0

n_neg <- sum(denoms_sub < counts_sub)
if (n_neg > 0) {
  cat(sprintf("  NOTE: %d intron/sample cells had denominator < count; clamped to 0\n", n_neg))
}

counts_int <- as_integer_counts(counts_sub, "counts")
other_int <- as_integer_counts(other_sub, "alternative counts")

# ---------------------------------------------------------------------------
# DEXSeq sample metadata and model
# ---------------------------------------------------------------------------
dex_samples <- data.frame(
  condition = factor(make.names(groups_sub_raw), levels=c(group2, group1)),
  stringsAsFactors = FALSE
)
rownames(dex_samples) <- samples_sub$sample_id

has_batch <- !is.null(args$batch_col) && args$batch_col %in% colnames(samples_sub)
if (has_batch) {
  dex_samples$batch <- factor(make.names(as.character(samples_sub[[args$batch_col]])))
  full_model <- ~ sample + exon + batch:exon + condition:exon
  reduced_model <- ~ sample + exon + batch:exon
  cat("  Batch correction: ~ sample + exon + batch:exon + condition:exon\n")
} else {
  full_model <- ~ sample + exon + condition:exon
  reduced_model <- ~ sample + exon
  cat("  Design: ~ sample + exon + condition:exon\n")
}

cat("Creating DEXSeqDataSet...\n")
dxd <- DEXSeqDataSet(
  countData = counts_int,
  sampleData = dex_samples,
  design = full_model,
  featureID = rownames(counts_int),
  groupID = rownames(counts_int),
  alternativeCountData = other_int
)

cat("Estimating size factors...\n")
dxd <- estimateSizeFactors(dxd)

cat("Estimating dispersions...\n")
dxd <- estimateDispersions(dxd, formula=full_model, quiet=FALSE)

cat("Running LRT for condition:exon interaction...\n")
dxd <- testForDEU(dxd, fullModel=full_model, reducedModel=reduced_model)

cat("Estimating exon fold changes...\n")
dxd <- estimateExonFoldChanges(
  dxd,
  fitExpToVar = "condition",
  denominator = group2,
  independentFiltering = FALSE
)

dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
dxr_df <- as.data.frame(dxr)

logfc_col <- find_logfc_column(dxr_df, group1, group2)
cat(sprintf("  DEXSeq logFC column: %s\n", logfc_col))

# ---------------------------------------------------------------------------
# Collect DSF-compatible results
# ---------------------------------------------------------------------------
results <- data.frame(
  intron_id = as.character(dxr_df$featureID),
  logFC = dxr_df[[logfc_col]],
  logCPM = dxr_df$exonBaseMean,
  LR = dxr_df$stat,
  PValue = dxr_df$pvalue,
  FDR = dxr_df$padj,
  stringsAsFactors = FALSE
)

contrast_label <- sprintf("%s_vs_%s", group1, group2)
results$contrast <- contrast_label
results$stat_mode <- "interaction"
results$intx_engine <- "DEXSeq"

annotation_cols <- c(
  "chr", "start", "end", "strand", "donor", "acceptor",
  "splice_pair", "splice_flag",
  "donor_cluster", "acceptor_cluster", "donor_cluster_size", "acceptor_cluster_size",
  "both_splice_sites_singleton", "offset_mode", "offset_source", "site_depth_fallback_used",
  "gene_name", "intron_status", "overlapping_genes"
)
for (col in annotation_cols) {
  if (col %in% colnames(annotations)) {
    results[[col]] <- annotations[results$intron_id, col]
  }
}

results$significant <- !is.na(results$FDR) &
  (results$FDR < args$fdr_threshold) &
  (abs(results$logFC) >= args$min_logFC)

g1_samples <- samples_sub$sample_id[groups_sub_raw == group1_raw]
g2_samples <- samples_sub$sample_id[groups_sub_raw == group2_raw]
results$contrast_group1_mean_count <- rowMeans(counts_sub[results$intron_id, g1_samples, drop=FALSE])
results$contrast_group2_mean_count <- rowMeans(counts_sub[results$intron_id, g2_samples, drop=FALSE])

logcpm_focal <- edgeR::cpm(counts_sub, log=TRUE)
results$contrast_group1_mean_logCPM <- rowMeans(logcpm_focal[results$intron_id, g1_samples, drop=FALSE])
results$contrast_group2_mean_logCPM <- rowMeans(logcpm_focal[results$intron_id, g2_samples, drop=FALSE])

results <- results[order(is.na(results$PValue), results$PValue), ]

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("\n=== DEXSeq Interaction Analysis Results Summary ===\n")
cat(sprintf("Total introns tested: %d\n", nrow(results)))
cat(sprintf("Significant (FDR < %.2f, |logFC| >= %.2f): %d (%.1f%%)\n",
            args$fdr_threshold, args$min_logFC,
            sum(results$significant), 100 * mean(results$significant)))
if (any(!is.na(results$logFC))) {
  cat(sprintf("\nlogFC range: %.2f to %.2f  median: %.2f\n",
              min(results$logFC, na.rm=TRUE),
              max(results$logFC, na.rm=TRUE),
              median(results$logFC, na.rm=TRUE)))
}

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------
output_file <- paste0(args$output, ".intron_results.tsv")
cat(sprintf("\nWriting results to: %s\n", output_file))
write.table(results, output_file, sep="\t", quote=FALSE, row.names=FALSE)

sig_file <- paste0(args$output, ".significant_introns.tsv")
sig_results <- results[results$significant, , drop=FALSE]
if (nrow(sig_results) > 0) {
  write.table(sig_results, sig_file, sep="\t", quote=FALSE, row.names=FALSE)
  cat(sprintf("Significant introns: %s\n", sig_file))
}

rdata_file <- paste0(args$output, ".DEXSeq_interaction.RData")
save(dxd, dxr, results, file=rdata_file)
cat(sprintf("R objects saved to:  %s\n", rdata_file))
