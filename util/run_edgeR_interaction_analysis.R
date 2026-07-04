#!/usr/bin/env Rscript

# DEXSeq-style stacked NB interaction analysis for differential intron usage.
#
# Model (per intron i, doubled pseudo-sample design):
#
#   log(E[Y_i,s,t]) = sample_s + exon_type_t + group_s:exon_type_t  [+ batch_s:exon_type_t]
#
# where t in {focal, other}, Y_i,s,focal = intron_count, Y_i,s,other = denominator - count.
# The sample effect absorbs per-sample depth variation without a fixed offset.
# The group:exon_type interaction is tested by LRT.
#
# logFC: log2 of the focal/other odds-ratio change between groups.
#   Positive = focal count increases relative to other in the treatment group.
#   NOTE: this is a log-odds-ratio, not a PSI log-ratio (differs from offset mode).
#
# PSI and delta_PSI are computed upstream by run_diff_splice_analysis.py and
# added to the output there; this script reports the interaction logFC and FDR.

suppressPackageStartupMessages({
  library(edgeR)
  library(optparse)
})

option_list <- list(
  make_option(c("--counts"),       type="character", help="Filtered intron count matrix (introns x samples)"),
  make_option(c("--denominators"), type="character", help="Filtered raw denominator matrix (introns x samples)"),
  make_option(c("--annotations"),  type="character", help="Filtered intron annotation file"),
  make_option(c("--samples"),      type="character", help="Sample metadata (sample_id, group, [batch])"),
  make_option(c("--output"),       type="character", help="Output prefix"),
  make_option(c("--group_col"),    type="character", default="group",  help="Group column name [default: group]"),
  make_option(c("--batch_col"),    type="character", default=NULL,     help="Optional batch column name"),
  make_option(c("--contrast"),     type="character", default=NULL,     help="GroupA,GroupB (logFC = GroupA/GroupB)"),
  make_option(c("--fdr_threshold"),type="double",    default=0.05,     help="FDR threshold [default: 0.05]"),
  make_option(c("--min_logFC"),    type="double",    default=0,        help="Minimum |log2FC| [default: 0]")
)

parser <- OptionParser(option_list=option_list)
args   <- parse_args(parser)

if (is.null(args$counts) || is.null(args$denominators) || is.null(args$annotations) ||
    is.null(args$samples) || is.null(args$output) || is.null(args$contrast)) {
  stop("Missing required arguments. Use --help for usage.")
}

cat("=== edgeR DEXSeq-style Interaction Analysis ===\n\n")
cat("Model: ~ sample_id + exon_type + group:exon_type\n")
cat("Test:  LRT on group:exon_type interaction\n\n")

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
# Parse contrast
# ---------------------------------------------------------------------------
contrast_parts <- strsplit(args$contrast, ",")[[1]]
if (length(contrast_parts) != 2) {
  stop("Contrast must be GroupA,GroupB")
}
group1 <- trimws(contrast_parts[1])   # treatment / numerator
group2 <- trimws(contrast_parts[2])   # reference  / denominator

groups <- as.character(samples[[args$group_col]])
if (!group1 %in% groups) stop(sprintf("Group '%s' not found in metadata", group1))
if (!group2 %in% groups) stop(sprintf("Group '%s' not found in metadata", group2))

# Subset to contrast groups only
use  <- groups %in% c(group1, group2)
counts_sub  <- as.matrix(counts[,  use, drop=FALSE])
denoms_sub  <- as.matrix(denoms[,  use, drop=FALSE])
samples_sub <- samples[use, ]
groups_sub  <- groups[use]
N <- ncol(counts_sub)

cat(sprintf("\nContrast: %s vs %s  (%d samples)\n", group1, group2, N))

# ---------------------------------------------------------------------------
# Build other-count matrix: clamp to 0 so denominators never go negative
# ---------------------------------------------------------------------------
other_sub <- pmax(0L, denoms_sub - counts_sub)

n_neg <- sum(denoms_sub < counts_sub)
if (n_neg > 0) {
  cat(sprintf("  NOTE: %d intron/sample cells had denominator < count; clamped to 0\n", n_neg))
}

# ---------------------------------------------------------------------------
# Build doubled count matrix: n_introns x 2N
#   Columns 1..N  : focal (Y_i,s)
#   Columns N+1..2N: other (D_i,s - Y_i,s)
# ---------------------------------------------------------------------------
stacked <- cbind(counts_sub, other_sub)
colnames(stacked) <- c(
  paste0(colnames(counts_sub), "__focal"),
  paste0(colnames(counts_sub), "__other")
)

# ---------------------------------------------------------------------------
# Build doubled sample metadata: 2N rows
# exon_type "other" is reference level so the focal:group coefficient tests
# the change in focal/other ratio for the treatment vs reference.
# ---------------------------------------------------------------------------
exon_type  <- factor(c(rep("focal", N), rep("other", N)), levels = c("other", "focal"))
sample_id  <- factor(rep(samples_sub$sample_id, 2))
group_fac  <- factor(rep(groups_sub, 2), levels = c(group2, group1))

sample_meta <- data.frame(
  sample_id = sample_id,
  group     = group_fac,
  exon_type = exon_type,
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Design matrices (full and reduced)
# ---------------------------------------------------------------------------
has_batch <- !is.null(args$batch_col) && args$batch_col %in% colnames(samples_sub)

if (has_batch) {
  batch_vec  <- factor(rep(as.character(samples_sub[[args$batch_col]]), 2))
  sample_meta$batch <- batch_vec
  cat(sprintf("  Batch correction: ~ sample_id + exon_type + batch:exon_type + group:exon_type\n"))
  design_full    <- model.matrix(~ sample_id + exon_type + batch:exon_type + group:exon_type,
                                 data = sample_meta)
  design_reduced <- model.matrix(~ sample_id + exon_type + batch:exon_type,
                                 data = sample_meta)
} else {
  cat(sprintf("  Design: ~ sample_id + exon_type + group:exon_type\n"))
  design_full    <- model.matrix(~ sample_id + exon_type + group:exon_type,
                                 data = sample_meta)
  design_reduced <- model.matrix(~ sample_id + exon_type,
                                 data = sample_meta)
}

# Drop aliased columns (safety check)
qr_full <- qr(design_full)
if (qr_full$rank < ncol(design_full)) {
  keep <- qr_full$pivot[seq_len(qr_full$rank)]
  cat(sprintf("  Dropped %d aliased design columns\n", ncol(design_full) - length(keep)))
  design_full <- design_full[, keep, drop=FALSE]
}

# Locate interaction coefficient: group<treatment>:exon_typefocal
interaction_col <- grep(
  paste0("group", make.names(group1), ":exon_typefocal"),
  colnames(design_full)
)
if (length(interaction_col) == 0) {
  # Fallback: any column with both group1 and "focal"
  interaction_col <- intersect(
    grep(make.names(group1), colnames(design_full)),
    grep("exon_typefocal",   colnames(design_full))
  )
}
if (length(interaction_col) != 1) {
  stop(sprintf(
    "Could not identify a unique interaction coefficient. Candidates: %s",
    paste(colnames(design_full)[interaction_col], collapse=", ")
  ))
}
cat(sprintf("  Interaction coefficient: %s (col %d)\n\n",
            colnames(design_full)[interaction_col], interaction_col))

# ---------------------------------------------------------------------------
# edgeR: fit and test
# ---------------------------------------------------------------------------
cat("Creating DGEList (doubled matrix)...\n")
cat(sprintf("  Dimensions: %d introns x %d pseudo-samples\n", nrow(stacked), ncol(stacked)))

dge <- DGEList(counts = stacked)
dge$samples$norm.factors <- 1   # no TMM; sample_id effects absorb depth

cat("Estimating dispersions...\n")
dge <- estimateDisp(dge, design = design_full, robust = TRUE)
cat(sprintf("  Common BCV: %.1f%%\n", sqrt(dge$common.dispersion) * 100))

cat("Fitting NB GLM...\n")
fit <- glmFit(dge, design = design_full)

cat("Running LRT for interaction term...\n")
lrt <- glmLRT(fit, coef = interaction_col)

# ---------------------------------------------------------------------------
# Collect results
# ---------------------------------------------------------------------------
results <- topTags(lrt, n = Inf, sort.by = "PValue")$table

results$intron_id <- rownames(results)
contrast_label    <- sprintf("%s_vs_%s", make.names(group1), make.names(group2))
results$contrast  <- contrast_label
results$stat_mode <- "interaction"

# Annotations
annotation_cols <- c(
  "chr", "start", "end", "strand", "donor", "acceptor",
  "splice_pair", "splice_flag",
  "donor_cluster", "acceptor_cluster", "donor_cluster_size", "acceptor_cluster_size",
  "both_splice_sites_singleton", "offset_mode", "offset_source", "site_depth_fallback_used",
  "gene_name", "intron_status", "overlapping_genes"
)
for (col in annotation_cols) {
  if (col %in% colnames(annotations)) {
    results[[col]] <- annotations[rownames(results), col]
  }
}

# Significance
results$significant <- (results$FDR < args$fdr_threshold) & (abs(results$logFC) >= args$min_logFC)

# Per-group mean focal counts (for interpretability alongside PSI)
g1_focal_cols <- paste0(samples_sub$sample_id[groups_sub == group1], "__focal")
g2_focal_cols <- paste0(samples_sub$sample_id[groups_sub == group2], "__focal")
results$contrast_group1_mean_count <- rowMeans(stacked[rownames(results), g1_focal_cols, drop=FALSE])
results$contrast_group2_mean_count <- rowMeans(stacked[rownames(results), g2_focal_cols, drop=FALSE])

# logCPM from focal pseudo-samples only (interpretable as intron expression level)
focal_cols     <- paste0(colnames(counts_sub), "__focal")
logcpm_focal   <- cpm(stacked[, focal_cols, drop=FALSE], log=TRUE)
results$contrast_group1_mean_logCPM <- rowMeans(logcpm_focal[rownames(results), g1_focal_cols, drop=FALSE])
results$contrast_group2_mean_logCPM <- rowMeans(logcpm_focal[rownames(results), g2_focal_cols, drop=FALSE])

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("\n=== Interaction Analysis Results Summary ===\n")
cat(sprintf("Total introns tested: %d\n", nrow(results)))
cat(sprintf("Significant (FDR < %.2f, |logFC| >= %.2f): %d (%.1f%%)\n",
            args$fdr_threshold, args$min_logFC,
            sum(results$significant), 100 * mean(results$significant)))
cat(sprintf("  Increased focal usage (%s > %s): %d\n",
            group1, group2, sum(results$significant & results$logFC > 0)))
cat(sprintf("  Decreased focal usage (%s < %s): %d\n",
            group1, group2, sum(results$significant & results$logFC < 0)))
cat(sprintf("\nlogFC (log2 odds-ratio) range: %.2f to %.2f  median: %.2f\n",
            min(results$logFC), max(results$logFC), median(results$logFC)))

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

rdata_file <- paste0(args$output, ".interaction.RData")
save(dge, fit, lrt, results, file=rdata_file)
cat(sprintf("R objects saved to:  %s\n", rdata_file))
