#!/usr/bin/env Rscript

# Differential splicing analysis using edgeR with cluster-total offsets
#
# This script implements the core statistical model:
#   log(μ_i,s) = X_s * β_i + log(T_C,s)
#
# where:
# - μ_i,s = expected intron count
# - X_s = design matrix (sample groups)
# - β_i = coefficients (intron-specific effects)
# - T_C,s = cluster total offset (provided, not estimated)
#
# Key features:
# - NO library size normalization (norm.factors = 1)
# - Cluster-total offsets handle compositional nature
# - QL GLMs with robust dispersion estimation
# - Tests intron usage proportions, not expression

suppressPackageStartupMessages({
  library(edgeR)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--counts"), type="character", 
              help="Input count matrix (introns x samples)"),
  make_option(c("--offsets"), type="character",
              help="Input offset matrix (log-transformed cluster totals)"),
  make_option(c("--annotations"), type="character",
              help="Intron annotation file with cluster assignments"),
  make_option(c("--samples"), type="character",
              help="Sample metadata file (sample_id, group, [batch])"),
  make_option(c("--output"), type="character",
              help="Output prefix for results files"),
  make_option(c("--output-prefix"), type="character", default=NULL,
              help="Optional prefix for final output filenames"),
  make_option(c("--group_col"), type="character", default="group",
              help="Column name for sample groups [default: group]"),
  make_option(c("--batch_col"), type="character", default=NULL,
              help="Optional column name for batch effects"),
  make_option(c("--contrast"), type="character", default=NULL,
              help="Contrast to test (format: 'GroupA-GroupB'). Required for this script; orchestration of multiple contrasts handled by main pipeline"),
  make_option(c("--fdr_threshold"), type="double", default=0.05,
              help="FDR threshold for significance [default: 0.05]"),
  make_option(c("--min_logFC"), type="double", default=0,
              help="Minimum absolute log2FC for significance [default: 0]")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# Validate required arguments
if (is.null(args$counts) || is.null(args$offsets) || 
    is.null(args$samples) || is.null(args$output)) {
  stop("Missing required arguments. Use --help for usage information.")
}

cat("=== edgeR Differential Splicing Analysis ===\n\n")

# Load data
cat(sprintf("Loading count matrix from: %s\n", args$counts))
counts <- read.table(args$counts, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
cat(sprintf("  %d introns x %d samples\n", nrow(counts), ncol(counts)))

cat(sprintf("Loading offset matrix from: %s\n", args$offsets))
offsets <- read.table(args$offsets, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

cat(sprintf("Loading annotations from: %s\n", args$annotations))
annotations <- read.table(args$annotations, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

cat(sprintf("Loading sample metadata from: %s\n", args$samples))
samples <- read.table(args$samples, header=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char="#", fill=TRUE)

# Verify sample order matches
if (!all(colnames(counts) == colnames(offsets))) {
  stop("Sample order mismatch between counts and offsets")
}

# Match samples metadata to count matrix columns
samples <- samples[match(colnames(counts), samples$sample_id), ]
if (any(is.na(samples$sample_id))) {
  stop("Sample metadata missing for some samples in count matrix")
}

cat(sprintf("Sample groups: %s\n", paste(unique(samples[[args$group_col]]), collapse=", ")))

# Create DGEList
cat("\nCreating DGEList object...\n")
dge <- DGEList(counts=counts, samples=samples)

# CRITICAL: Set normalization factors to 1 (no library size normalization)
dge$samples$norm.factors <- 1
cat("  Normalization factors set to 1 (no lib-size normalization)\n")

# Add cluster-total offsets
cat("  Adding cluster-total offsets...\n")
dge$offset <- as.matrix(offsets)

# Design matrix
cat("\nBuilding design matrix...\n")
if (!is.null(args$batch_col) && args$batch_col %in% colnames(samples)) {
  formula_str <- sprintf("~ 0 + %s + %s", args$group_col, args$batch_col)
  cat(sprintf("  Formula: %s\n", formula_str))
  design <- model.matrix(as.formula(formula_str), data=samples)
} else {
  formula_str <- sprintf("~ 0 + %s", args$group_col)
  cat(sprintf("  Formula: %s\n", formula_str))
  design <- model.matrix(as.formula(formula_str), data=samples)
}

# Clean up design matrix column names
colnames(design) <- gsub(paste0("^", args$group_col), "", colnames(design))
colnames(design) <- make.names(colnames(design))
cat(sprintf("  Design matrix: %d samples x %d coefficients\n", nrow(design), ncol(design)))
cat(sprintf("  Coefficients: %s\n", paste(colnames(design), collapse=", ")))

# Estimate dispersion
cat("\nEstimating dispersions...\n")
dge <- estimateDisp(dge, design, robust=TRUE)
cat(sprintf("  Common dispersion: %.3f\n", dge$common.dispersion))
cat(sprintf("  Trended dispersion range: %.3f - %.3f\n", 
            min(dge$trended.dispersion), max(dge$trended.dispersion)))
cat(sprintf("  BCV range: %.1f%% - %.1f%%\n", 
            sqrt(min(dge$tagwise.dispersion))*100, 
            sqrt(max(dge$tagwise.dispersion))*100))

# Fit QL GLM
cat("\nFitting quasi-likelihood GLM...\n")
fit <- glmQLFit(dge, design, robust=TRUE)
cat(sprintf("  QL dispersion range: %.3f - %.3f\n", 
            min(fit$df.residual.zeros), max(fit$df.residual.zeros)))

# Determine contrast
group_levels <- colnames(design)[!grepl("^batch", colnames(design), ignore.case=TRUE)]

if (!is.null(args$contrast)) {
  # Parse user-specified contrast (e.g., "TDP43-control" or "TDP43-control1,control2")
  contrast_parts <- strsplit(args$contrast, "-")[[1]]
  if (length(contrast_parts) != 2) {
    stop("Contrast must be in format 'GroupA-GroupB' or 'GroupA-GroupB1,GroupB2,...'")
  }
  
  # Group1 is the treatment/numerator
  group1 <- make.names(contrast_parts[1])
  
  # Group2 might be a single group or comma-separated control groups
  group2_raw <- contrast_parts[2]
  group2_list <- strsplit(group2_raw, ",")[[1]]
  group2_list <- sapply(group2_list, function(x) make.names(trimws(x)))
  
  # Validate all groups exist
  if (!group1 %in% group_levels) {
    stop(sprintf("Group '%s' not found in design. Available: %s", 
                 group1, paste(group_levels, collapse=", ")))
  }
  
  missing_controls <- group2_list[!group2_list %in% group_levels]
  if (length(missing_controls) > 0) {
    stop(sprintf("Control group(s) not found in design: %s. Available: %s", 
                 paste(missing_controls, collapse=", "),
                 paste(group_levels, collapse=", ")))
  }
  
  group2 <- group2_list
} else {
  # Default: first vs second group
  if (length(group_levels) < 2) {
    stop("Need at least 2 groups for comparison")
  }
  group1 <- group_levels[1]
  group2 <- group_levels[2]
}

# Create contrast label for output
if (length(group2) == 1) {
  cat(sprintf("\nTesting contrast: %s vs %s\n", group1, group2))
  contrast_label <- sprintf("%s_vs_%s", group1, group2)
} else {
  cat(sprintf("\nTesting contrast: %s vs %s\n", group1, paste(group2, collapse=",")))
  contrast_label <- sprintf("%s_vs_%s", group1, paste(group2, collapse="_"))
}

# Create contrast vector
contrast_vec <- rep(0, ncol(design))
names(contrast_vec) <- colnames(design)
contrast_vec[group1] <- 1

# For multiple control groups, distribute the weight equally
if (length(group2) > 1) {
  for (ctrl_group in group2) {
    contrast_vec[ctrl_group] <- -1 / length(group2)
  }
} else {
  contrast_vec[group2] <- -1
}

# Perform QL F-test
cat("Running quasi-likelihood F-test...\n")
qlf <- glmQLFTest(fit, contrast=contrast_vec)

# Get results
results <- topTags(qlf, n=Inf, sort.by="PValue")$table

# Add intron annotations
results$intron_id <- rownames(results)
results$contrast <- contrast_label

# Add gene_name and intron_status if available in annotations
# Note: donor_cluster and acceptor_cluster are available in the annotations file
# but not included in results since analysis uses shared offsets from both
if ("gene_name" %in% colnames(annotations)) {
  results$gene_name <- annotations[rownames(results), "gene_name"]
}
if ("intron_status" %in% colnames(annotations)) {
  results$intron_status <- annotations[rownames(results), "intron_status"]
}
if ("overlapping_genes" %in% colnames(annotations)) {
  results$overlapping_genes <- annotations[rownames(results), "overlapping_genes"]
}

# Add significance flags
results$significant <- (results$FDR < args$fdr_threshold) & 
                       (abs(results$logFC) >= args$min_logFC)

# Summary statistics
cat("\n=== Results Summary ===\n")
cat(sprintf("Total introns tested: %d\n", nrow(results)))
cat(sprintf("Significant introns (FDR < %.2f, |logFC| >= %.2f): %d (%.1f%%)\n",
            args$fdr_threshold, args$min_logFC, 
            sum(results$significant),
            100 * mean(results$significant)))

sig_up <- sum(results$significant & results$logFC > 0)
sig_down <- sum(results$significant & results$logFC < 0)
cat(sprintf("  Increased usage (%s > %s): %d\n", group1, group2, sig_up))
cat(sprintf("  Decreased usage (%s < %s): %d\n", group1, group2, sig_down))

# LogFC distribution
cat(sprintf("\nLogFC distribution (all introns):\n"))
cat(sprintf("  Range: %.2f to %.2f\n", min(results$logFC), max(results$logFC)))
cat(sprintf("  Median: %.2f\n", median(results$logFC)))

# Determine output prefix for files
if (!is.null(args$output_prefix)) {
  final_prefix <- args$output_prefix
} else {
  final_prefix <- args$output
}

# Write results
output_file <- paste0(final_prefix, ".intron_results.tsv")
cat(sprintf("\nWriting results to: %s\n", output_file))
write.table(results, output_file, sep="\t", quote=FALSE, row.names=FALSE)

# Write significant hits only
sig_file <- paste0(final_prefix, ".significant_introns.tsv")
sig_results <- results[results$significant, ]
if (nrow(sig_results) > 0) {
  write.table(sig_results, sig_file, sep="\t", quote=FALSE, row.names=FALSE)
  cat(sprintf("Significant introns written to: %s\n", sig_file))
}

# Save R objects for downstream analysis
rdata_file <- paste0(final_prefix, ".RData")
save(dge, fit, qlf, results, file=rdata_file)
cat(sprintf("R objects saved to: %s\n", rdata_file))

# Create diagnostic plots
pdf_file <- paste0(final_prefix, ".diagnostics.pdf")
cat(sprintf("\nGenerating diagnostic plots: %s\n", pdf_file))
pdf(pdf_file, width=10, height=8)

# BCV plot
plotBCV(dge, main="Biological Coefficient of Variation")

# QL dispersion plot  
plotQLDisp(fit, main="Quasi-Likelihood Dispersion")

# MA plot
plotMD(qlf, main=sprintf("MA Plot: %s vs %s", group1, group2))
abline(h=c(-1, 1), col="blue", lty=2)

# Volcano plot
plot(results$logFC, -log10(results$PValue),
     pch=20, cex=0.5, col=ifelse(results$significant, "red", "gray"),
     xlab="log2 Fold Change (Intron Usage)",
     ylab="-log10(P-value)",
     main=sprintf("Volcano Plot: %s vs %s", group1, group2))
abline(v=c(-args$min_logFC, args$min_logFC), col="blue", lty=2)
abline(h=-log10(args$fdr_threshold), col="blue", lty=2)
legend("topright", legend=c("Significant", "Not significant"),
       col=c("red", "gray"), pch=20)

dev.off()

cat("\n=== Analysis Complete! ===\n")
