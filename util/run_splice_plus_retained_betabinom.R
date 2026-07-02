#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--counts"), type="character", help="Filtered intron count matrix"),
  make_option(c("--denominators"), type="character", help="Filtered raw trial denominator matrix"),
  make_option(c("--annotations"), type="character", help="Filtered intron annotations"),
  make_option(c("--samples"), type="character", help="Sample metadata file"),
  make_option(c("--output"), type="character", help="Output prefix"),
  make_option(c("--group_col"), type="character", default="group", help="Group column"),
  make_option(c("--batch_col"), type="character", default=NULL, help="Optional batch column"),
  make_option(c("--contrast"), type="character", default=NULL, help="GroupA,GroupB"),
  make_option(c("--fdr_threshold"), type="double", default=0.05, help="FDR threshold"),
  make_option(c("--min_logFC"), type="double", default=0, help="Minimum absolute log2 odds ratio")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

if (is.null(args$counts) || is.null(args$denominators) || is.null(args$annotations) ||
    is.null(args$samples) || is.null(args$output) || is.null(args$contrast)) {
  stop("Missing required arguments. Use --help for usage information.")
}

cat("=== Splice-plus-retained beta-binomial-style analysis ===\n\n")
cat("Model engine: base R quasibinomial GLM over focal/rest counts\n")

counts <- read.table(args$counts, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
denoms <- read.table(args$denominators, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
annotations <- read.table(args$annotations, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
samples <- read.table(args$samples, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                      comment.char="#", fill=TRUE)

if (!all(colnames(counts) == colnames(denoms))) {
  stop("Sample order mismatch between counts and denominators")
}
if (!all(rownames(counts) == rownames(denoms))) {
  stop("Intron order mismatch between counts and denominators")
}

samples <- samples[match(colnames(counts), samples$sample_id), ]
if (any(is.na(samples$sample_id))) {
  stop("Sample metadata missing for some samples in count matrix")
}

contrast_parts <- strsplit(args$contrast, ",")[[1]]
if (length(contrast_parts) != 2) {
  stop("Contrast must be in format 'GroupA,GroupB'")
}
group1 <- trimws(contrast_parts[1])
group2 <- trimws(contrast_parts[2])
if (grepl(";", group1, fixed=TRUE) || grepl(";", group2, fixed=TRUE)) {
  stop("Semicolon-pooled controls are not supported by splice_plus_retained_betabinom")
}
groups <- as.character(samples[[args$group_col]])
if (!all(c(group1, group2) %in% groups)) {
  stop(sprintf("Contrast groups not found in metadata. Requested: %s,%s; available: %s",
               group1, group2, paste(unique(groups), collapse=", ")))
}

counts_mat <- as.matrix(counts)
denom_mat <- as.matrix(denoms)
failures <- denom_mat - counts_mat
if (any(failures < 0, na.rm=TRUE)) {
  bad <- which(failures < 0, arr.ind=TRUE)
  examples <- apply(head(bad, 5), 1, function(idx) {
    sprintf("%s/%s", rownames(failures)[idx[1]], colnames(failures)[idx[2]])
  })
  stop(sprintf("Negative beta-binomial failures: denominator < count for %d intron/sample values. Examples: %s",
               nrow(bad), paste(examples, collapse=", ")))
}

sample_data <- samples
sample_data$group_factor <- factor(groups, levels=c(group2, group1))
if (any(is.na(sample_data$group_factor))) {
  keep_contrast_samples <- !is.na(sample_data$group_factor)
} else {
  keep_contrast_samples <- rep(TRUE, nrow(sample_data))
}

if (!is.null(args$batch_col) && args$batch_col %in% colnames(sample_data)) {
  formula_alt <- as.formula(sprintf("cbind(success, failure) ~ group_factor + %s", args$batch_col))
} else {
  formula_alt <- cbind(success, failure) ~ group_factor
}

safe_logcpm <- function(mat) {
  lib <- colSums(mat)
  sweep(mat + 0.5, 2, lib + 1, "/") * 1e6
}
logcpm <- log2(safe_logcpm(counts_mat))

results <- vector("list", nrow(counts_mat))
names(results) <- rownames(counts_mat)

for (i in seq_len(nrow(counts_mat))) {
  success <- as.numeric(counts_mat[i, ])
  failure <- as.numeric(failures[i, ])
  total <- success + failure
  keep <- keep_contrast_samples & total > 0
  row_out <- list(
    intron_id=rownames(counts_mat)[i],
    logFC=NA_real_,
    logCPM=mean(logcpm[i, ], na.rm=TRUE),
    F=NA_real_,
    PValue=NA_real_
  )

  if (sum(keep) >= 3 && length(unique(sample_data$group_factor[keep])) == 2) {
    df <- sample_data[keep, , drop=FALSE]
    df$success <- success[keep]
    df$failure <- failure[keep]
    fit <- tryCatch(glm(formula_alt, family=quasibinomial(), data=df),
                    error=function(e) NULL)
    if (!is.null(fit)) {
      co <- tryCatch(summary(fit)$coefficients, error=function(e) NULL)
      coef_name <- paste0("group_factor", group1)
      if (!is.null(co) && coef_name %in% rownames(co)) {
        estimate <- co[coef_name, "Estimate"]
        stat <- co[coef_name, "t value"]
        pval <- co[coef_name, "Pr(>|t|)"]
        row_out$logFC <- estimate / log(2)
        row_out$F <- stat^2
        row_out$PValue <- pval
      }
    }
  }

  results[[i]] <- as.data.frame(row_out, stringsAsFactors=FALSE)
}

res <- do.call(rbind, results)
res$contrast <- sprintf("%s_vs_%s", make.names(group1), make.names(group2))
res$stat_engine <- "quasibinomial"
res$FDR <- p.adjust(res$PValue, "BH")

annotation_cols_to_add <- c(
  "chr", "start", "end", "strand", "donor", "acceptor",
  "splice_pair", "splice_flag", "donor_cluster", "acceptor_cluster",
  "donor_cluster_size", "acceptor_cluster_size", "both_splice_sites_singleton",
  "offset_mode", "offset_source", "site_depth_fallback_used",
  "gene_name", "intron_status", "overlapping_genes"
)
for (annotation_col in annotation_cols_to_add) {
  if (annotation_col %in% colnames(annotations)) {
    res[[annotation_col]] <- annotations[res$intron_id, annotation_col]
  }
}

group1_samples <- samples$sample_id[samples[[args$group_col]] == group1]
group2_samples <- samples$sample_id[samples[[args$group_col]] == group2]
res$contrast_group1_mean_count <- rowMeans(counts_mat[res$intron_id, group1_samples, drop=FALSE])
res$contrast_group2_mean_count <- rowMeans(counts_mat[res$intron_id, group2_samples, drop=FALSE])
res$contrast_group1_mean_logCPM <- rowMeans(logcpm[res$intron_id, group1_samples, drop=FALSE])
res$contrast_group2_mean_logCPM <- rowMeans(logcpm[res$intron_id, group2_samples, drop=FALSE])
res$significant <- !is.na(res$FDR) & (res$FDR <= args$fdr_threshold) &
                   !is.na(res$logFC) & (abs(res$logFC) >= args$min_logFC)

res <- res[order(res$PValue, na.last=TRUE), ]

cat(sprintf("Total introns tested: %d\n", nrow(res)))
cat(sprintf("Fitted introns with p-values: %d\n", sum(!is.na(res$PValue))))
cat(sprintf("Significant introns (FDR <= %.2f, |logFC| >= %.2f): %d\n",
            args$fdr_threshold, args$min_logFC, sum(res$significant, na.rm=TRUE)))

output_file <- paste0(args$output, ".intron_results.tsv")
cat(sprintf("Writing results to: %s\n", output_file))
write.table(res, output_file, sep="\t", quote=FALSE, row.names=FALSE)

sig_file <- paste0(args$output, ".significant_introns.tsv")
sig_results <- res[res$significant, ]
if (nrow(sig_results) > 0) {
  write.table(sig_results, sig_file, sep="\t", quote=FALSE, row.names=FALSE)
  cat(sprintf("Significant introns written to: %s\n", sig_file))
}
