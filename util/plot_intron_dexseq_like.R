#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(rtracklayer)
  library(patchwork)
  library(ggrepel)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    y
  } else {
    x
  }
}

option_list <- list(
  make_option(
    c("--gene"),
    type = "character",
    help = "Target gene symbol to plot"
  ),
  make_option(
    c("--gtf"),
    type = "character",
    help = "Reference GTF used to draw transcript structure"
  ),
  make_option(
    c("--output"),
    type = "character",
    help = "Output PDF path"
  ),
  make_option(
    c("--results"),
    type = "character",
    default = NULL,
    help = "PSI-enhanced results TSV/TSV.GZ. If omitted, resolved from --output_dir"
  ),
  make_option(
    c("--output_dir"),
    type = "character",
    default = NULL,
    help = "Diff-Splice-Finder analysis output directory"
  ),
  make_option(
    c("--annotations"),
    type = "character",
    default = NULL,
    help = "Optional annotation TSV to merge intron_status/gene_name columns"
  ),
  make_option(
    c("--contrast"),
    type = "character",
    default = NULL,
    help = "Contrast to plot. Accepts either 'GroupA_vs_GroupB' or 'GroupA,GroupB1;GroupB2'"
  ),
  make_option(
    c("--fdr_threshold"),
    type = "double",
    default = 0.05,
    help = "FDR threshold used for highlighting"
  ),
  make_option(
    c("--min_delta_psi"),
    type = "double",
    default = 0.05,
    help = "Minimum absolute delta PSI used for highlighting"
  ),
  make_option(
    c("--intron_display_width"),
    type = "integer",
    default = NA,
    help = "Optional fixed display width for compressed introns"
  ),
  make_option(
    c("--use_fitted_logcpm"),
    action = "store_true",
    default = FALSE,
    help = "Use fitted logCPM summary columns instead of observed logCPM"
  )
)

parse_args_safe <- function() {
  parser <- OptionParser(option_list = option_list)
  args <- parse_args(parser)

  required <- c("gene", "gtf", "output")
  missing_args <- required[vapply(required, function(x) is.null(args[[x]]) || identical(args[[x]], ""), logical(1))]
  if (length(missing_args) > 0) {
    stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
  }

  args
}

normalize_contrast_label <- function(x) {
  if (is.null(x) || is.na(x) || x == "") {
    return(NULL)
  }
  if (str_detect(x, "_vs_")) {
    return(x)
  }

  parts <- str_split(x, ",", n = 2, simplify = TRUE)
  if (ncol(parts) != 2 || any(parts == "")) {
    stop("Contrast must be in format 'GroupA_vs_GroupB' or 'GroupA,GroupB1;GroupB2'")
  }

  group1 <- str_trim(parts[1])
  group2 <- str_replace_all(str_trim(parts[2]), ";", "_")
  paste0(group1, "_vs_", group2)
}

resolve_results_file <- function(output_dir, explicit_results) {
  if (!is.null(explicit_results)) {
    return(explicit_results)
  }
  if (is.null(output_dir)) {
    stop("Provide either --results or --output_dir")
  }

  candidates <- c(
    file.path(output_dir, "workdir", "edgeR_results.intron_results_with_psi.tsv"),
    file.path(output_dir, "edgeR_results.intron_results_with_psi.tsv"),
    file.path(output_dir, "workdir", "edgeR_results.intron_results_with_psi.tsv.gz"),
    file.path(output_dir, "edgeR_results.intron_results_with_psi.tsv.gz")
  )

  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop("Could not resolve PSI-enhanced results file from --output_dir")
  }

  existing[1]
}

resolve_annotations_file <- function(output_dir, explicit_annotations) {
  if (!is.null(explicit_annotations)) {
    return(explicit_annotations)
  }
  if (is.null(output_dir)) {
    return(NULL)
  }

  candidates <- c(
    file.path(output_dir, "workdir", "edgeR_input.annotations.tsv"),
    file.path(output_dir, "edgeR_input.annotations.tsv")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    return(NULL)
  }

  existing[1]
}

parse_intron_coords <- function(intron_ids) {
  motif_match <- str_match(intron_ids, "^([^:]+):(\\d+)-(\\d+)\\^([^\\^]+)\\^[^\\^]+$")
  standard_match <- str_match(intron_ids, "^([^:]+):(\\d+)-(\\d+):([+-])$")

  chrom <- ifelse(!is.na(motif_match[, 1]), motif_match[, 2], standard_match[, 2])
  start <- ifelse(!is.na(motif_match[, 1]), motif_match[, 3], standard_match[, 3])
  end <- ifelse(!is.na(motif_match[, 1]), motif_match[, 4], standard_match[, 4])
  strand <- ifelse(!is.na(motif_match[, 1]), NA_character_, standard_match[, 5])
  splice_sites <- ifelse(!is.na(motif_match[, 1]), motif_match[, 5], NA_character_)

  tibble(
    intron_id = intron_ids,
    chrom = chrom,
    intron_start = suppressWarnings(as.integer(start)),
    intron_end = suppressWarnings(as.integer(end)),
    strand = strand,
    splice_sites = splice_sites
  )
}

make_intron_compressed_transform <- function(exon_starts, exon_ends, intron_display_width = NULL) {
  ex <- tibble(start = exon_starts, end = exon_ends) |>
    arrange(start)

  cur_start <- ex$start[1]
  cur_end <- ex$end[1]
  merged <- list()

  if (nrow(ex) > 1) {
    for (i in seq_len(nrow(ex))[-1]) {
      if (ex$start[i] <= cur_end) {
        cur_end <- max(cur_end, ex$end[i])
      } else {
        merged[[length(merged) + 1]] <- c(cur_start, cur_end)
        cur_start <- ex$start[i]
        cur_end <- ex$end[i]
      }
    }
  }
  merged[[length(merged) + 1]] <- c(cur_start, cur_end)

  merged_df <- tibble(
    start = vapply(merged, `[[`, numeric(1), 1),
    end = vapply(merged, `[[`, numeric(1), 2)
  )

  exon_widths <- merged_df$end - merged_df$start
  if (is.null(intron_display_width) || is.na(intron_display_width)) {
    intron_display_width <- max(50, round(stats::median(exon_widths) * 0.2))
  }

  n <- nrow(merged_df)
  disp_start <- cumsum(c(0, head(exon_widths + intron_display_width, -1)))
  disp_end <- disp_start + exon_widths

  transform_pos <- function(p) {
    vapply(p, function(x) {
      in_ex <- which(x >= merged_df$start & x <= merged_df$end)
      if (length(in_ex) > 0) {
        k <- in_ex[1]
        disp_start[k] + (x - merged_df$start[k])
      } else {
        left <- which(merged_df$end < x)
        if (length(left) == 0L) {
          disp_start[1] - (merged_df$start[1] - x)
        } else {
          k <- max(left)
          if (k >= n) {
            disp_end[n] + (x - merged_df$end[n])
          } else {
            frac <- (x - merged_df$end[k]) / max(1L, merged_df$start[k + 1L] - merged_df$end[k])
            disp_end[k] + frac * intron_display_width
          }
        }
      }
    }, numeric(1))
  }

  intron_bands <- if (n > 1L) {
    tibble(xmin = disp_end[-n], xmax = disp_start[-1])
  } else {
    tibble(xmin = numeric(0), xmax = numeric(0))
  }

  list(transform = transform_pos, intron_bands = intron_bands, merged_exons = merged_df)
}

extract_gene_gtf <- function(gtf_file, gene_symbol, output_pdf) {
  cache_dir <- dirname(output_pdf)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(cache_dir, paste0(gene_symbol, ".subset.gtf"))

  if (!file.exists(cache_file) || file.info(cache_file)$size == 0) {
    pattern <- paste0("gene_name \"", gene_symbol, "\"")
    cmd <- sprintf("grep -F %s %s > %s",
                   shQuote(pattern), shQuote(gtf_file), shQuote(cache_file))
    status <- system(cmd)
    if (!identical(status, 0L) || !file.exists(cache_file) || file.info(cache_file)$size == 0) {
      stop("Failed to extract GTF rows for gene ", gene_symbol, " from ", gtf_file)
    }
  }

  gtf_tbl <- as_tibble(as.data.frame(rtracklayer::import(cache_file)))
  if (!"gene_name" %in% colnames(gtf_tbl)) {
    gtf_tbl$gene_name <- NA_character_
  }
  if (!"transcript_id" %in% colnames(gtf_tbl)) {
    gtf_tbl$transcript_id <- NA_character_
  }

  gtf_tbl |>
    mutate(
      seqnames = as.character(seqnames),
      strand = as.character(strand),
      type = as.character(type),
      gene_name = as.character(gene_name),
      transcript_id = as.character(transcript_id)
    )
}

merge_annotations <- function(results_df, annotations_file) {
  if (is.null(annotations_file) || !file.exists(annotations_file)) {
    return(results_df)
  }

  ann <- read_tsv(annotations_file, show_col_types = FALSE, name_repair = "minimal")
  if (!"intron_id" %in% colnames(ann) && ncol(ann) > 0) {
    colnames(ann)[1] <- "intron_id"
  }

  ann_keep <- ann |>
    select(any_of(c("intron_id", "gene_name", "intron_status", "overlapping_genes")))

  if (!"intron_id" %in% colnames(ann_keep)) {
    return(results_df)
  }

  merged <- results_df |>
    left_join(ann_keep, by = "intron_id", suffix = c("", ".ann"))

  if ("gene_name.ann" %in% colnames(merged)) {
    merged <- merged |>
      mutate(gene_name = coalesce(gene_name, gene_name.ann)) |>
      select(-gene_name.ann)
  }
  if ("intron_status.ann" %in% colnames(merged)) {
    merged <- merged |>
      mutate(intron_status = coalesce(intron_status, intron_status.ann)) |>
      select(-intron_status.ann)
  }
  if ("overlapping_genes.ann" %in% colnames(merged)) {
    merged <- merged |>
      mutate(overlapping_genes = coalesce(overlapping_genes, overlapping_genes.ann)) |>
      select(-overlapping_genes.ann)
  }

  merged
}

filter_introns_for_gene <- function(results_df, gene_symbol, gtf_df) {
  gene_match <- if ("gene_name" %in% colnames(results_df)) {
    results_df |>
      filter(!is.na(gene_name), gene_name == gene_symbol)
  } else {
    results_df[0, , drop = FALSE]
  }

  if (nrow(gene_match) > 0) {
    return(gene_match)
  }

  gene_span <- gtf_df |>
    filter(type %in% c("gene", "transcript", "exon")) |>
    summarise(
      seqnames = dplyr::first(seqnames),
      gene_start = min(start, na.rm = TRUE),
      gene_end = max(end, na.rm = TRUE)
    )

  results_df |>
    filter(
      chrom == gene_span$seqnames[1],
      intron_end >= gene_span$gene_start[1],
      intron_start <= gene_span$gene_end[1]
    )
}

pick_contrast <- function(results_df, contrast_arg) {
  if (!"contrast" %in% colnames(results_df)) {
    if (!is.null(contrast_arg)) {
      warning("Ignoring --contrast because results file has no contrast column")
    }
    return(results_df)
  }

  available <- unique(results_df$contrast)
  if (length(available) == 1 && is.null(contrast_arg)) {
    return(results_df)
  }

  contrast_label <- normalize_contrast_label(contrast_arg)
  if (is.null(contrast_label)) {
    stop("Results contain multiple contrasts. Provide --contrast. Available: ",
         paste(available, collapse = ", "))
  }
  if (!contrast_label %in% available) {
    stop("Requested contrast not found. Available: ", paste(available, collapse = ", "))
  }

  results_df |>
    filter(contrast == contrast_label)
}

build_label_from_splice_sites <- function(splice_sites, intron_start, intron_end) {
  splice_label <- ifelse(is.na(splice_sites), "NA", str_replace(splice_sites, "--", "-"))
  paste0(intron_start, "-", intron_end, " [", splice_label, "]")
}

derive_group_labels <- function(contrast_label) {
  if (is.null(contrast_label) || !str_detect(contrast_label, "_vs_")) {
    return(list(group1 = "group1", group2 = "group2"))
  }

  parts <- str_split(contrast_label, "_vs_", n = 2, simplify = TRUE)
  list(group1 = parts[1], group2 = parts[2])
}

args <- parse_args_safe()

results_file <- resolve_results_file(args$output_dir, args$results)
annotations_file <- resolve_annotations_file(args$output_dir, args$annotations)

message("Loading results: ", results_file)
results_df <- read_tsv(results_file, show_col_types = FALSE, name_repair = "minimal")

required_cols <- c("intron_id", "logFC", "FDR")
missing_cols <- required_cols[!required_cols %in% colnames(results_df)]
if (length(missing_cols) > 0) {
  stop("Results file is missing required columns: ", paste(missing_cols, collapse = ", "))
}
if (!"delta_PSI" %in% colnames(results_df)) {
  stop("Results file is missing delta_PSI. Use the PSI-enhanced results TSV.")
}

results_df <- merge_annotations(results_df, annotations_file)
results_df <- pick_contrast(results_df, args$contrast)

coord_df <- parse_intron_coords(results_df$intron_id)
if (!"gene_name" %in% colnames(results_df)) {
  results_df$gene_name <- NA_character_
}
if (!"intron_status" %in% colnames(results_df)) {
  results_df$intron_status <- NA_character_
}
if (!"overlapping_genes" %in% colnames(results_df)) {
  results_df$overlapping_genes <- NA_character_
}
results_df <- results_df |>
  left_join(coord_df, by = "intron_id") |>
  mutate(
    neg_log_fdr = -log10(pmax(FDR, 1e-300)),
    intron_status = coalesce(intron_status, "unknown"),
    intron_label = build_label_from_splice_sites(splice_sites, intron_start, intron_end)
  )

if (any(is.na(results_df$intron_start) | is.na(results_df$intron_end))) {
  stop("Failed to parse intron coordinates from at least one intron_id")
}

message("Extracting transcript structure for gene: ", args$gene)
gtf_df <- extract_gene_gtf(args$gtf, args$gene, args$output)

if (nrow(gtf_df) == 0) {
  stop("No GTF features found for gene ", args$gene)
}

gene_introns <- filter_introns_for_gene(results_df, args$gene, gtf_df)
if (nrow(gene_introns) == 0) {
  stop("No introns found for gene ", args$gene, " in the selected results")
}

exons_all <- gtf_df |>
  filter(type == "exon")
if (nrow(exons_all) == 0) {
  stop("No exon features found in GTF subset for gene ", args$gene)
}

coord <- make_intron_compressed_transform(
  exon_starts = exons_all$start,
  exon_ends = exons_all$end,
  intron_display_width = if (is.na(args$intron_display_width)) NULL else args$intron_display_width
)
tf <- coord$transform

gene_introns <- gene_introns |>
  mutate(
    start_disp = tf(intron_start),
    end_disp = tf(intron_end),
    mid_disp = (start_disp + end_disp) / 2,
    is_highlight = FDR < args$fdr_threshold & abs(delta_PSI) >= args$min_delta_psi
  )

txs <- gtf_df |>
  filter(type == "transcript") |>
  transmute(
    transcript_id = coalesce(transcript_id, args$gene),
    start = tf(start),
    end = tf(end),
    strand = as.character(strand)
  ) |>
  distinct() |>
  mutate(y = row_number())

if (nrow(txs) == 0) {
  txs <- exons_all |>
    summarise(start = min(start), end = max(end)) |>
    mutate(
      transcript_id = args$gene,
      start = tf(start),
      end = tf(end),
      strand = as.character(dplyr::first(gtf_df$strand)),
      y = 1
    ) |>
    select(transcript_id, start, end, strand, y)
}

exons <- exons_all |>
  select(transcript_id, start, end) |>
  mutate(transcript_id = coalesce(transcript_id, args$gene)) |>
  left_join(txs |> select(transcript_id, y), by = "transcript_id") |>
  mutate(start = tf(start), end = tf(end))

cds_data <- gtf_df |>
  filter(type == "CDS") |>
  select(transcript_id, start, end) |>
  mutate(transcript_id = coalesce(transcript_id, args$gene)) |>
  left_join(txs |> select(transcript_id, y), by = "transcript_id") |>
  mutate(start = tf(start), end = tf(end))

contrast_label <- if ("contrast" %in% colnames(gene_introns)) unique(gene_introns$contrast) else NULL
if (length(contrast_label) > 1) {
  stop("Filtered rows still contain multiple contrasts; this should not happen")
}
contrast_label <- if (length(contrast_label) == 0) NULL else contrast_label[1]
group_labels <- derive_group_labels(contrast_label)

group1_col <- if (args$use_fitted_logcpm) "contrast_group1_mean_fitted_logCPM" else "contrast_group1_mean_logCPM"
group2_col <- if (args$use_fitted_logcpm) "contrast_group2_mean_fitted_logCPM" else "contrast_group2_mean_logCPM"
missing_cpm_cols <- c(group1_col, group2_col)[!c(group1_col, group2_col) %in% colnames(gene_introns)]
if (length(missing_cpm_cols) > 0) {
  stop("Results file is missing required logCPM columns: ", paste(missing_cpm_cols, collapse = ", "))
}

shared_x <- range(
  c(gene_introns$start_disp, gene_introns$end_disp, txs$start, txs$end, exons$start, exons$end),
  na.rm = TRUE
)
x_pad <- diff(shared_x) * 0.02
shared_x_lim <- shared_x + c(-x_pad, x_pad)

fdr_range <- range(gene_introns$neg_log_fdr, na.rm = TRUE)
if (fdr_range[1] == fdr_range[2]) {
  fdr_range <- fdr_range + c(0, 1)
}

fitted_long <- gene_introns |>
  select(
    intron_id, intron_label, intron_status, significant, FDR, delta_PSI,
    mid_disp, start_disp, end_disp, is_highlight,
    group1_value = all_of(group1_col),
    group2_value = all_of(group2_col)
  ) |>
  pivot_longer(
    cols = c(group1_value, group2_value),
    names_to = "group",
    values_to = "logCPM_value"
  ) |>
  mutate(
    group_label = if_else(group == "group1_value", group_labels$group1, group_labels$group2),
    group_label = factor(group_label, levels = c(group_labels$group1, group_labels$group2))
  )

cpm_range <- range(fitted_long$logCPM_value, na.rm = TRUE)
cpm_pad <- if (diff(cpm_range) == 0) 1 else diff(cpm_range) * 0.15
y_lim_cpm <- cpm_range + c(-cpm_pad, cpm_pad)

strand_dir <- unique(txs$strand)
strand_dir <- strand_dir[!is.na(strand_dir)]
strand_dir <- if (length(strand_dir) == 0) "." else strand_dir[1]

top_panel_label <- if (args$use_fitted_logcpm) "Mean fitted logCPM" else "Mean logCPM"

p_cpm <- ggplot() +
  geom_rect(
    data = coord$intron_bands,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = "grey90", alpha = 0.7, inherit.aes = FALSE
  ) +
  geom_rect(
    data = gene_introns |> filter(is_highlight),
    aes(
      xmin = start_disp,
      xmax = end_disp,
      ymin = pmin(.data[[group1_col]], .data[[group2_col]]),
      ymax = pmax(.data[[group1_col]], .data[[group2_col]])
    ),
    fill = "yellow", alpha = 0.35, inherit.aes = FALSE
  ) +
  geom_segment(
    data = gene_introns,
    aes(
      x = mid_disp, xend = mid_disp,
      y = .data[[group1_col]], yend = .data[[group2_col]],
      alpha = is_highlight
    ),
    color = "grey70", linewidth = 0.4
  ) +
  geom_segment(
    data = fitted_long,
    aes(
      x = start_disp, xend = end_disp,
      y = logCPM_value, yend = logCPM_value,
      color = group_label, linetype = intron_status, linewidth = intron_status, alpha = is_highlight
    ),
    lineend = "butt"
  ) +
  geom_segment(
    data = fitted_long,
    aes(
      x = start_disp, xend = start_disp,
      y = logCPM_value - 0.12, yend = logCPM_value + 0.12,
      color = group_label, alpha = is_highlight
    ),
    linewidth = 0.8
  ) +
  geom_segment(
    data = fitted_long,
    aes(
      x = end_disp, xend = end_disp,
      y = logCPM_value - 0.12, yend = logCPM_value + 0.12,
      color = group_label, alpha = is_highlight
    ),
    linewidth = 0.8
  ) +
  geom_text_repel(
    data = fitted_long |>
      filter(FDR < args$fdr_threshold, group == "group1_value", abs(delta_PSI) >= args$min_delta_psi),
    aes(x = mid_disp, y = logCPM_value, label = intron_label, color = group_label),
    size = 2.1, fontface = "bold",
    box.padding = 0.4, point.padding = 0.2,
    min.segment.length = 0.2, segment.color = "grey55", segment.size = 0.3,
    direction = "y", force = 2, nudge_y = cpm_pad * 0.5,
    show.legend = FALSE, max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(setNames("steelblue", group_labels$group1), setNames("firebrick", group_labels$group2)),
    name = "Group"
  ) +
  scale_linetype_manual(
    values = c("known" = "solid", "novel" = "22", "unknown" = "dotted"),
    name = "Intron status"
  ) +
  scale_linewidth_manual(
    values = c("known" = 1.0, "novel" = 0.5, "unknown" = 0.6),
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c("TRUE" = 1.0, "FALSE" = 0.15),
    guide = "none"
  ) +
  scale_x_continuous(limits = shared_x_lim, expand = c(0, 0)) +
  scale_y_continuous(limits = y_lim_cpm) +
  labs(
    title = paste0(args$gene, "  |  ", strand_dir, " strand  |  ", contrast_label %||% "single contrast",
                   "  |  intron ", top_panel_label),
    subtitle = paste0(
      nrow(gene_introns), " introns tested; ",
      sum(gene_introns$FDR < args$fdr_threshold, na.rm = TRUE),
      " significant at FDR < ", args$fdr_threshold
    ),
    x = NULL,
    y = top_panel_label
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

logfc_range <- range(gene_introns$logFC, na.rm = TRUE)
logfc_pad <- if (diff(logfc_range) == 0) 1 else diff(logfc_range) * 0.55
y_lim_dtu <- logfc_range + c(-logfc_pad, logfc_pad)
dtu_guides <- c(-1, 1)
dtu_guides <- dtu_guides[dtu_guides >= y_lim_dtu[1] & dtu_guides <= y_lim_dtu[2]]

p_dtu <- ggplot() +
  geom_rect(
    data = coord$intron_bands,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = "grey90", alpha = 0.7, inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.5) +
  geom_segment(
    data = gene_introns,
    aes(
      x = start_disp, xend = end_disp,
      y = logFC, yend = logFC,
      color = neg_log_fdr, linetype = intron_status, linewidth = intron_status
    ),
    lineend = "butt"
  ) +
  geom_segment(
    data = gene_introns,
    aes(x = start_disp, xend = start_disp, y = logFC - 0.12, yend = logFC + 0.12, color = neg_log_fdr),
    linewidth = 0.8
  ) +
  geom_segment(
    data = gene_introns,
    aes(x = end_disp, xend = end_disp, y = logFC - 0.12, yend = logFC + 0.12, color = neg_log_fdr),
    linewidth = 0.8
  ) +
  geom_text_repel(
    data = gene_introns |> filter(is_highlight),
    aes(
      x = mid_disp, y = logFC,
      label = sprintf(
        "%s\nlogFC=%.1f  FDR=%s\nDeltaPSI=%.3f",
        intron_label, logFC, formatC(FDR, format = "e", digits = 1), delta_PSI
      ),
      color = neg_log_fdr
    ),
    size = 2.2, fontface = "bold", lineheight = 0.95,
    box.padding = 0.5, point.padding = 0.2,
    min.segment.length = 0.2, segment.color = "grey55", segment.size = 0.35,
    direction = "both", force = 2,
    show.legend = FALSE, max.overlaps = Inf
  ) +
  scale_color_gradient(low = "grey80", high = "red3", name = expression(-log[10](FDR)), limits = fdr_range) +
  scale_linetype_manual(values = c("known" = "solid", "novel" = "22", "unknown" = "dotted"), name = "Intron status") +
  scale_linewidth_manual(values = c("known" = 1.0, "novel" = 0.5, "unknown" = 0.6), guide = "none") +
  scale_x_continuous(limits = shared_x_lim, expand = c(0, 0)) +
  scale_y_continuous(limits = y_lim_dtu) +
  labs(x = NULL, y = "logFC (DTU)") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

if (length(dtu_guides) > 0) {
  p_dtu <- p_dtu +
    geom_hline(yintercept = dtu_guides, linetype = "dashed", color = "grey65", linewidth = 0.35)
}
if (1 >= y_lim_dtu[1] && 1 <= y_lim_dtu[2]) {
  p_dtu <- p_dtu +
    annotate("text", x = Inf, y = 1.12, label = "logFC = +1", hjust = 1.05, size = 2.3, color = "grey55")
}
if (-1 >= y_lim_dtu[1] && -1 <= y_lim_dtu[2]) {
  p_dtu <- p_dtu +
    annotate("text", x = Inf, y = -1.12, label = "logFC = -1", hjust = 1.05, size = 2.3, color = "grey55")
}

p_struct <- ggplot() +
  geom_rect(
    data = coord$intron_bands,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = "grey90", alpha = 0.7, inherit.aes = FALSE
  ) +
  geom_segment(
    data = txs,
    aes(x = start, xend = end, y = y, yend = y),
    color = "grey55", linewidth = 0.4
  ) +
  geom_rect(
    data = exons,
    aes(xmin = start, xmax = end, ymin = y - 0.20, ymax = y + 0.20),
    fill = "grey72", color = "grey40", linewidth = 0.2
  ) +
  geom_rect(
    data = cds_data,
    aes(xmin = start, xmax = end, ymin = y - 0.38, ymax = y + 0.38),
    fill = "grey55", color = "grey35", linewidth = 0.2
  ) +
  scale_x_continuous(limits = shared_x_lim, expand = c(0, 0)) +
  scale_y_continuous(
    breaks = txs$y,
    labels = txs$transcript_id,
    limits = c(min(txs$y) - 0.8, max(txs$y) + 0.8),
    expand = c(0, 0)
  ) +
  labs(x = "Exon structure (introns compressed)", y = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank()
  )

p_combined <- p_cpm / p_dtu / p_struct +
  plot_layout(heights = c(2, 2, 1), guides = "collect")

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
n_txs <- nrow(txs)
n_intr <- nrow(gene_introns)
plot_h <- min(50, max(8, n_txs * 0.35 + n_intr * 0.45 + 7))

ggsave(
  filename = args$output,
  plot = p_combined,
  width = 14,
  height = plot_h,
  limitsize = FALSE
)

message("Saved: ", args$output)
