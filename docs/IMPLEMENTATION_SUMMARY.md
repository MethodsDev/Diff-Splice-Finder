# Implementation Summary

## Overview

Diff-Splice-Finder is currently a single intron-level edgeR pipeline. Each
intron is tested once with a negative-binomial GLM and a user-selected exposure
offset. The main production mode uses read-depth denominators rather than
donor/acceptor cluster-total offsets.

The older `max(donor_cluster_total, acceptor_cluster_total)` implementation is
kept in utility code and historical docs, but it is not the main analysis path
used by `run_diff_splice_analysis.py`.

## Current Modes

### DSF exon-adjacent depth

```bash
--offset_mode exon_adjacent_depth
```

- Numerator: read-level split-junction counts from `count_introns_from_bam.py`
- PSI denominator: `max_adjacent_depth`
- edgeR offset: `log(max_adjacent_depth + 0.5)`
- Denominator definition:

```text
max_adjacent_depth = max(left_adjacent_depth, right_adjacent_depth)
```

Adjacent depths are exon-side windows outside the intron and are computed with
`samtools depth`. `site_depth_offset` remains a compatibility alias for
`max_adjacent_depth`.

### DSF splice-plus-retained

```bash
--offset_mode splice_plus_retained
```

- Numerator: read-level split-junction counts
- PSI denominator: `max_splice_plus_retained_depth`
- edgeR offset: `log(max_splice_plus_retained_depth + 0.5)`
- Denominator definition:

```text
max_splice_plus_retained_depth =
    max(left_splice_depth + left_retained_depth,
        right_splice_depth + right_retained_depth)
```

Splice depths are sums of canonical intron counts sharing the boundary. Retained
depths are intron-interior windows computed with `samtools depth`. By default,
retained windows start 20 bases inside the intron, controlled by
`--retained_depth_inner_offset`, and use the same width as
`--site_depth_window_radius` unless `--retained_depth_window_radius` is set.

### DSF splice-plus-retained beta-binomial-style

```bash
--offset_mode splice_plus_retained_betabinom
```

- Numerator/successes: read-level split-junction counts
- Trials: `max_splice_plus_retained_depth`
- Failures: `max_splice_plus_retained_depth - count`
- PSI denominator: `max_splice_plus_retained_depth`
- Statistical engine: base-R quasibinomial GLM over focal/rest counts

This mode is the depth-proxy analogue of the DSF-beta fate-pool idea. It does
not do fragment-name collapsing. The reported `logFC` is the model log2 odds
ratio, not an edgeR NB-offset coefficient.

### DSF gene-median splice-plus-retained

```bash
--offset_mode gene_median_splice_plus_retained
--gtf annotation.gtf
```

- Numerator: read-level split-junction counts
- PSI denominator: `max_splice_plus_retained_depth`
- edgeR offset: `log(per-gene median max_splice_plus_retained_depth + 0.5)`
- Fallback: introns without a gene assignment use their own
  `max_splice_plus_retained_depth`

In this mode, `delta_PSI` remains a splice-plus-retained PSI difference, while
`logFC` is the edgeR coefficient under the gene-median exposure. It should not be
read as the exact log-ratio of the reported PSI values.

## Input Modes

### Matrix mode

Matrix mode requires:

- `--matrix`: intron x sample count matrix
- `--offset_matrix`: intron x sample raw denominator matrix
- `--samples`: metadata with `sample_id` and the configured group column

`--site_depth_offsets` is retained as a deprecated compatibility alias for
`--offset_matrix` in `exon_adjacent_depth` mode.

### BAM-manifest mode

When `--matrix` is omitted, `--samples` is interpreted as a BAM manifest with:

```text
sample_type    replicate_id    bam_file
```

The pipeline runs a discovery pass, a targeted pass for complete depth
offsets, builds matrices under `workdir/bam_inputs/`, writes downstream
`sample_id/group` metadata, and then runs the selected statistical model.

## Filtering and FDR

Filtering happens before statistical testing:

- canonical splice filter; noncanonical introns are not supported by this
  offset-mode refactor
- total count and nonzero-sample filters
- selected denominator depth filter via `--min_offset_depth` and `--min_offset_samples`
- optional pre-test `--min_delta_psi` filter

edgeR's FDR is computed over the introns that pass these prefilters. The current
pipeline does not recompute FDR after edgeR. Final significant results are rows
marked significant by edgeR using the configured `--fdr_threshold` and
`--min_logFC`; if `--min_delta_psi` was nonzero, all tested rows already passed
that threshold.

## Strandedness

`--site_depth_strand_mode` can be `unstranded`, `F`, `R`, `FR`, or `RF`.

This option controls strand filtering for depth
denominators. Junction counts are discovered from split alignments without an
explicit read-orientation filter, then annotated by splice motif. Canonical
splice junctions are usually strand-resolved by their reference dinucleotides,
but overlapping antisense transcription can still matter for local coverage,
which is why stranded depth modes are useful.

## Core Files

- `run_diff_splice_analysis.py`: main orchestrator, filtering, mode selection,
  final output promotion
- `util/count_introns_from_bam.py`: read-level junction counts and depth
  denominator reporting
- `util/depth_windows.py`: samtools-depth adjacent/retained window depths
- `util/site_depth.py`: stranded read partitioning helpers used for stranded
  depth modes
- `util/run_edgeR_analysis.R`: edgeR GLM with supplied log offsets
- `util/run_splice_plus_retained_betabinom.R`: focal/rest quasibinomial GLM
- `util/build_intron_count_matrix.py`: count and depth matrix construction

## Output Structure

Primary outputs:

```text
results/
├── edgeR_results.all.tsv
└── edgeR_results.significant_introns.tsv
```

Key intermediates:

```text
results/workdir/
├── introns_filtered.tsv
├── site_depth_offsets.filtered.tsv
├── edgeR_input.counts.tsv
├── edgeR_input.offsets.tsv
├── edgeR_input.annotations.tsv
├── edgeR_results.intron_results.tsv
├── edgeR_results.intron_results_with_psi.tsv
├── psi.psi_values.tsv
└── bam_inputs/
    ├── intron_counts.matrix
    ├── intron_counts.offsets.matrix
    ├── intron_counts.max_adjacent_depth.matrix
    └── intron_counts.max_splice_plus_retained_depth.matrix
```
