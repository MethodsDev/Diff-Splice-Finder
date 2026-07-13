# Implementation Summary

## Overview

Diff-Splice-Finder is a single intron-level statistical pipeline. Each intron is
tested once using a user-selected depth denominator and a user-selected
statistical model. Two offset modes and two statistical modes are supported.

## Offset Modes (`--offset_mode`)

### splice_plus_retained (default)

```bash
--offset_mode splice_plus_retained
```

- Numerator: read-level split-junction counts from `count_introns_from_bam.py`
- PSI denominator: `max_splice_plus_retained_depth`
- Offset (offset stat mode): `log(max_splice_plus_retained_depth + 0.5)`
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

### splice_vs_rest (requires `--gtf`)

```bash
--offset_mode splice_vs_rest \
--gtf annotation.gtf
```

- Numerator: read-level split-junction counts
- PSI denominator: `D^svr` (see below)
- Offset (offset stat mode): `log(D^svr + 0.5)`
- Denominator definition (computed in-pipeline from count matrix + GTF):

```text
D^svr_i,s = gene_total_junction_count_i,s
          + max(0, max_splice_plus_retained_depth_i,s - Y_i,s)
```

where `gene_total_junction_count_i,s = sum of all canonical intron counts in
gene(i) for sample s`. The second term is the local spr remainder (competing
site splices + retained depth after removing the focal count). Introns with no
gene assignment fall back to `splice_plus_retained` automatically.

Input `--offset_matrix` should be the same `max_splice_plus_retained_depth.matrix`
used for `splice_plus_retained`; the `D^svr` computation uses the count matrix
and GTF internally.

## Statistical Modes (`--stat_mode`)

### offset (default)

```bash
--stat_mode offset
```

The log-transformed denominator is supplied as a fixed edgeR offset. Each
intron is tested with edgeR’s quasi-likelihood NB GLM:

```text
log(E[Y_i,s]) = group_s + log(D_i,s)
```

`logFC` ≈ `log2(PSI_A / PSI_B)`. Dispersion is estimated with `glmQLFit`;
testing uses `glmQLFTest`.

### interaction

```bash
--stat_mode interaction
```

DEXSeq-style doubled count matrix. For each intron `i`, two pseudo-sample
columns are created per original sample `s`:

```text
focal column:  Y_i,s
other column:  max(0, D_i,s - Y_i,s)
```

Design: `~ sample_id + exon_type + group:exon_type`

The `sample_id` blocking factor absorbs per-sample depth variation; no fixed
offset is needed. The `group:exon_type` interaction is tested by LRT
(`glmFit` + `glmLRT`). `logFC` is a log2 focal/other odds-ratio change—larger
in magnitude than offset-mode `logFC` for the same PSI shift when baseline PSI
is high. `delta_PSI` is computed from raw counts before statistical testing and
is unaffected by `--stat_mode`.

Batch correction (`--batch_col`) adds `batch:exon_type` to both full and reduced
models.

## Input Modes

### Matrix mode

Matrix mode requires:

- `--matrix`: intron x sample count matrix
- `--offset_matrix`: intron x sample raw `max_splice_plus_retained_depth` matrix
- `--samples`: metadata with `sample_id` and the configured group column
- `--gtf`: required for `splice_vs_rest` mode

### BAM-manifest mode

When `--matrix` is omitted, `--samples` is interpreted as a BAM manifest with:

```text
sample_type    replicate_id    bam_file
```

The pipeline runs a discovery pass, a targeted pass for complete depth offsets,
builds matrices under `workdir/bam_inputs/`, writes downstream `sample_id/group`
metadata, and then runs the selected statistical model.

## Filtering and FDR

Filtering happens before statistical testing:

- canonical splice filter
- total count and nonzero-sample filters (`--min_intron_count`, `--min_intron_samples`)
- selected denominator depth filter (`--min_offset_depth`, `--min_offset_samples`)
- optional pre-test `--min_delta_psi` filter

FDR is computed over the introns that pass these prefilters. The pipeline does
not recompute FDR after statistical testing.

## Core Files

- `DSF.py`: main orchestrator, mode selection, filtering,
  svr denominator computation, output promotion
- `util/count_introns_from_bam.py`: read-level junction counts and depth
  denominator reporting
- `util/depth_windows.py`: samtools-depth adjacent/retained window depths
- `util/site_depth.py`: stranded read partitioning helpers
- `util/run_edgeR_analysis.R`: edgeR QL GLM with supplied log offsets
- `util/run_edgeR_interaction_analysis.R`: DEXSeq-style stacked NB interaction model
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
    └── intron_counts.max_splice_plus_retained_depth.matrix
```
