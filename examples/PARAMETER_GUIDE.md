# Differential Splicing Analysis - Parameter Guide

## Input Modes

The pipeline can run in either matrix mode or BAM-manifest mode.

### Matrix Mode

Required inputs:

- `--matrix`: intron count matrix with introns as rows and samples as columns
- `--offset_matrix`: raw denominator matrix for the selected `--offset_mode`
- `--samples`: metadata TSV with `sample_id` and `group`
- `--contrast`: one explicit contrast, formatted `Treatment,Reference`

`--site_depth_offsets` is not supported; use `--offset_matrix` with the
`max_splice_plus_retained_depth.matrix` file.

Example:

```bash
python3 DSF.py \
    --matrix intron_counts.matrix \
    --offset_matrix intron_counts.max_splice_plus_retained_depth.matrix \
    --samples sample_metadata.tsv \
    --output_dir results \
    --contrast perturb,control
```

### BAM-Manifest Mode

If `--matrix` is omitted, `--samples` is interpreted as a BAM manifest:

```tsv
sample_type	replicate_id	bam_file
perturb	perturb_1	/path/to/perturb_1.bam
control	control_1	/path/to/control_1.bam
```

In this mode `--genome_fa` is required. The pipeline counts introns, creates
count and depth denominator matrices under `workdir/bam_inputs/`, writes a
downstream `sample_id/group` metadata file, and runs the analysis.

```bash
python3 DSF.py \
    --samples samples.tsv \
    --genome_fa reference.fa \
    --output_dir results \
    --contrast perturb,control \
    --site_depth_strand_mode RF
```

`--site_depth_strand_mode` can be `unstranded`, `F`, `R`, `FR`, or `RF`.
The default is `unstranded`.

This option controls strand filtering for **depth denominators**. Junction
counts are discovered from split alignments and then annotated by splice motif;
canonical junctions are effectively strand-resolved by their reference
dinucleotides, but the counting loop does not explicitly filter split reads by
library orientation. Use a stranded mode when local coverage from antisense or
overlapping genes could contaminate the denominator; use `unstranded` for
unstranded libraries.

When the intron-counting WDL is run with a stranded mode, it also emits
`*.transcript_plus.bam(.bai)` and `*.transcript_minus.bam(.bai)` outputs for the
alignments assigned to each transcript strand.

## Model

The pipeline tests each canonical intron once using a selected denominator:

```text
log(expected intron count) = group effect + log(selected_denominator)
```

Supported modes:

- `splice_plus_retained` (default): PSI/filtering and the model use
  `max_splice_plus_retained_depth`, the max of boundary splice-depth plus
  intron-side retained depth.
- `splice_vs_rest` (requires `--gtf`): PSI/filtering still use
  `max_splice_plus_retained_depth`, while the model uses a gene-scoped extension:
  `D = gene_total_junction_count + max(0, D^spr - Y_i)`. Introns with no gene
  assignment fall back to `splice_plus_retained`.

### Offset Modes

The default mode:

```bash
--offset_mode splice_plus_retained
```

Gene-scoped mode (requires `--gtf`):

```bash
--offset_mode splice_vs_rest \
--gtf annotation.gtf
```

Both modes use `max_splice_plus_retained_depth.matrix` as `--offset_matrix`.
The `splice_vs_rest` denominator is computed in-pipeline from the count matrix
and GTF.

### Statistical Modes

```bash
--stat_mode offset        # default: edgeR NB GLM with log-depth fixed offset
--stat_mode interaction   # focal/other interaction model (LRT)
--intx_engine edgeR       # default interaction engine
--intx_engine DEXSeq      # optional DEXSeq interaction engine
```

In `offset` mode with the SPR model denominator, `logFC` approximates
`log2(reported_PSI_A / reported_PSI_B)`. SVR `logFC` uses its gene-wide model
denominator instead.
In `interaction` mode, `logFC` is a log2 focal/other odds-ratio-like change.
`delta_PSI` is computed from raw counts and is unaffected by `--stat_mode`.

BAM-manifest mode computes `max_splice_plus_retained_depth` during the targeted
intron pass. Matrix mode requires `--offset_matrix` pointing to the appropriate
pre-computed denominator.

## Filtering

Filtering happens before statistical testing, so FDR is computed only over the
tested intron set. There is no post-test delta-PSI filtering or FDR
recomputation.

### Count Filters

- `--min_intron_count 10`: minimum total intron count across selected samples
- `--min_intron_samples 2`: minimum samples with nonzero intron count
Only canonical introns are tested in this offset-mode refactor.

### Denominator Depth Filters

- `--min_offset_depth 20`: minimum splice-plus-retained PSI denominator in a sample
- `--min_offset_samples 3`: minimum samples meeting `--min_offset_depth`

### Delta PSI Filter

- `--min_delta_psi 0.05`: minimum absolute group mean PSI difference required before testing
- `--min_delta_psi 0`: disable this prefilter

PSI is computed as:

```text
PSI = numerator_count / max_splice_plus_retained_depth
delta_PSI = mean_PSI_treatment - mean_PSI_reference
```

Only one contrast is supported per pipeline invocation.

## Statistical Mode Parameters

- `--stat_mode offset` (default): edgeR NB QL GLM; SPR-model `logFC` ≈ `log2(reported_PSI_A/reported_PSI_B)`
- `--stat_mode interaction --intx_engine edgeR`: edgeR stacked focal/other LRT; `logFC` is a log2 odds ratio
- `--stat_mode interaction --intx_engine DEXSeq`: DEXSeq focal/other LRT; `logFC` is a log2 odds-ratio-like effect

## edgeR Parameters

- `--contrast "GroupA,GroupB"`: required; logFC is GroupA / GroupB
- `--group_col group`: metadata column containing sample groups
- `--batch_col batch`: optional metadata column for batch correction
- `--fdr_threshold 0.05`: FDR cutoff for significance
- `--min_logFC 0.0`: minimum absolute edgeR log2 fold-change for significance

Run one explicit contrast per invocation.

## Depth Denominator Generation

Per-sample `.introns` files generated by `count_introns_from_bam.py` include a
set of depth columns:

- `left_adjacent_depth`, `right_adjacent_depth`, `max_adjacent_depth`
- `left_retained_depth`, `right_retained_depth`
- `left_splice_depth`, `right_splice_depth`
- `left_splice_plus_retained_depth`, `right_splice_plus_retained_depth`,
  `max_splice_plus_retained_depth`

Adjacent and retained depths are strand-aware read-depth windows from
`samtools depth`; adjacent windows are outside the intron and retained windows
are inside the intron. Retained windows start after
`--retained_depth_inner_offset` bases inside the intron, default `20`, and use
`--retained_depth_window_radius` bases, defaulting to
`--site_depth_window_radius`. Splice depths are sums of canonical intron counts
sharing the left or right boundary.

`build_intron_count_matrix.py` writes the count matrix plus one matrix for each
numeric depth column. `site_depth_offset` is retained as a compatibility alias
for `max_adjacent_depth`.

Use `--site_depth_strand_mode F`, `R`, `FR`, or `RF` with
`count_introns_from_bam.py` for strand-specific depth denominators.

For complete offsets across all samples, use a two-pass count:

```bash
# Discovery pass
python3 util/count_introns_from_bam.py \
    --genome_fa reference.fa \
    --bam sample1.bam > sample1.discovery.introns

python3 util/build_intron_count_matrix.py \
    --intron_files *.discovery.introns \
    --output_matrix discovery.matrix

# Targeted pass
python3 util/count_introns_from_bam.py \
    --genome_fa reference.fa \
    --bam sample1.bam \
    --target_introns discovery.matrix \
    --site_depth_strand_mode RF > sample1.targeted.introns

python3 util/build_intron_count_matrix.py \
    --intron_files *.targeted.introns \
    --output_matrix intron_counts.matrix
```

This writes `intron_counts.max_adjacent_depth.matrix`,
`intron_counts.max_splice_plus_retained_depth.matrix`, and the other
column-specific denominator matrices next to `intron_counts.matrix`.

## Outputs

Primary outputs:

- `edgeR_results.all.tsv`: all introns tested by edgeR, with PSI summaries
- `edgeR_results.significant_introns.tsv`: tested introns passing FDR and logFC thresholds

Key workdir outputs:

- `workdir/introns_filtered.tsv`: introns passing count, denominator-depth, and delta-PSI prefilters
- `workdir/site_depth_offsets.filtered.tsv`: raw selected denominators for tested introns
- `workdir/edgeR_input.counts.tsv`: edgeR count matrix
- `workdir/edgeR_input.offsets.tsv`: log selected-denominator matrix
- `workdir/edgeR_input.annotations.tsv`: intron annotations
- `workdir/psi.psi_values.tsv`: per-sample PSI, group means, and delta PSI
- `workdir/edgeR_results.intron_results.tsv`: raw model result table

Useful result columns:

- `logFC`: log2 effect size. Offset mode: ≈ PSI log-ratio. Interaction mode: log2 odds ratio.
- `FDR`: Benjamini-Hochberg adjusted p-value over the prefiltered tested set
- `delta_PSI`: group mean PSI difference (PSI-scale effect size for both stat modes)
- `*_mean_PSI`: group mean PSI values
- `offset_mode`: selected offset mode (`splice_plus_retained` or `splice_vs_rest`)
- `stat_mode`: statistical mode used
- `intx_engine`: interaction engine used (present in interaction-mode output)

## Suggested Settings

Balanced:

```bash
--min_intron_count 10 \
--min_intron_samples 2 \
--min_offset_depth 20 \
--min_offset_samples 3 \
--min_delta_psi 0.05
```

More lenient:

```bash
--min_intron_count 5 \
--min_intron_samples 2 \
--min_offset_depth 10 \
--min_offset_samples 2 \
--min_delta_psi 0.02
```

More conservative:

```bash
--min_intron_count 20 \
--min_offset_depth 50 \
--min_delta_psi 0.1 \
--min_logFC 0.5 \
--fdr_threshold 0.01
```
