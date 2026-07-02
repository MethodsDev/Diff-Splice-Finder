# Differential Splicing Analysis - Parameter Guide

## Input Modes

The pipeline can run in either matrix mode or BAM-manifest mode.

### Matrix Mode

Required inputs:

- `--matrix`: intron count matrix with introns as rows and samples as columns
- `--offset_matrix`: raw denominator matrix for the selected `--offset_mode`
- `--samples`: metadata TSV with `sample_id` and `group`
- `--contrast`: one explicit contrast, formatted `Treatment,Reference`

`--site_depth_offsets` remains available as a deprecated alias for
`--offset_matrix` in `exon_adjacent_depth` mode.

Example:

```bash
python3 run_diff_splice_analysis.py \
    --matrix intron_counts.matrix \
    --offset_matrix intron_counts.offsets.matrix \
    --samples sample_metadata.tsv \
    --output_dir results \
    --contrast perturb,control \
    --offset_mode exon_adjacent_depth
```

### BAM-Manifest Mode

If `--matrix` is omitted, `--samples` is interpreted as a BAM manifest:

```tsv
sample_type	replicate_id	bam_file
perturb	perturb_1	/path/to/perturb_1.bam
control	control_1	/path/to/control_1.bam
```

In this mode `--genome_fa` is required. The pipeline counts introns, creates
count and site-depth offset matrices under `workdir/bam_inputs/`, writes a
downstream `sample_id/group` metadata file, and runs the analysis.

```bash
python3 run_diff_splice_analysis.py \
    --samples samples.tsv \
    --genome_fa reference.fa \
    --output_dir results \
    --contrast perturb,control \
    --site_depth_strand_mode RF
```

`--site_depth_strand_mode` can be `unstranded`, `F`, `R`, `FR`, or `RF`.
The default is `unstranded`.

This option controls strand filtering for **site-depth offsets**. Junction
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

- `exon_adjacent_depth`: PSI and edgeR offset use `max_adjacent_depth`, the max
  of left/right exon-side adjacent read-depth windows.
- `splice_plus_retained`: PSI and edgeR offset use
  `max_splice_plus_retained_depth`, the max of boundary splice-depth plus
  intron-side retained depth.
- `gene_median_splice_plus_retained`: PSI uses
  `max_splice_plus_retained_depth`; edgeR uses the per-gene median of that
  denominator. Requires `--gtf`.

Adjacent and retained depths are computed with `samtools depth`. Splice depths
are derived from counts of canonical introns sharing each boundary.

### Offset Modes

The default mode is:

```bash
--offset_mode exon_adjacent_depth
```

Alternative modes:

```bash
--offset_mode splice_plus_retained

--offset_mode gene_median_splice_plus_retained \
--gtf annotation.gtf
```

BAM-manifest mode computes all denominator matrices during the targeted intron
pass. Matrix mode can use any mode if the corresponding denominator is supplied
with `--offset_matrix`.

## Filtering

Filtering happens before edgeR, so FDR is computed only over the tested intron
set. There is no post-edgeR delta-PSI filtering or FDR recomputation.

### Count Filters

- `--min_intron_count 10`: minimum total intron count across selected samples
- `--min_intron_samples 2`: minimum samples with nonzero intron count
Noncanonical introns are not supported in this offset-mode refactor.

### Denominator Depth Filters

- `--min_offset_depth 20`: minimum selected denominator in a sample
- `--min_offset_samples 3`: minimum samples meeting `--min_offset_depth`

### Delta PSI Filter

- `--min_delta_psi 0.05`: minimum absolute group mean PSI difference required before edgeR
- `--min_delta_psi 0`: disable this prefilter

PSI is computed as:

```text
PSI = numerator_count / selected_denominator
delta_PSI = mean_PSI_treatment - mean_PSI_reference
```

For `gene_median_splice_plus_retained`, `selected_denominator` for PSI is still
`max_splice_plus_retained_depth`; only the edgeR exposure is replaced by the
gene median.

Only one contrast is supported per pipeline invocation.

## edgeR Parameters

- `--contrast "GroupA,GroupB"`: required; logFC is GroupA / GroupB
- `--group_col group`: metadata column containing sample groups
- `--batch_col batch`: optional metadata column for batch correction
- `--fdr_threshold 0.05`: FDR cutoff for significance
- `--min_logFC 0.0`: minimum absolute edgeR log2 fold-change for significance

Run one explicit contrast per invocation.

## Site-Depth Offset Generation

Per-sample `.introns` files generated by `count_introns_from_bam.py` include a
`site_depth_offset` column. `build_intron_count_matrix.py` writes both count and
offset matrices when that column is present.

Use `--site_depth_strand_mode F`, `R`, `FR`, or `RF` with
`count_introns_from_bam.py` for strand-specific site-depth offsets.

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
    --output_matrix intron_counts.matrix \
    --output_offset_matrix intron_counts.offsets.matrix
```

## Outputs

Primary outputs:

- `edgeR_results.all.tsv`: all introns tested by edgeR, with PSI summaries
- `edgeR_results.significant_introns.tsv`: tested introns passing FDR and logFC thresholds

Key workdir outputs:

- `workdir/introns_filtered.tsv`: introns passing count, site-depth, and delta-PSI prefilters
- `workdir/site_depth_offsets.filtered.tsv`: raw offsets for tested introns
- `workdir/edgeR_input.counts.tsv`: edgeR count matrix
- `workdir/edgeR_input.offsets.tsv`: log site-depth offset matrix
- `workdir/edgeR_input.annotations.tsv`: intron annotations
- `workdir/psi.psi_values.tsv`: per-sample PSI, group means, and delta PSI
- `workdir/edgeR_results.intron_results.tsv`: raw edgeR result table

Useful result columns:

- `logFC`: edgeR-estimated log2 fold-change in offset-adjusted intron usage
- `FDR`: Benjamini-Hochberg adjusted p-value over the prefiltered tested set
- `delta_PSI`: group mean PSI difference used during prefiltering
- `*_mean_PSI`: group mean PSI values
- `offset_mode`, `offset_source`: `site_depth`

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
