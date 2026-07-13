# Resume Functionality

## Overview

The pipeline supports automatic resume. If a run is interrupted, rerun the same command and completed steps will be skipped based on existing output files.

## Current Checkpoints

The current single-workflow pipeline checkpoints these files under `workdir/`:

1. `introns_filtered.tsv`
2. `site_depth_offsets.filtered.tsv`
3. `psi.psi_values.tsv`
4. `edgeR_input.counts.tsv`
5. `edgeR_input.offsets.tsv`
6. `edgeR_input.annotations.tsv`
7. `edgeR_input.filter_params.json`
8. `edgeR_results.intron_results.tsv`
9. `edgeR_results.intron_results_with_psi.tsv`
10. `edgeR_results.params.json`

In BAM-manifest mode, input preparation also checkpoints files under
`workdir/bam_inputs/`, including `intron_counts.matrix`,
`intron_counts.offsets.matrix`, and column-specific depth matrices such as
`intron_counts.max_adjacent_depth.matrix` and
`intron_counts.max_splice_plus_retained_depth.matrix`.

Final user-facing outputs are written in the main output directory:

- `edgeR_results.all.tsv`
- `edgeR_results.significant_introns.tsv`

A step is skipped if its checkpoint file already exists and is non-empty.

## Usage

```bash
python3 DSF.py \
    --matrix data/intron_counts.matrix \
    --offset_matrix data/intron_counts.max_splice_plus_retained_depth.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/analysis \
    --contrast perturb,control \
    --offset_mode splice_plus_retained
```

If the run is interrupted, rerun the same command:

```bash
python3 DSF.py \
    --matrix data/intron_counts.matrix \
    --offset_matrix data/intron_counts.max_splice_plus_retained_depth.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/analysis \
    --contrast perturb,control \
    --offset_mode splice_plus_retained
```

To force a full rerun:

```bash
python3 DSF.py \
    --matrix data/intron_counts.matrix \
    --offset_matrix data/intron_counts.max_splice_plus_retained_depth.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/analysis \
    --contrast perturb,control \
    --offset_mode splice_plus_retained \
    --force_rerun
```

## Limitations

Resume checks only:

- file existence
- file size greater than zero

It does not verify:

- parameter consistency
- file integrity
- whether code changed since the previous run

Use `--force_rerun` after changing inputs, thresholds, contrasts, or code.

## Manual Recovery

To rerun downstream analysis from the edgeR stage onward, remove the relevant files from `workdir/` and rerun the pipeline.

Example:

```bash
rm results/analysis/workdir/edgeR_results.intron_results.tsv
rm results/analysis/workdir/psi.psi_values.tsv
rm results/analysis/workdir/edgeR_results.intron_results_with_psi.tsv

python3 DSF.py \
    --matrix data/intron_counts.matrix \
    --offset_matrix data/intron_counts.max_splice_plus_retained_depth.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/analysis \
    --contrast perturb,control \
    --offset_mode splice_plus_retained
```

## Best Practice

- Reuse the exact same command when resuming.
- Check the logs to confirm which steps were skipped.
- Use a new output directory or `--force_rerun` after substantive changes.
