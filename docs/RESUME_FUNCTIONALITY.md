# Resume Functionality

## Overview

The pipeline supports automatic resume. If a run is interrupted, rerun the same command and completed steps will be skipped based on existing output files.

## Current Checkpoints

The current single-workflow pipeline checkpoints these files under `workdir/`:

1. `introns_clustered.tsv`
2. `shared_offsets.tsv`
3. `introns_filtered.tsv`
4. `edgeR_input.counts.tsv`
5. `edgeR_input.offsets.tsv`
6. `edgeR_input.annotations.tsv`
7. `edgeR_results.intron_results.tsv`
8. `psi.psi_values.tsv`
9. `edgeR_results.intron_results_with_psi.tsv`
10. `edgeR_results.intron_results_with_psi.psi_filtered.tsv`

Final user-facing outputs are written in the main output directory:

- `edgeR_results.all.tsv`
- `edgeR_results.significant_introns.tsv`

A step is skipped if its checkpoint file already exists and is non-empty.

## Usage

```bash
python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/analysis
```

If the run is interrupted, rerun the same command:

```bash
python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/analysis
```

To force a full rerun:

```bash
python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/analysis \
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
rm results/analysis/workdir/edgeR_results.intron_results_with_psi.psi_filtered.tsv

python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/analysis
```

## Best Practice

- Reuse the exact same command when resuming.
- Check the logs to confirm which steps were skipped.
- Use a new output directory or `--force_rerun` after substantive changes.
