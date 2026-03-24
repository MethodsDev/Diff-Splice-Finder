# Session Memory

This file is a reusable brain dump for future work on `Diff-Splice-Finder`.

## Current Pipeline Model

- The pipeline is now a single intron-level workflow.
- Introns are clustered by both donor and acceptor.
- Shared offsets are computed as `max(donor_total, acceptor_total)` per intron per sample.
- edgeR is run once per intron with those shared offsets.
- PSI uses the same shared denominators as edgeR.

Main entrypoint:
- `run_diff_splice_analysis.py`

Key utility modules:
- `util/cluster_introns.py`
- `util/filter_introns.py`
- `util/compute_offsets.py`
- `util/run_edgeR_analysis.R`
- `util/compute_psi.py`

## Important Output Conventions

Final outputs go in the user-specified output directory:
- `edgeR_results.intron_results.tsv`
- `edgeR_results.significant_introns.tsv`
- `edgeR_results.diagnostics.pdf`
- `edgeR_results.RData`

Intermediate outputs go in:
- `<output_dir>/workdir/`

Important intermediate files:
- `workdir/introns_clustered.tsv`
- `workdir/shared_offsets.tsv`
- `workdir/introns_filtered.tsv`
- `workdir/edgeR_input.counts.tsv`
- `workdir/edgeR_input.offsets.tsv`
- `workdir/edgeR_input.annotations.tsv`
- `workdir/psi.psi_values.tsv`
- `workdir/edgeR_results.intron_results_with_psi.tsv`
- `workdir/edgeR_results.intron_results_with_psi.psi_filtered.tsv` if `--min_delta_psi` is used

## Contrast Contract

Use comma-delimited contrasts:
- Single control: `GroupA,GroupB`
- Pooled controls: `GroupA,GroupB1;GroupB2`

This format is now aligned across:
- `run_diff_splice_analysis.py`
- `util/run_edgeR_analysis.R`
- `util/compute_psi.py`

Default all-pairwise comparisons also now generate comma-delimited contrasts.

## Recent Fixes

These were implemented and pushed in commit `6c9e462`:

1. Fixed contrast format drift.
- The all-pairwise path in `run_diff_splice_analysis.py` had been generating `A-B`.
- edgeR and PSI parsing expect `A,B`.
- That mismatch is fixed.

2. Fixed `delta_PSI` for explicit contrasts.
- `util/compute_psi.py` had been parsing hyphen-delimited contrasts only.
- It now handles:
  - `GroupA,GroupB`
  - `GroupA,GroupB1;GroupB2`
- For pooled controls, `delta_PSI` is computed as:
  - `GroupA_mean_PSI - mean(control_group_means)`

3. Added per-contrast expression summary columns to edgeR results.
- Current columns added by `util/run_edgeR_analysis.R`:
  - `contrast_group1_mean_count`
  - `contrast_group2_mean_count`
  - `contrast_group1_mean_logCPM`
  - `contrast_group2_mean_logCPM`
  - `contrast_group1_mean_fitted_logCPM`
  - `contrast_group2_mean_fitted_logCPM`

Observed logCPM is currently computed from the count matrix.
Fitted logCPM is currently computed from `fit$fitted.values`.

## Verified Behavior

The following path was rerun successfully after the recent fixes:

```bash
cd testing
../run_diff_splice_analysis.py \
  --matrix test_intron_counts.matrix \
  --samples test_metadata_control.tsv \
  --output_dir quick_test_contrast_review \
  --contrast 'TDP43,control' \
  --min_intron_count 5 \
  --min_intron_samples 2 \
  --min_cluster_count 10 \
  --min_cluster_samples 2 \
  --fdr_threshold 0.05 \
  --cpu 1 \
  --force_rerun
```

This verified:
- explicit comma-delimited contrasts work
- `delta_PSI` is present for explicit contrasts
- mean count columns are present
- mean logCPM columns are present
- mean fitted logCPM columns are present

## Testing Notes

- `testing/run_quick_test.sh` is intended to be run from within `testing/`, not from the repo root.
- The small test dataset is useful for fast end-to-end verification after code changes.
- If outputs already exist, the pipeline resumes unless `--force_rerun` is provided.

Useful command:

```bash
cd testing
bash run_quick_test.sh
```

## Known Gaps / Watch List

These items were observed during review and may still need cleanup:

1. Some docs and example scripts appear stale.
- Example scripts in `examples/` still reference older CLI flags and older contrast syntax.
- Some docs reference old output filenames or pre-workdir layouts.

2. The quick-test script assumptions are local.
- It is meant to run from `testing/`.
- Do not assume repo-root invocation is supported unless the script is updated.

3. edgeR emits this warning on the small dataset:
- `QL dispersion range: Inf - -Inf`
- This came from `fit$df.residual.zeros` on the test case.
- The run still completed successfully, but the warning may be worth revisiting if diagnostics are cleaned up later.

## Practical Guidance For Future Sessions

If a user reports a contrast issue:
- First check whether the contrast string is comma-delimited.
- Then check whether PSI output contains `delta_PSI`.
- Then confirm the `contrast` column in the result TSV matches the intended direction.

If a user wants more expression-like summaries in output:
- `util/run_edgeR_analysis.R` is the right place to add them.
- The result table is assembled there before writing `edgeR_results.intron_results.tsv`.

If a user wants PSI changes:
- `util/compute_psi.py` is the place to edit.

If a user wants output layout changes:
- `run_diff_splice_analysis.py` controls the final vs intermediate file placement.

## Recommended First Reads In A New Session

1. `README.md`
2. `docs/SESSION_MEMORY.md`
3. `run_diff_splice_analysis.py`
4. `util/run_edgeR_analysis.R`
5. `util/compute_psi.py`

## Current Branch Context

At the time this file was written:
- branch: `devel`
- recent commit of note: `6c9e462`

