# Session Memory

This file is a reusable brain dump for future work on `Diff-Splice-Finder`.

## Current Pipeline Model

- The pipeline is a single intron-level workflow.
- The current main path uses `--offset_mode`, not donor/acceptor cluster-total
  shared offsets.
- Default mode, `exon_adjacent_depth`:
  - numerator: read-level junction counts
  - PSI denominator: `max_adjacent_depth`
  - edgeR offset: `log(max_adjacent_depth + 0.5)`
- Experimental modes:
  - `splice_plus_retained`: PSI and edgeR exposure use
    `max_splice_plus_retained_depth`
  - `gene_median_splice_plus_retained`: PSI uses
    `max_splice_plus_retained_depth`, while edgeR uses the per-gene median
    exposure and requires `--gtf`
  - `splice_plus_retained_betabinom`: PSI uses
    `max_splice_plus_retained_depth`; the model uses focal/rest trials with a
    base-R quasibinomial GLM
- `--site_depth_strand_mode` controls stranded filtering for depth denominators.
  Junction discovery is not explicitly orientation-filtered, but canonical
  junctions are strand-resolved by splice motif.
- Retained-depth windows are shifted into the intron by
  `--retained_depth_inner_offset`, default `20`.

Main entrypoint:
- `run_diff_splice_analysis.py`

Key utility modules:
- `util/count_introns_from_bam.py`
- `util/depth_windows.py`
- `util/site_depth.py`
- `util/run_edgeR_analysis.R`

## Important Output Conventions

Final outputs go in the user-specified output directory:
- `edgeR_results.all.tsv`
- `edgeR_results.significant_introns.tsv`

Intermediate outputs go in:
- `<output_dir>/workdir/`

Important intermediate files:
- `workdir/introns_filtered.tsv`
- `workdir/site_depth_offsets.filtered.tsv`
- `workdir/edgeR_input.counts.tsv`
- `workdir/edgeR_input.offsets.tsv`
- `workdir/edgeR_input.annotations.tsv`
- `workdir/psi.psi_values.tsv`
- `workdir/edgeR_results.intron_results.tsv`
- `workdir/edgeR_results.intron_results_with_psi.tsv`
- `workdir/bam_inputs/intron_counts.matrix` and
  `workdir/bam_inputs/intron_counts.offsets.matrix` in BAM-manifest mode
- `workdir/bam_inputs/intron_counts.max_adjacent_depth.matrix`
- `workdir/bam_inputs/intron_counts.max_splice_plus_retained_depth.matrix`

## Contrast Contract

Use comma-delimited contrasts:
- Single control: `GroupA,GroupB`
- Pooled controls: `GroupA,GroupB1;GroupB2`

This format is now aligned across:
- `run_diff_splice_analysis.py`
- `util/run_edgeR_analysis.R`
- `util/compute_psi.py`

The streamlined main pipeline expects one explicit contrast per invocation.

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
  --offset_matrix test_site_depth_offsets.matrix \
  --samples test_metadata_control.tsv \
  --output_dir quick_test_contrast_review \
  --contrast 'perturb,control' \
  --min_intron_count 5 \
  --min_intron_samples 2 \
  --min_offset_depth 10 \
  --min_offset_samples 2 \
  --offset_mode exon_adjacent_depth \
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

- Prefer the local Makefile entrypoints:
  - `make -C testing test`
  - `make -C testing test_viz`
  - `make -C testing clean`
- The small test dataset is useful for fast end-to-end verification after code changes.
- If outputs already exist, the pipeline resumes unless `--force_rerun` is provided.

## Known Gaps / Watch List

These items were observed during review and may still need cleanup:

1. Some docs and example scripts appear stale.
- Historical docs under `docs/` still describe older shared-offset utilities.
- `--site_depth_offsets` remains as a deprecated alias for the default
  exon-adjacent matrix-mode denominator.

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

At the time this file was updated:
- branch: `dsf-mode-refactor`
- recent commits of note:
  - `3374320` refactored offset modes around read-depth denominators
  - `02cf293` added the tiny offset-mode fixture
