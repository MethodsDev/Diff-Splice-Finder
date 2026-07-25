# Session Memory

This file is a reusable brain dump for future work on `Diff-Splice-Finder`.

## Current Pipeline Model

- The pipeline is a single intron-level workflow.
- The current main path uses `--offset_mode`, not donor/acceptor cluster-total
  shared offsets.
- Default mode, `splice_plus_retained`:
  - numerator: read-level junction counts
  - PSI denominator: `max_splice_plus_retained_depth`
  - edgeR offset: exact `log(max_splice_plus_retained_depth)`; requires `D >= 1` in every replicate and `D - Y > 0` in at least 25% of replicates by default
- Gene-scoped mode, `splice_vs_rest`:
  - requires `--gtf`
  - computes a gene-total junction denominator plus the local
    `splice_plus_retained` remainder
  - introns without a gene assignment fall back to `splice_plus_retained`
- Statistical modes:
  - `offset`: edgeR NB QL GLM with supplied log denominator offsets
  - `interaction`: DEXSeq-style stacked NB interaction model
- `--site_depth_strand_mode` controls stranded filtering for depth denominators.
  Junction discovery is not explicitly orientation-filtered, but canonical
  junctions are strand-resolved by splice motif.
- Retained-depth windows are shifted into the intron by
  `--retained_depth_inner_offset`, default `20`.

Main entrypoint:
- `DSF.py`

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
- `workdir/bam_inputs/intron_counts.matrix`
- `workdir/bam_inputs/intron_counts.max_adjacent_depth.matrix`
- `workdir/bam_inputs/intron_counts.max_splice_plus_retained_depth.matrix`

## Contrast Contract

Use comma-delimited contrasts:
- Single control: `GroupA,GroupB`
- Pooled controls: `GroupA,GroupB1;GroupB2`

This format is now aligned across:
- `DSF.py`
- `util/run_edgeR_analysis.R`
- `util/compute_psi.py`

The streamlined main pipeline expects one explicit contrast per invocation.

## Recent Fixes

These were implemented and pushed in commit `6c9e462`:

1. Fixed contrast format drift.
- The all-pairwise path in `DSF.py` had been generating `A-B`.
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
../DSF.py \
  --matrix mode_refactor_inputs/intron_counts.matrix \
  --offset_matrix mode_refactor_inputs/intron_counts.max_splice_plus_retained_depth.matrix \
  --samples mode_refactor_inputs/sample_metadata.tsv \
  --output_dir quick_test_contrast_review \
  --contrast 'A,B' \
  --min_intron_count 1 \
  --min_intron_samples 1 \
  --min_offset_depth 1 \
  --min_offset_samples 1 \
  --offset_mode splice_plus_retained \
  --fdr_threshold 1 \
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

## Known Notes / Watch List

These items are useful context for future maintenance:

1. Some historical docs are intentionally retained.
- `docs/SHARED_OFFSETS.md` describes the older shared-offset strategy and starts
  with a historical note pointing readers to current docs.
- `--site_depth_offsets` is no longer supported; matrix mode uses
  `--offset_matrix`.

2. The quick-test script resolves local paths.
- It resolves paths relative to its own location and can be run from the repo
  root or from `testing/`.

3. edgeR emits this warning on the small dataset:
- `QL dispersion range: Inf - -Inf`
- This came from `fit$df.residual.zeros` on the test case.
- The run still completed successfully, but the warning may be worth revisiting if diagnostics are cleaned up later.

4. Zero denominators and constitutive introns behave differently by stat mode.
- Offset mode: `log(D)` is a fixed input; `D = 0` → `log(0)` = undefined → hard drop.
  `--min_other_fraction` also drops introns where `D - Y > 0` holds in fewer than
  25% of samples, which can preemptively remove a real constitutive→skipped event
  when the changed group is a small fraction of total samples.
- Interaction mode: `other = max(0, D - Y)` zeros are valid NB observations.
  edgeR uses `glmLRT` (LRT, not Wald) — detection (LRT p-value and FDR)
  remains valid under separation; only Wald SEs blow up. A large `logFC` is
  correct in sign and direction (true odds ratio → ∞, like a zero-cell 2×2
  table), not spurious. Only the magnitude is not precisely determined under
  separation; use `delta_PSI` for a bounded effect size. No preemptive drop.
- See `docs/IMPLEMENTATION_SUMMARY.md` "Zero denominators and constitutively
  expressed introns" for the full explanation.

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
- `DSF.py` controls the final vs intermediate file placement.

## Recommended First Reads In A New Session

1. `README.md`
2. `docs/SESSION_MEMORY.md`
3. `DSF.py`
4. `util/run_edgeR_analysis.R`
5. `util/compute_psi.py`

## Current Branch Context

At the time this file was updated:
- branch: `master`
- recent commits of note:
  - `52cea34` fixed interaction mode edgeR design handling
  - `ea06f05` renamed the driver script to `DSF.py`
