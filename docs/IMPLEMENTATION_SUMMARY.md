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
- Offset (offset stat mode): exact `log(max_splice_plus_retained_depth)` after requiring `D >= 1` in every sample
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
- PSI/filtering denominator: `max_splice_plus_retained_depth`, identical to the default mode
- Model denominator: `D^svr` (see below)
- Offset (offset stat mode): exact `log(D^svr)` after requiring `D^svr >= 1` in every sample
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

Offset mode rejects an intron if its model denominator is zero in any biological
replicate or if `D - Y > 0` in fewer than `ceil(0.25 * n_samples)` replicates.
Zero denominators are not replaced with a pseudocount. The alternative-evidence
fraction is configurable with `--min_other_fraction`.

In `splice_plus_retained` mode, `logFC` approximates
`log2(reported_PSI_A / reported_PSI_B)`. In `splice_vs_rest` mode it instead
describes change relative to `D^svr`; reported PSI remains SPR-based. Dispersion
is estimated with `glmQLFit`; testing uses `glmQLFTest`.

### interaction

```bash
--stat_mode interaction
```

Focal/other interaction model. For each intron `i`, DSF uses:

```text
focal:  Y_i,s
other:  max(0, D_i,s - Y_i,s)
```

Choose the interaction engine with `--intx_engine`:

```bash
--intx_engine edgeR    # default: doubled count matrix, edgeR glmFit + glmLRT
--intx_engine DEXSeq   # optional: DEXSeq package with alternativeCountData
```

Both engines test a group-by-feature-type interaction by LRT. `logFC` is a log2
focal/other odds-ratio-like change. With an SPR model denominator, it is larger
in magnitude than offset-mode `logFC` for the same PSI shift when baseline PSI
is high. `delta_PSI` is computed from SPR depth before statistical testing and
is unaffected by `--stat_mode`.

Batch correction (`--batch_col`) adds the corresponding batch-by-feature-type
interaction to both full and reduced interaction models.

### Zero denominators and constitutively expressed introns

The two statistical modes handle zero or near-zero denominator evidence
differently, with important consequences for constitutively expressed introns.

**Offset mode — zero D is a hard failure.**
`log(D)` is supplied as a fixed numerical input. `log(0) = -∞`, so any intron
with `D = 0` in any replicate is dropped before testing. The
`--min_other_fraction` filter (default 25% of samples must have `D - Y > 0`)
is an additional guard: it removes introns where the denominator is not
meaningfully larger than the numerator in enough samples.

**Interaction mode — zero other-counts are valid NB observations.**
The NB GLM with a log link requires the fitted mean `μ > 0`, which the link
guarantees, not that every observed count `Y > 0`. A zero "other" count is a
valid draw from `NB(μ, φ)` for any positive `μ`; it contributes a
well-defined log-likelihood term and causes no singularity. edgeR tests via
`glmLRT` (likelihood ratio test), not a Wald test — the LRT compares
log-likelihoods of the full and reduced models and remains valid even when
separation pushes a coefficient toward a boundary. Only Wald standard errors
blow up under separation; the LRT p-value is unaffected. The offset mode's
`D = 0` and `--min_other_fraction` filters are not applied in interaction mode
(see `DSF.py:502-519`).

**Constitutive-in-one-group events.**
When an intron is constitutively used in group A (`D ≈ Y`, so `other ≈ 0`)
but partially skipped in group B (`other > 0`), the two modes behave
differently:

- Interaction mode correctly detects the event. The asymmetry — `other = 0` in
  group A, `other > 0` in group B — is exactly the group:exon_type interaction
  signal the model is designed to measure. Detection (LRT p-value and FDR)
  remains valid under separation; only Wald standard errors blow up, and edgeR
  uses `glmLRT`, not a Wald test. A large `logFC` reflects a genuinely large
  odds ratio (approaching infinity as `other → 0`), analogous to a zero-cell
  2×2 table — correct in sign and direction, not spurious. Only the magnitude
  is not precisely determined under separation; use `delta_PSI` for a bounded
  effect size.
- Offset mode's `--min_other_fraction` filter counts `D - Y > 0` samples
  across **all** samples. If the skipped group is a small fraction of total
  samples, the intron may not reach the threshold and is dropped before
  testing, losing a real event.

Truly constitutive in both groups (`other ≈ 0` everywhere) presents no
meaningful problem in either mode: there is no differential event, and both
modes return non-significant results.

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
- `util/run_edgeR_interaction_analysis.R`: edgeR stacked NB interaction model
- `util/run_DEXSeq_interaction_analysis.R`: DEXSeq focal/other interaction model
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
