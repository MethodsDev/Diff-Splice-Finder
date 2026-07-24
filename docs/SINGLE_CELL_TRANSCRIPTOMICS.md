# Single-Cell Transcriptomics Compatibility

## Status

Diff-Splice-Finder (DSF) does not currently support cell-level statistical
inference. Its statistical model treats each matrix column as an independent
biological sample. Supplying individual cells as columns would therefore treat
cells from the same donor as independent replicates.

The recommended route for single-cell or single-nucleus RNA-seq is
**cell-type-resolved donor pseudobulk**. Cells are aggregated within each donor
and cell type, and the resulting donor-level aggregates are analyzed as DSF
samples. This retains cell-type specificity while preserving the biological
replicate structure required by edgeR and DEXSeq.

This document describes a proposed extension. It has not yet been implemented or
validated in DSF.

## What Can Be Reused

DSF models a focal intron count relative to a local or gene-scoped denominator:

```text
log(E[Y_i,s]) = X_s * beta_i + log(D_i,s)
```

where:

- `Y_i,s` is the focal junction count for intron `i` in sample `s`;
- `D_i,s` is the selected denominator;
- `X_s` contains the condition and optional batch terms.

This formulation does not require the source assay to be bulk RNA-seq. It does
require each column to represent an independent experimental unit with enough
coverage for dispersion estimation. A donor-level pseudobulk satisfies that
requirement when donors are the biological replicates.

The matrix-mode interface is the cleanest extension point. It already accepts:

- an intron-by-sample count matrix;
- a matching raw denominator matrix;
- sample metadata containing `sample_id` and the condition column.

A single-cell adapter can generate these files without changing the downstream
DSF model.

## Why Cells Cannot Be Used as Replicates

Cells from the same donor share genetic, environmental, preparation, and batch
effects. Treating them as independent samples underestimates uncertainty and
inflates the effective degrees of freedom. Adding donor as a batch covariate does
not generally repair this problem, particularly when condition is assigned at
the donor level.

The current implementation also conflicts with cell-level sparsity:

1. Offset mode requires `D > 0` in every selected sample. Per-cell junction and
   local-depth denominators are usually zero for most introns.
2. The edgeR interaction engine estimates median-ratio size factors from introns
   with positive counts in every column. Such introns may not exist in a
   single-cell matrix.
3. Zero-denominator observations are reported as PSI zero, which would conflate
   missing coverage with observed zero usage.
4. Python and R load dense matrices. The interaction model doubles the sample
   columns and includes a blocking coefficient for each sample, which is not
   suitable for thousands of cells.
5. BAM-mode counting is read based. It does not inspect cell barcodes or collapse
   junction evidence by UMI.

These constraints make direct cell input both statistically invalid for
multi-donor comparisons and operationally impractical.

## Recommended Pseudobulk Unit

Run one cell type or cell state at a time. Within that stratum, create one column
for every independent donor and condition combination.

For a case-control study in which each donor belongs to one condition:

```text
sample_id          donor_id    cell_type    group
D01_CD4            D01         CD4          case
D02_CD4            D02         CD4          case
D03_CD4            D03         CD4          control
D04_CD4            D04         CD4          control
```

The DSF columns are `D01_CD4`, `D02_CD4`, and so on. The group contrast is case
versus control. Donors, rather than cells, determine the replicate count.

For a paired or within-donor perturbation, construct donor-by-condition
pseudobulks:

```text
sample_id              donor_id    cell_type    group
D01_CD4_control        D01         CD4          control
D01_CD4_perturb        D01         CD4          perturb
D02_CD4_control        D02         CD4          control
D02_CD4_perturb        D02         CD4          perturb
```

The offset engine may be able to use `donor_id` through `--batch_col` as a fixed
blocking factor in this balanced design. DSF does not currently validate design
rank or distinguish technical batch from donor blocking, so paired designs need
specific tests before they are declared supported.

Do not pool all donors into one cell-type aggregate. That produces one
observation per condition and removes the biological variance needed for
replicate-aware inference. Randomly splitting cells into pseudoreplicates also
does not create independent biological samples.

## Aggregation Invariants

### Preserve donor identity

Aggregate cells only within the same donor, cell type, and modeled condition.
Record the number of contributing cells and the aggregation criteria for quality
control, but do not use cell count as the replicate count.

### Sum counts, not PSI

For each donor-level pseudobulk, sum focal junction evidence and denominator
components. Compute PSI after aggregation:

```text
Y_pseudobulk = sum(Y_cell)
PSI_pseudobulk = Y_pseudobulk / D_pseudobulk
```

Averaging per-cell PSI would give unstable, low-coverage cells the same weight as
well-covered cells and would require a defensible missing-coverage rule.

### Compute the maximum after aggregation

The default splice-plus-retained denominator is:

```text
D = max(left_splice_depth + left_retained_depth,
        right_splice_depth + right_retained_depth)
```

Aggregate the primitive left and right components first, then compute `max`:

```text
left_pseudobulk  = sum(left_cell)
right_pseudobulk = sum(right_cell)
D_pseudobulk     = max(left_pseudobulk, right_pseudobulk)
```

Do not sum denominators that were already calculated per cell, because in
general:

```text
sum(max(left_cell, right_cell)) != max(sum(left_cell), sum(right_cell))
```

Filtering a donor and cell-type BAM before running DSF would preserve this
ordering automatically. A matrix adapter must implement it explicitly.

### Keep numerator and denominator units coherent

Current DSF junction counts are read-level split-alignment counts, while retained
components come from read depth. For UMI assays, a new adapter must define how
junction UMIs and local denominator evidence represent the same molecular
exposure. Mixing UMI-deduplicated numerators with non-deduplicated read-depth
denominators would change the estimand and may create protocol-dependent effects.

The adapter should record at least:

- count unit, such as reads, fragments, or UMIs;
- donor, condition, cell type, and sequencing batch;
- cell inclusion rules and contributing cell count;
- denominator construction method;
- reference annotation and alignment source.

## Protocol Considerations

### Full-length single-cell RNA-seq

Smart-seq2 and Smart-seq3 are the strongest candidates for donor-level DSF
pseudobulk. Their broader transcript coverage makes split junctions and local
depth windows more informative. PCR amplification and molecule recovery still
need attention, particularly when read counts rather than UMIs are used.

### Droplet 3-prime or 5-prime RNA-seq

10x-style assays have sparse, position-biased junction coverage. DSF may retain
power for abundant genes and events near the captured transcript end, but a
negative result away from that end may only indicate lack of information.

Neither `splice_plus_retained` nor `splice_vs_rest` should be assumed to remove
this bias. Both modes should be evaluated against matched bulk, full-length, or
targeted measurements. Coverage thresholds will need to be tuned at the donor
pseudobulk level rather than inherited unchanged from bulk defaults.

### Single-nucleus RNA-seq

Nuclear libraries contain abundant unspliced and partially processed RNA. The
retained-depth component may therefore measure pre-mRNA abundance or nuclear RNA
processing rather than regulated intron retention. A junction-only denominator
could be more appropriate for some analyses, but DSF does not currently provide
such a mode. Any new denominator needs simulation and orthogonal validation.

### Single-cell long-read RNA-seq

Long reads can supply informative junction chains, but cells should still be
aggregated by donor and cell type for condition-level inference. UMI and molecule
handling must prevent multiple reads from the same molecule from contributing
as independent junction evidence.

## Possible Cell-Level Backend

Studies of cell-to-cell splicing heterogeneity, continuous cell state, or
within-donor perturbations may require a true cell-level model. Such a model
would be a separate backend rather than a switch on the current edgeR analysis.
Possible formulations include:

```text
log(E[Y_i,d,c]) = log(D_i,d,c) + X_d,c * beta_i + b_i,d
```

where `b_i,d` is a donor random effect, or a beta-binomial mixed model for focal
versus alternative evidence when `0 <= Y <= D` is guaranteed. A hurdle model
could separate the probability of observing informative molecules from the
conditional intron-usage change.

A cell-level backend would also require:

- sparse matrix input and processing;
- barcode-aware and UMI-aware counting;
- explicit donor nesting and repeated-measure support;
- a missing-coverage state distinct from PSI zero;
- cell-level technical covariates;
- donor-level permutation or mixed-model inference;
- methods for continuous trajectories if group contrasts are insufficient.

The present DSF dispersion estimates and FDR calculations cannot be reused with
cells treated as independent observations.

## Proposed Implementation Sequence

1. Add strict validation for unique samples, finite nonnegative counts,
   denominator alignment, replicate counts, design rank, and count units.
2. Add a donor-by-cell-type pseudobulk adapter that writes DSF matrix-mode inputs
   and aggregation provenance.
3. Add tests showing that pooled alignments and component-wise matrix aggregation
   produce the same numerator and denominator matrices.
4. Benchmark `splice_plus_retained` and `splice_vs_rest` separately for
   full-length, droplet, and nuclear assays.
5. Calibrate count and denominator filters using null simulations and matched
   truth data. Check type I error at the donor level, power, effect-size bias, and
   sensitivity to unequal cell numbers.
6. Add paired-donor design support with explicit formula and rank validation.
7. Consider a separate cell-level backend only if the scientific question cannot
   be answered by donor pseudobulk.

## Interpretation

A pseudobulk extension would make DSF a cell-type-resolved differential splicing
method for single-cell data. It would not estimate reliable PSI for every cell or
establish cell-to-cell splicing heterogeneity. Results would describe condition
or cell-type differences among independent donors after aggregating the cells in
each analysis stratum.

## References

1. Zimmerman KD, Espeland MA, Langefeld CD. A practical solution to
   pseudoreplication bias in single-cell studies. *Nature Communications*.
   2021;12:738. <https://doi.org/10.1038/s41467-021-21038-1>
2. Lee H, Han B. Pseudobulk with proper offsets has the same statistical
   properties as generalized linear mixed models in single-cell case-control
   studies. *Bioinformatics*. 2024;40:btae498.
   <https://doi.org/10.1093/bioinformatics/btae498>
3. Kubota N, Chen L, Zheng S. Shiba: a versatile computational method for
   systematic identification of differential RNA splicing across platforms.
   *Nucleic Acids Research*. 2025;53:gkaf098.
   <https://doi.org/10.1093/nar/gkaf098>
4. Buen Abad Najar CF, Yosef N, Lareau LF. Coverage-dependent bias creates the
   appearance of binary splicing in single cells. *eLife*. 2020;9:e54603.
   <https://doi.org/10.7554/eLife.54603>
