# Single-Cell Transcriptomics Compatibility

## Status

Diff-Splice-Finder (DSF) does not currently support cell-level or cell-pool
statistical inference. Its statistical model treats each matrix column as an
independent sample. Supplying individual cells as columns would therefore treat
cells from the same donor as independent biological replicates.

For population-level comparisons, the recommended route is
**cell-type-resolved donor pseudobulk**. Cells are aggregated within each donor
and cell type, and the resulting donor-level aggregates are analyzed as DSF
samples. This retains cell-type specificity and the biological replicate
structure required by edgeR and DEXSeq.

A single donor supports a narrower analysis: conditional comparison of cell
types or states within the observed specimen. Section
"Single-Donor, Within-Specimen Analysis" describes a proposed micro-pseudobulk
mode for that purpose. It must not be interpreted as population-level inference.

Neither extension has been implemented or validated in DSF.

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

These constraints make direct cell input incompatible with the current DSF
backends, even when the intended inference is limited to one donor.

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

Do not pool all donors into one cell-type aggregate for a population-level
comparison. That produces one observation per condition and removes the
biological variance needed for replicate-aware inference. Randomly splitting
cells also does not create independent biological samples. Such splits can,
however, support the conditional single-donor analysis described below.

## Single-Donor, Within-Specimen Analysis

### Estimand

A one-donor analysis can address whether intron usage differs between the
sampled cell populations in that specimen. It cannot establish that the same
difference recurs across donors or that a condition causes the difference in a
population. Tissue niche, lineage, cell cycle, processing batch, and cell-type
classification can all contribute to the observed association.

The analysis should therefore use language such as:

> Within this donor and these recovered cells, intron usage differed between
> the annotated cell populations and was stable across cell-pool partitions.

It should not describe cell pools as biological replicates or report the result
as a population-level cell-type effect.

### Precedents in single-cell splicing methods

Published single-cell methods use several strategies for sparse splicing data:

| Method | Observation and model | Relevance to DSF |
| --- | --- | --- |
| `scASfind` | Randomly combines cells of the same type into small pools, calculates event PSI per pool, and searches for cell-type-specific patterns. The default is five cells per pool; confidently quantified nodes require at least 10 reads. | The closest precedent for DSF micro-pseudobulks. Pooling improves junction coverage but does not create biological replication. |
| `Sierra` | Creates artificial pseudobulk profiles from UMI peak counts and fits DEXSeq. | Shows that an artificial-pool count GLM is computationally feasible for single-cell feature usage. Its main target is polyadenylation and transcript-end usage rather than general intron usage. |
| `satuRn` | Fits a quasi-binomial feature-usage GLM directly to cells and moderates overdispersion across transcripts. Its single-cell workflow includes empirical-null recalibration because the theoretical null can be liberal. | Demonstrates that within-donor count-GLM inference is possible, while showing that single-cell calibration cannot be assumed from a bulk model. |
| `SCATS` | Fits a hierarchical model to per-cell inclusion and exclusion evidence, with terms for capture, amplification, dropout, and group inclusion level. | Models single-cell technical noise directly instead of manufacturing count-library replicates. |
| `BRIE2` | Uses a multinomial likelihood and logit-normal regression for latent per-cell PSI. Cell type, state, or pseudotime can enter as covariates. | Suited to annotated, predominantly two-isoform events in full-length data. Its ELBO gain is model-selection evidence, not a frequentist p-value. |
| `scQuint` | Fits a Dirichlet-multinomial GLM to sparse alternative-intron counts that share a splice site. | A close compositional alternative for full-length data, although the reported p-values still condition on the sampled cells when only one donor is present. |
| `SpliZ` | Computes a gene-level junction-deviation score for each cell and uses permutation tests for group or state associations. | Better suited than current DSF to UMI droplet data and complex local splicing, but it tests a gene-level score rather than one DSF intron. |
| `MARVEL` | Tests per-cell PSI distributions for plate data. For droplet data it pools cells by population and uses cell-label permutation on pooled junction usage. | Provides a precedent for within-donor pooled junction analysis and for separating plate and droplet workflows. |

These methods do not estimate donor-to-donor variance from one donor. Their
single-donor results are conditional on the observed cells, regardless of
whether the software reports a p-value, empirical score, Bayes factor, or pool
enrichment statistic.

DSF also tests intron usage rather than transcript usage. Comparisons with
`satuRn`, `BRIE2`, or transcript-level DTU tools concern the statistical strategy,
not an identical biological estimand.

### Proposed DSF micro-pseudobulk workflow

The most compatible design for the existing edgeR backend is a set of disjoint
cell pools within each annotated cell type or state:

```text
sample_id       donor_id    cell_type    pool_id    group
D01_CD4_P01     D01         CD4          P01        CD4
D01_CD4_P02     D01         CD4          P02        CD4
D01_B_P01       D01         B             P01        B
D01_B_P02       D01         B             P02        B
```

Construct one partition without replacement for the primary analysis. Use the
same pool size in both groups and keep cells from different library batches or
processing strata separate unless the design accounts for those strata. Do not
form pools using the tested splicing measurements.

For deep full-length assays, five cells per pool is a reasonable starting point
based on `scASfind`, not a DSF default. Increase the pool size until the selected
denominators are positive and adequately covered for the retained introns. At
least three usable pools per group are needed to fit a replicated design; five
or more gives a more useful diagnostic sample. These are pilot criteria, not a
claim that the pools are independent biological units.

For each pool:

1. Sum focal junction counts across its cells.
2. Sum the primitive left and right denominator components across its cells.
3. Calculate the selected denominator after aggregation, including the `max`
   operation used by `splice_plus_retained`.
4. Write one matrix column and one metadata row for the pool.

Use the edgeR offset engine for an initial implementation. The interaction
engine requires positive counts across all columns for its size-factor
calculation and creates two modeled columns plus a blocking coefficient per
pool, making it a poor first choice for sparse micro-pseudobulks. Run pre-test
coverage filters, but set `--min_delta_psi 0` during model fitting so that an
effect-size threshold does not also select the tested hypotheses.

Designate one pool size and one random partition before testing. Repeat the
analysis across additional seeds and pool sizes as a stability assessment.
Report sign agreement, effect-size variation, rank stability, and the fraction
of partitions in which an event passes coverage and significance criteria. The
repeated partitions reuse cells, so their p-values are dependent and must not be
combined as independent evidence.

### Interpretation of edgeR results

The fitted edgeR dispersion would describe variation among cell aggregates from
this donor. It would not include between-donor biological variation. Random
splitting can improve coverage and measure sensitivity to the sampled cells,
but increasing the number of pools cannot increase the biological sample size
above one.

Until null calibration has been completed, report edgeR p-values and FDR as
exploratory, donor-conditional statistics. Fixed-dispersion or no-replicate
edgeR analyses are weaker alternatives because the data cannot estimate the
required variability.

### Validation gate

Before DSF exposes this mode, validation should cover:

1. Randomly split one annotated cell type into pseudo-groups, then run the full
   adapter and DSF workflow. Repeat this null experiment across seeds, pool
   sizes, coverage levels, and unequal group sizes to measure false-positive
   rates.
2. Simulate known intron-usage changes while preserving realistic cell-level
   sparsity and denominator structure. Measure power and effect-size bias.
3. Downsample cells and molecules to test whether discoveries depend on a small
   number of high-coverage cells.
4. Compare full-length results with `scQuint`, `SCATS`, or `BRIE2`; compare UMI
   droplet results with `SpliZ` or pooled-junction methods such as `MARVEL`.
5. Confirm selected events with targeted RT-PCR, matched bulk or full-length
   sequencing, or another independent measurement.

Nominal FDR should only be retained if the null experiments show acceptable
calibration. Otherwise, the mode should report effect sizes and partition
stability without inferential FDR claims.

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
3. Test pooled alignments against component-wise matrix aggregation for exact
   agreement of numerator and denominator matrices.
4. Validate the donor-pseudobulk workflow across full-length, droplet, and
   nuclear assays before making population-level claims.
5. Add an explicitly named single-donor cell-pool mode. Keep its metadata,
   output labels, and documentation separate from donor-pseudobulk inference.
6. Run the null calibration, simulation, downsampling, and cross-method checks
   listed above. Do not expose inferential FDR unless those checks support it.
7. Add paired-donor design support with explicit formula and rank validation.
8. Consider a separate cell-level backend only if the scientific question needs
   cell-to-cell heterogeneity, continuous state, or covariates that pooling would
   discard.

## Interpretation

A donor-pseudobulk extension would make DSF a cell-type-resolved differential
splicing method for replicated single-cell studies. A separate cell-pool mode
could identify donor-conditional differences between cell types or states in
one specimen. Neither mode would estimate reliable PSI for every cell or test
cell-to-cell splicing heterogeneity.

Population-level claims require independent donors. Single-donor results should
emphasize effect size, direct junction support, denominator coverage, and
stability across cell-pool partitions.

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
5. Hu Y, Wang K, Li M. Detecting differential alternative splicing events in
   scRNA-seq with or without unique molecular identifiers. *PLOS Computational
   Biology*. 2020;16:e1007925.
   <https://doi.org/10.1371/journal.pcbi.1007925>
6. Huang Y, Sanguinetti G. BRIE2: computational identification of splicing
   phenotypes from single-cell transcriptomic experiments. *Genome Biology*.
   2021;22:251. <https://doi.org/10.1186/s13059-021-02461-5>
7. Benegas G, Fischer J, Song YS, Eyras E. Robust and annotation-free analysis
   of alternative splicing across diverse cell types in mice. *eLife*.
   2022;11:e73520. <https://doi.org/10.7554/eLife.73520>
8. Olivieri JE et al. RNA splicing programs define tissue compartments and cell
   types at single-cell resolution. *eLife*. 2021;10:e70692.
   <https://doi.org/10.7554/eLife.70692>
9. Song Y, Parada G, Lee JTH, Hemberg M. Mining alternative splicing patterns in
   scRNA-seq data using scASfind. *Genome Biology*. 2024;25:196.
   <https://doi.org/10.1186/s13059-024-03323-6>
10. Gilis J, Vitting-Seerup K, Van den Berge K, Clement L. satuRn: scalable
    analysis of differential transcript usage for bulk and single-cell
    RNA-sequencing applications. *F1000Research*. 2022;10:374.
    <https://doi.org/10.12688/f1000research.51749.2>
11. Wen WX, Mead AJ, Thongjuea S. MARVEL: an integrated alternative splicing
    analysis platform for single-cell RNA sequencing data. *Nucleic Acids
    Research*. 2023;51:e29. <https://doi.org/10.1093/nar/gkac1260>
12. Patrick R, Humphreys DT, Janbandhu V, Oshlack A, Ho JWK, Harvey RP, Lo KK.
    Sierra: discovery of differential transcript usage from polyA-captured
    single-cell RNA-seq data. *Genome Biology*. 2020;21:167.
    <https://doi.org/10.1186/s13059-020-02071-7>
