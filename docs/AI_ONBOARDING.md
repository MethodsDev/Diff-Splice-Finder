# Intron-Focused Differential Splicing Analysis

## AI Onboarding / Design Context

Diff-Splice-Finder is a bulk RNA-seq differential splicing pipeline built around
intron-level evidence and edgeR GLMs with supplied offsets. The goal is to detect
differential intron usage while reducing confounding from local expression and
coverage changes.

## Current Model

Two independent axes of choice: **offset mode** (denominator) and **statistical
mode** (test engine).

### Offset modes

The default denominator is:

```text
max_splice_plus_retained_depth =
    max(left_splice_depth + left_retained_depth,
        right_splice_depth + right_retained_depth)
```

Splice depths are sums of canonical intron counts sharing each boundary.
Retained depths are intron-interior read-depth windows computed by `samtools depth`.

The older donor/acceptor cluster-total denominator
`max(donor_cluster_total, acceptor_cluster_total)` is historical background and
no longer supported in `DSF.py`.

## Feature Definition

Each feature is an observed intron/junction:

```text
(chrom, intron_start, intron_end, inferred_strand)
```

Strand is inferred from splice-site motifs. Canonical plus-strand motifs include
`GT--AG`, `GC--AG`, and `AT--AC`; canonical minus-strand motifs appear on the
reference as their reverse-complement pairs.

Default junction counts are read-level split-junction counts. A long read that
spans multiple introns contributes one count to each intron it supports.

## Supported Modalities

### Splice-Plus-Retained (default)

```bash
--offset_mode splice_plus_retained
```

Uses read-level junction counts and `max_splice_plus_retained_depth` for both
PSI and the edgeR/interaction denominator.

### Splice-vs-Rest (requires `--gtf`)

```bash
--offset_mode splice_vs_rest --gtf annotation.gtf
```

Extends `splice_plus_retained` with gene-total junction evidence:

```text
D^svr_i,s = gene_total_junction_count_i,s
          + max(0, max_splice_plus_retained_depth_i,s - Y_i,s)
```

Gene totals are computed in-pipeline from the count matrix + GTF. Introns with
no gene assignment fall back to `splice_plus_retained` automatically.

## Statistical Modes (`--stat_mode`)

### offset (default)

```bash
--stat_mode offset
```

Fixed log-offset edgeR NB QL GLM:

```text
log(E[Y_i,s]) = group_s + log(D_i,s)
```

`logFC` ≈ `log2(PSI_A / PSI_B)`.

### interaction

```bash
--stat_mode interaction
```

Focal/other interaction model. For each intron, DSF compares focal (`Y_i,s`) to
other (`max(0, D_i,s - Y_i,s)`). `--intx_engine edgeR` is the default and uses
the current edgeR stacked NB LRT. `--intx_engine DEXSeq` uses the DEXSeq package
with `alternativeCountData`. `logFC` is a log2 focal/other odds-ratio-like
effect—not a PSI log-ratio. `delta_PSI` is unaffected by `--stat_mode`.

## Strandedness

`--site_depth_strand_mode` controls stranded filtering for depth denominators:

- `unstranded`
- `F` or `R` for single-end protocols
- `FR` or `RF` for paired-end protocols

Junction discovery itself does not explicitly filter split reads by read
orientation; canonical junction strand is inferred from the splice motif. This
usually makes the numerator effectively strand-resolved, but denominator depth
can still be contaminated by antisense or overlapping local coverage, so
stranded depth modes matter for stranded libraries.

## Filtering

Filtering happens before statistical testing:

- canonical splice filter; noncanonical introns are not supported by the
  offset-mode refactor
- minimum total intron count
- minimum number of samples with nonzero intron counts
- minimum selected denominator depth in a configurable number of samples
- optional pre-test `--min_delta_psi`

FDR is computed over the tested post-filter intron set. The current pipeline
does not do post-test delta-PSI filtering or FDR recomputation.

## Non-Goals

- transcript assembly
- isoform EM or abundance estimation
- exon-inclusion PSI from annotation-defined events
- single-cell modeling
- UMI deduplication

## First Files To Read

1. `README.md`
2. `examples/PARAMETER_GUIDE.md`
3. `DSF.py`
4. `util/depth_windows.py`
5. `util/site_depth.py`
6. `util/run_edgeR_analysis.R`
7. `util/run_edgeR_interaction_analysis.R`
8. `util/run_DEXSeq_interaction_analysis.R`
