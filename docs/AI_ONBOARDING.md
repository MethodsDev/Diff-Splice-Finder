# Intron-Focused Differential Splicing Analysis

## AI Onboarding / Design Context

Diff-Splice-Finder is a bulk RNA-seq differential splicing pipeline built around
intron-level evidence and edgeR GLMs with supplied offsets. The goal is to detect
differential intron usage while reducing confounding from local expression and
coverage changes.

## Current Model

The current mainline model is offset-mode based:

```text
log(expected intron count) = design effect + log(selected_denominator)
```

The default denominator mode is:

```text
exon_adjacent_depth = max(left_adjacent_depth, right_adjacent_depth)
```

Adjacent depths are exon-side read-depth windows outside the intron, computed by
`samtools depth`. `site_depth_offset` remains as a compatibility alias for
`max_adjacent_depth`.

The older donor/acceptor cluster-total denominator
`max(donor_cluster_total, acceptor_cluster_total)` is historical background and
still exists in utility code, but it is not the main `run_diff_splice_analysis.py`
path.

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

### Exon-Adjacent DSF

```bash
--offset_mode exon_adjacent_depth
```

This is the default mode. It uses read-level junction counts and
`max_adjacent_depth` for both PSI and the edgeR exposure.

### Splice-Plus-Retained DSF

```bash
--offset_mode splice_plus_retained
```

This experimental mode uses read-level junction counts and:

```text
max_splice_plus_retained_depth =
    max(left_splice_depth + left_retained_depth,
        right_splice_depth + right_retained_depth)
```

Splice depths are derived from canonical intron counts sharing each boundary.
Retained depths are intron-side read-depth windows computed by `samtools depth`.

### Gene-Median Splice-Plus-Retained DSF

```bash
--offset_mode gene_median_splice_plus_retained
--gtf annotation.gtf
```

This mode uses `max_splice_plus_retained_depth` for PSI but substitutes a
per-gene median of that denominator as the edgeR exposure. Introns without gene
assignment fall back to their own splice-plus-retained denominator.

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

Filtering happens before edgeR:

- canonical splice filter; noncanonical introns are not supported by the
  offset-mode refactor
- minimum total intron count
- minimum number of samples with nonzero intron counts
- minimum selected denominator depth in a configurable number of samples
- optional pre-edgeR `--min_delta_psi`

FDR is computed by edgeR over the tested post-filter intron set. The current
pipeline does not do post-edgeR delta-PSI filtering or FDR recomputation.

## Non-Goals

- transcript assembly
- isoform EM or abundance estimation
- exon-inclusion PSI from annotation-defined events
- single-cell modeling
- UMI deduplication

## First Files To Read

1. `README.md`
2. `examples/PARAMETER_GUIDE.md`
3. `run_diff_splice_analysis.py`
4. `util/depth_windows.py`
5. `util/site_depth.py`
6. `util/run_edgeR_analysis.R`
