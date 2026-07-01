# Intron-Focused Differential Splicing Analysis

## AI Onboarding / Design Context

Diff-Splice-Finder is a bulk RNA-seq differential splicing pipeline built around
intron-level evidence and edgeR GLMs with supplied offsets. The goal is to detect
differential intron usage while reducing confounding from local expression and
coverage changes.

## Current Model

The current mainline model is site-depth based:

```text
log(expected intron count) = design effect + log(selected_denominator)
```

The default denominator is:

```text
site_depth_offset = max(donor_site_window_depth, acceptor_site_window_depth)
```

The depth window includes aligned reference blocks around each splice-site
coordinate. It captures local coverage available at either end of the intron, so
edgeR tests whether the intron count changes relative to that local coverage.

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

### Default DSF

```bash
--count_unit read
--psi_denominator_mode site_depth
--test_offset_mode site_depth
```

This is the normal production mode. It can run from precomputed matrices or from
a BAM manifest.

### Strict-local DSF

```bash
--count_unit fragment
--psi_denominator_mode strict_local_depth
--test_offset_mode strict_local_depth
```

This experimental mode uses focal fragment junction counts and a strict local
splice-decision denominator:

```text
strict_local_depth = max(donor_decision_depth, acceptor_decision_depth)
decision_depth = split fragments sharing the boundary + unspliced fragments crossing the boundary
```

### Gene-median-strict DSF

```bash
--count_unit fragment
--psi_denominator_mode strict_local_depth
--test_offset_mode gene_median_strict_depth
--gtf annotation.gtf
```

This mode uses strict-local PSI but substitutes a per-gene median strict depth
as the edgeR exposure. It is intended as a stabilization experiment. In this
mode, `delta_PSI` is a strict-local PSI difference, while `logFC` is the GLM
coefficient under the gene-median exposure.

Strict modes require BAM-manifest input because focal fragment counts and strict
depths are generated from BAM files.

## Strandedness

`--site_depth_strand_mode` controls stranded filtering for site-depth offsets:

- `unstranded`
- `F` or `R` for single-end protocols
- `FR` or `RF` for paired-end protocols

Junction discovery itself does not explicitly filter split reads by read
orientation; canonical junction strand is inferred from the splice motif. This
usually makes the numerator effectively strand-resolved, but denominator depth
can still be contaminated by antisense or overlapping local coverage, so
stranded site-depth offsets matter for stranded libraries.

The strict-depth helper currently falls back to genomic/unstranded counting if a
stranded mode is supplied.

## Filtering

Filtering happens before edgeR:

- canonical splice filter unless `--keep_noncanonical` is supplied
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
4. `util/site_depth.py`
5. `util/strict_splice_depth.py`
6. `util/run_edgeR_analysis.R`
