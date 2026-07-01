# Implementation Summary

## Overview

Diff-Splice-Finder is currently a single intron-level edgeR pipeline. Each
intron is tested once with a negative-binomial GLM and a user-selected exposure
offset. The main production mode uses splice-site depth denominators rather than
donor/acceptor cluster-total offsets.

The older `max(donor_cluster_total, acceptor_cluster_total)` implementation is
kept in utility code and historical docs, but it is not the main analysis path
used by `run_diff_splice_analysis.py`.

## Current Modes

### DSF default

```bash
--count_unit read
--psi_denominator_mode site_depth
--test_offset_mode site_depth
```

- Numerator: read-level split-junction counts from `count_introns_from_bam.py`
- PSI denominator: `site_depth_offset`
- edgeR offset: `log(site_depth_offset + 0.5)`
- Denominator definition:

```text
site_depth_offset = max(donor_site_window_depth, acceptor_site_window_depth)
```

Depth is computed over aligned reference blocks in windows around the two splice
site coordinates. Paired-end overlaps are collapsed so an overlapping mate pair
does not contribute twice to the same depth position.

### DSF strict-local

```bash
--count_unit fragment
--psi_denominator_mode strict_local_depth
--test_offset_mode strict_local_depth
```

- Numerator: focal fragment junction counts
- PSI denominator: `strict_local_depth`
- edgeR offset: `log(strict_local_depth + 0.5)`
- Denominator definition:

```text
strict_local_depth = max(donor_decision_depth, acceptor_decision_depth)
decision_depth = split fragments sharing the boundary + unspliced fragments crossing the boundary
```

### DSF gene-median-strict

```bash
--count_unit fragment
--psi_denominator_mode strict_local_depth
--test_offset_mode gene_median_strict_depth
--gtf annotation.gtf
```

- Numerator: focal fragment junction counts
- PSI denominator: `strict_local_depth`
- edgeR offset: `log(per-gene median strict_local_depth + 0.5)`
- Fallback: introns without a gene assignment use their local strict depth

In this mode, `delta_PSI` remains a strict-local PSI difference, while `logFC`
is the edgeR coefficient under the gene-median exposure. It should not be read
as the exact log-ratio of the reported strict-local mean PSI values.

## Input Modes

### Matrix mode

Matrix mode requires:

- `--matrix`: intron x sample count matrix
- `--site_depth_offsets`: intron x sample raw site-depth denominator matrix
- `--samples`: metadata with `sample_id` and the configured group column

Strict fragment-depth modes are not available in matrix mode because the strict
focal fragment counts and strict local depths are computed directly from BAMs.

### BAM-manifest mode

When `--matrix` is omitted, `--samples` is interpreted as a BAM manifest with:

```text
sample_type    replicate_id    bam_file
```

The pipeline runs a discovery pass, a targeted pass for complete site-depth
offsets, builds matrices under `workdir/bam_inputs/`, writes downstream
`sample_id/group` metadata, and then runs edgeR. Strict-depth matrices are
generated only when a strict mode is requested.

## Filtering and FDR

Filtering happens before edgeR:

- canonical splice filter unless `--keep_noncanonical` is supplied
- total count and nonzero-sample filters
- selected denominator depth filter via `--min_offset_depth` and `--min_offset_samples`
- optional pre-edgeR `--min_delta_psi` filter

edgeR's FDR is computed over the introns that pass these prefilters. The current
pipeline does not recompute FDR after edgeR. Final significant results are rows
marked significant by edgeR using the configured `--fdr_threshold` and
`--min_logFC`; if `--min_delta_psi` was nonzero, all tested rows already passed
that threshold.

## Strandedness

`--site_depth_strand_mode` can be `unstranded`, `F`, `R`, `FR`, or `RF`.

In default DSF mode this option controls strand filtering for site-depth
denominators. Junction counts are discovered from split alignments without an
explicit read-orientation filter, then annotated by splice motif. Canonical
splice junctions are usually strand-resolved by their reference dinucleotides,
but overlapping antisense transcription can still matter for local coverage,
which is why stranded site-depth offsets are useful.

The strict-depth helper currently counts strict fragment depths in
genomic/unstranded mode. If a stranded `--site_depth_strand_mode` is passed while
strict modes are enabled, strict-depth counting warns and falls back to
unstranded behavior.

## Core Files

- `run_diff_splice_analysis.py`: main orchestrator, filtering, mode selection,
  final output promotion
- `util/count_introns_from_bam.py`: read-level junction counts and site-depth
  offset reporting
- `util/site_depth.py`: site-depth denominator calculation and stranded depth
  handling
- `util/strict_splice_depth.py`: focal fragment counts and strict local
  splice-decision depths
- `util/run_edgeR_analysis.R`: edgeR GLM with supplied log offsets
- `util/build_intron_count_matrix.py`: count and site-depth matrix construction

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
    ├── intron_counts.offsets.matrix
    ├── focal_fragment_counts.matrix          # strict modes only
    └── strict_local_depth.matrix             # strict modes only
```
