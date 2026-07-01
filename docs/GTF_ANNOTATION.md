# GTF Annotation and PSI Notes

## Overview

The pipeline supports optional GTF annotation to enrich intron-level results with:

- `gene_name`
- `intron_status`
- `overlapping_genes`

The same GTF can also be used by the DEXSeq-like plotting utility to draw transcript structure for a target gene.

## Main Pipeline Usage

```bash
python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/my_analysis \
    --gtf /path/to/annotation.gtf
```

When `--gtf` is provided, annotation is added during the main pipeline run and propagated into the edgeR result tables.

## Result Columns

When annotation is available, these columns may appear in the results:

1. `gene_name`
2. `intron_status`
   - `known`
   - `novel`
   - `unknown`
3. `overlapping_genes`

`known` means the intron exactly matches an intron derived from the GTF. `novel` means the intron was observed in the data but not found as an exact annotated intron.

## Intron ID Formats

Supported intron IDs include:

1. `chr:start-end:strand`
2. `chr:start-end^DONOR--ACCEPTOR^STATUS`

For splice-motif IDs, strand is inferred from the motif.

## GTF Caching

Parsed annotation is cached alongside the input GTF as:

`<gtf_basename>.intron_cache.tsv`

The cache is reused if it is newer than the GTF and regenerated automatically if the GTF changes.

## PSI Outputs

The pipeline computes PSI using the selected denominator mode:

- Default mode: `site_depth_offset = max(donor_site_window_depth, acceptor_site_window_depth)`
- Strict PSI modes: `strict_local_depth = max(donor_decision_depth, acceptor_decision_depth)`

Current PSI-related outputs are:

- `workdir/psi.psi_values.tsv`
- `workdir/edgeR_results.intron_results_with_psi.tsv`
- `edgeR_results.all.tsv`

`delta_PSI` follows the requested contrast direction. For example, for `--contrast "perturb,control"`, `delta_PSI` is computed as:

`perturb_mean_PSI - control_mean_PSI`

## Plotting Utility

The plotting utility uses the PSI-enhanced results plus a GTF to generate a gene-centric PDF:

```bash
Rscript util/plot_intron_dexseq_like.R \
    --output_dir results/my_analysis \
    --gene RPL3 \
    --gtf /path/to/annotation.gtf \
    --contrast perturb_vs_control \
    --output results/my_analysis/RPL3_perturb_vs_control_intron_DEXseq_like.pdf
```

If `gene_name` is not present in the results, the plotting utility falls back to selecting introns by overlap with the gene locus in the GTF.
