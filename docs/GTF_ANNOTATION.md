# GTF Annotation and PSI Calculation Features

## Overview

The differential splicing analysis pipeline now supports optional GTF annotation to enrich intron results with gene information and known/novel intron status.

## Usage

### Main Pipeline

Add the `--gtf` parameter to include gene annotations in the integrated results:

```bash
python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/my_analysis \
    --gtf /path/to/annotation.gtf
```

### Integration Script Only

You can also add annotations when running the integration script directly:

```bash
python3 util/integrate_results.py \
    --donor_results results/donor/donor_edgeR_results.intron_results.tsv \
    --acceptor_results results/acceptor/acceptor_edgeR_results.intron_results.tsv \
    --output_prefix results/integrated \
    --gtf /path/to/annotation.gtf
```

## Output Columns

When a GTF file is provided, three additional columns are added to the integrated results:

1. **gene_name**: The best-matching gene name for the intron (based on genomic overlap)
2. **intron_status**: Classification of the intron:
   - `known`: Exact match to an annotated intron in the GTF
   - `novel`: Intron not found in the annotation (but may still be within a known gene)
   - `unknown`: Failed to parse intron coordinates
3. **overlapping_genes**: Comma-separated list of all genes overlapping the intron

## Intron ID Format Support

The annotation system supports multiple intron ID formats:

1. **Standard format**: `chr:start-end:strand`
   - Example: `chr1:1000-2000:+`

2. **Splice motif format**: `chr:start-end^DONOR--ACCEPTOR^STATUS`
   - Example: `chr1:1000-2000^GT--AG^OK`
   - Strand is automatically inferred from splice motif:
     - `GT--AG` → plus strand (+)
     - `CT--AC` → minus strand (-)

## Known vs Novel Introns

**Known introns** are identified by exact coordinate match:
- The intron must have the same chromosome, start position, end position, and strand as an intron derived from the GTF annotation
- Introns are extracted from GTF by identifying gaps between consecutive exons within transcripts

**Novel introns** are either:
- Introns within annotated genes but with different boundaries than known isoforms
- Introns in unannotated regions of the genome

## Gene Assignment

When multiple genes overlap an intron:
- The `gene_name` column contains the first (best) match
- The `overlapping_genes` column contains all overlapping genes (comma-separated)

Gene overlap is determined by:
1. Direct overlap with gene features in the GTF
2. Overlap with transcript spans (from first to last exon)

## GTF Format Requirements

The GTF file should follow standard GTF2.2 format with:
- `gene` features containing `gene_id` and `gene_name` attributes
- `exon` features with `gene_id`, `transcript_id`, and `gene_name` attributes

Example GTF entries:
```
chr1  HAVANA  gene        1000    5000    .  +  .  gene_id "ENSG001"; gene_name "GENE1";
chr1  HAVANA  transcript  1000    5000    .  +  .  gene_id "ENSG001"; transcript_id "ENST001"; gene_name "GENE1";
chr1  HAVANA  exon        1000    2000    .  +  .  gene_id "ENSG001"; transcript_id "ENST001"; gene_name "GENE1";
chr1  HAVANA  exon        3000    5000    .  +  .  gene_id "ENSG001"; transcript_id "ENST001"; gene_name "GENE1";
```

## Backward Compatibility

The `--gtf` parameter is optional. If not provided:
- The pipeline runs normally without gene annotations
- Output files maintain their original structure
- No breaking changes to existing workflows

## Performance Notes

- GTF parsing happens once during the integration step
- **Parsed data is cached** in a `.intron_cache.tsv` file alongside the GTF
- Subsequent runs automatically use the cache for instant loading
- Cache is automatically invalidated if the GTF file is modified
- For human genome annotations (e.g., GENCODE):
  - First run: ~1-2 minutes parsing time
  - Subsequent runs: < 5 seconds to load from cache
- Memory usage scales with GTF file size
- Cache files are human-readable TSV format

### Cache Files

When you provide a GTF file, a cache file is automatically created:
- **Location**: Same directory as GTF, named `<gtf_basename>.intron_cache.tsv`
- **Format**: Tab-separated values with two sections:
  1. Annotated introns with their gene assignments
  2. Gene regions for overlap checking
- **Automatic invalidation**: If GTF is newer than cache, cache is regenerated
- **Manual deletion**: Simply delete the `.intron_cache.tsv` file to force re-parsing

Example:
```bash
# Using GTF: gencode.v45.annotation.gtf
# Cache created: gencode.v45.annotation.intron_cache.tsv
# First run parses GTF and creates cache (~2 min for human genome)
# Second run loads from cache (~5 seconds)
```

---

## PSI (Percent Spliced In) Calculations

### Overview

In addition to statistical testing, you can compute PSI values to quantify the relative usage of each intron within its cluster. PSI provides an intuitive metric similar to those reported by tools like LeafCutter.

**PSI Definition:**
```
PSI_i,s = count_i,s / Σ(count_j,s for all j in cluster_i)
```

Where:
- `PSI_i,s` = Percent Spliced In for intron i in sample s
- `count_i,s` = Read count for intron i in sample s  
- `cluster_i` = All introns sharing a splice site with intron i

### Computing PSI

Use the `compute_psi.py` utility to add PSI values to your results:

```bash
python3 util/compute_psi.py \
    --counts results/donor/donor_edgeR_input.counts.tsv \
    --annotations results/donor/donor_edgeR_input.annotations.tsv \
    --samples examples/sample_metadata.tsv \
    --results results/donor/donor_edgeR_results.intron_results.tsv \
    --output results/donor/donor_with_PSI.tsv \
    --cluster_col donor_cluster
```

### Output Columns

The PSI-enhanced results include:

**Per-Group Statistics:**
- `{group}_mean_PSI`: Mean PSI across samples in the group
- `{group}_median_PSI`: Median PSI across samples  
- `{group}_std_PSI`: Standard deviation of PSI

**Differential Metrics (for 2-group comparisons):**
- `delta_PSI`: Difference in mean PSI between groups (group1 - group2)
- `abs_delta_PSI`: Absolute value of delta PSI

**Optional Per-Sample Values:**
Add `--keep_individual_samples` to include PSI for each individual sample.

### Interpretation

**PSI Values:**
- PSI = 0.75 means intron accounts for 75% of splicing at that site
- PSI = 0.25 means intron is a minor isoform (25% usage)

**Delta PSI:**
- |ΔPSI| > 0.1 suggests biologically relevant change
- |ΔPSI| > 0.2 indicates strong differential usage
- Sign indicates direction: positive = increased in group1

### Example Results

```
intron_id                          logFC    FDR      TDP43_PSI  control_PSI  delta_PSI
chr7:100650-100700^GT--AG^OK      -1.23    0.007    0.42       0.99         -0.57
chr22:42599794-42602769^CT--AC^OK  2.16    0.008    0.76       0.17          0.59
```

In the first example:
- Intron accounts for 42% of splicing in TDP43 samples
- Same intron accounts for 99% in controls
- ΔPSI = -0.57 indicates dramatically reduced usage in TDP43

### Combining with Statistical Results

PSI complements statistical testing:
- **logFC**: Tests if proportion differs (accounts for uncertainty)
- **PSI**: Quantifies the actual proportions (easier to interpret)
- **delta_PSI**: Direct measure of biological effect size

**Best practice:** Filter for both:
- Statistical significance (FDR < 0.05)
- Biological significance (|ΔPSI| > 0.1 or 0.2)

### Integration with Pipeline

You can chain PSI calculation into your workflow:

```bash
# Run main pipeline
python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples metadata.tsv \
    --output_dir results

# Add PSI to donor results
python3 util/compute_psi.py \
    --counts results/donor/donor_edgeR_input.counts.tsv \
    --annotations results/donor/donor_edgeR_input.annotations.tsv \
    --samples metadata.tsv \
    --results results/donor/donor_edgeR_results.intron_results.tsv \
    --output results/donor/donor_with_PSI.tsv \
    --cluster_col donor_cluster

# Add PSI to acceptor results  
python3 util/compute_psi.py \
    --counts results/acceptor/acceptor_edgeR_input.counts.tsv \
    --annotations results/acceptor/acceptor_edgeR_input.annotations.tsv \
    --samples metadata.tsv \
    --results results/acceptor/acceptor_edgeR_results.intron_results.tsv \
    --output results/acceptor/acceptor_with_PSI.tsv \
    --cluster_col acceptor_cluster
```
