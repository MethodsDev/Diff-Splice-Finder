# Diff-Splice-Finder

**Intron-focused differential splicing analysis for bulk RNA-seq**

A robust pipeline for detecting differential intron usage in both short-read (Illumina) and long-read (PacBio Kinnex) RNA-seq data using edgeR with compositional normalization.

## Overview

Diff-Splice-Finder identifies changes in splicing by testing **intron usage proportions within local splicing clusters**, rather than absolute intron counts. This approach:

- ✅ Works consistently across short and long-read technologies
- ✅ Separates splicing changes from expression changes
- ✅ Uses edgeR's robust statistical framework with custom offsets
- ✅ Provides both intron-level and cluster-level results
- ✅ Requires minimal annotation (splice junctions define features)

## Key Concepts

### Compositional Analysis
Splicing is inherently compositional - it reflects **choices among introns** within a locus. We model this using cluster-total offsets in edgeR GLMs:

```
log(μ_i,s) = X_s × β_i + log(T_C,s)
```

where `T_C,s` is the total intron support for cluster C in sample s.

### Clustering Strategies
- **Donor clusters**: Introns sharing the same 5' splice site (captures alternative acceptors)
- **Acceptor clusters**: Introns sharing the same 3' splice site (captures alternative donors)

Both are run by default to comprehensively detect splicing changes.

## Installation

### Requirements

**Python 3.7+** with packages:
```bash
pip install pandas numpy scipy statsmodels pysam
```

**R 4.0+** with Bioconductor:
```R
install.packages("BiocManager")
BiocManager::install("edgeR")
install.packages("optparse")
```

### Clone Repository
```bash
git clone https://github.com/MethodsDev/Diff-Splice-Finder.git
cd Diff-Splice-Finder
chmod +x run_diff_splice_analysis.py
chmod +x util/*.py
chmod +x util/*.R
chmod +x examples/*.sh
```

## Quick Start

### 1. Count Introns from BAM Files
```bash
# For each sample
python3 util/count_introns_from_bam.py \
    --genome_fa reference.fa \
    --bam sample1.bam > sample1.introns

# Build count matrix
python3 util/build_intron_count_matrix.py \
    --intron_files sample*.introns \
    --output_matrix intron_counts.matrix
```

### 2. Create Sample Metadata
```tsv
sample_id	group
sample1	TDP43
sample2	TDP43
sample3	control
sample4	control
```

### 3. Run Differential Splicing Analysis
```bash
python3 run_diff_splice_analysis.py \
    --matrix intron_counts.matrix \
    --samples sample_metadata.tsv \
    --output_dir results \
    --contrast "TDP43-control"
```

That's it! Results are in `results/donor/` and `results/acceptor/`.

### Resume on Crash

The pipeline automatically resumes where it left off if interrupted:

```bash
# Run pipeline
python3 run_diff_splice_analysis.py ...

# If it crashes or is interrupted, just rerun the same command
python3 run_diff_splice_analysis.py ...
# Will skip completed steps and resume

# To force complete rerun from scratch:
python3 run_diff_splice_analysis.py ... --force_rerun
```

The pipeline checks for existing output files and skips completed steps, saving time and avoiding redundant computation.

## Detailed Workflow

### Step-by-Step Pipeline

The main script orchestrates these steps:

1. **Cluster introns** by shared donor/acceptor sites
   - Groups introns into local splicing clusters
   - Both donor and acceptor clustering by default

2. **Filter low-confidence features**
   - Remove non-canonical splice sites (GT-AG, GC-AG, AT-AC only)
   - Filter introns with insufficient read support
   - Filter clusters with too few reads

3. **Compute cluster-total offsets**
   - Calculate total intron support per cluster per sample
   - Log-transform for use as GLM offsets

4. **Run edgeR analysis**
   - Negative binomial GLM with QL framework
   - NO library size normalization (norm.factors = 1)
   - All normalization via cluster-total offsets
   - Tests intron usage proportions, not expression

5. **Aggregate to cluster level**
   - Combine intron p-values per cluster (Cauchy method)
   - Report cluster-level FDR and summary statistics

### Running Individual Modules

For more control, run modules separately:

```bash
# 1. Cluster introns
python3 util/cluster_introns.py \
    --matrix intron_counts.matrix \
    --output_donor donor_clustered.tsv \
    --cluster_type donor

# 2. Filter
python3 util/filter_introns.py \
    --matrix donor_clustered.tsv \
    --output donor_filtered.tsv \
    --cluster_type donor \
    --min_intron_count 10 \
    --min_cluster_count 20

# 3. Compute offsets
python3 util/compute_offsets.py \
    --matrix donor_filtered.tsv \
    --output_prefix donor_edgeR \
    --cluster_type donor

# 4. Run edgeR
Rscript util/run_edgeR_analysis.R \
    --counts donor_edgeR.counts.tsv \
    --offsets donor_edgeR.offsets.tsv \
    --annotations donor_edgeR.annotations.tsv \
    --samples sample_metadata.tsv \
    --output donor_results \
    --contrast "TDP43-control"

# 5. Aggregate clusters
python3 util/aggregate_clusters.py \
    --intron_results donor_results.intron_results.tsv \
    --output_prefix donor_aggregated \
    --method cauchy
```

## Output Files

### Key Results Files

**Integrated results (combining donor + acceptor):**
- `integrated.integrated_results.tsv` - All introns from both analyses with comparative statistics
- `integrated.significant_integrated.tsv` - Significant introns only (in either or both analyses)
- `integrated.integration_summary.tsv` - High-level summary statistics

**Intron-level results (per cluster type):**
- `{cluster_type}_edgeR_results.intron_results.tsv` - All tested introns with statistics
- `{cluster_type}_edgeR_results.significant_introns.tsv` - Significant introns only

**Cluster-level results (per cluster type):**
- `{cluster_type}_aggregated.cluster_results.tsv` - Cluster-level significance
- `{cluster_type}_aggregated.significant_cluster_introns.tsv` - Introns in significant clusters

**Diagnostics:**
- `{cluster_type}_edgeR_results.diagnostics.pdf` - BCV, dispersion, MA, volcano plots

### Result Interpretation

**Integrated results columns:**
- `intron_id`: Intron coordinates and splice sites
- `tested_in`: Where the intron was tested (both/donor_only/acceptor_only)
- `significant_in`: Where the intron is significant (both/donor_only/acceptor_only/neither)
- `best_analysis`: Which clustering gave the most significant result (donor/acceptor)
- `best_FDR`: Minimum FDR across both analyses
- `best_logFC`: LogFC from the best (most significant) analysis
- `direction_consistent`: For introns significant in both - do they agree on direction?
- `donor_*`: Results from donor clustering
- `acceptor_*`: Results from acceptor clustering

**Intron-level columns (per cluster type):**
- `logFC`: Log2 fold-change in **intron usage proportion** (NOT expression)
  - Positive = increased usage in first group
  - Negative = decreased usage in first group
- `PValue`: P-value from quasi-likelihood F-test
- `FDR`: Benjamini-Hochberg adjusted p-value

**Cluster-level columns (per cluster type):**
- `cluster_pvalue`: Combined p-value across introns in cluster
- `cluster_FDR`: Cluster-level FDR correction
- `dominant_direction`: Overall splicing pattern
  - **increased**: Majority of significant introns show increased usage
  - **decreased**: Majority of significant introns show decreased usage
  - **mixed**: Equal counts of increased/decreased (complex alternative splicing)
  - **none**: No significant introns
- `n_increased`: Count of introns with increased usage
- `n_decreased`: Count of introns with decreased usage
- `n_significant_introns`: Number of introns with FDR < 0.05

## Parameter Tuning

### Default Parameters (balanced)
```bash
--min_intron_count 10 \
--min_intron_samples 2 \
--min_cluster_count 20 \
--min_cluster_samples 3
```

### For Long-Read Data (more lenient)
```bash
--min_intron_count 5 \
--min_intron_samples 2 \
--min_cluster_count 10 \
--min_cluster_samples 2
```

### For Conservative Analysis (fewer false positives)
```bash
--min_intron_count 20 \
--min_cluster_count 50 \
--min_logFC 1.0 \
--fdr_threshold 0.01
```

See [examples/PARAMETER_GUIDE.md](examples/PARAMETER_GUIDE.md) for detailed guidance.

## Example Analyses

### Basic Analysis
```bash
./examples/run_example_analysis.sh
```

### With Batch Correction
```bash
./examples/run_with_batch_correction.sh
```

### Custom Parameters
```bash
python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/custom \
    --min_intron_count 20 \
    --min_cluster_count 30 \
    --min_logFC 0.5 \
    --fdr_threshold 0.05 \
    --batch_col batch
```

## Design Principles

This pipeline implements specific design choices detailed in [AI_ONBOARDING.md](AI_ONBOARDING.md):

1. **Splicing is compositional** - Test relative usage within clusters, not absolute counts
2. **Counts remain counts** - No PSI/TPM transformation; normalization via GLM offsets
3. **Technology-agnostic** - Same framework for short and long reads
4. **Annotation-light** - Introns defined from observed junctions
5. **Robust statistics** - edgeR's QL framework with robust dispersion

## Non-Goals

This pipeline does NOT attempt to:
- Assemble or quantify transcript isoforms
- Detect isoform switching
- Calculate PSI/Ψ from exon inclusion
- Handle single-cell data
- Perform UMI deduplication

For isoform-level analysis, consider IsoformSwitchAnalyzeR, SUPPA2, or rMATS.

## Citation

If you use this pipeline, please cite:
- **edgeR**: Robinson MD, McCarthy DJ, Smyth GK (2010). Bioinformatics 26(1):139-140
- **ACAT** (if using Cauchy aggregation): Liu Y, Xie J (2020). JASA 115(529):393-402

## Support

For issues, questions, or feature requests:
- GitHub Issues: https://github.com/MethodsDev/Diff-Splice-Finder/issues
- See [AI_ONBOARDING.md](AI_ONBOARDING.md) for design rationale
- See [examples/PARAMETER_GUIDE.md](examples/PARAMETER_GUIDE.md) for parameter tuning

## License

MIT License - See LICENSE file for details
