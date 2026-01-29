# Diff-Splice-Finder

**Intron-focused differential splicing analysis for bulk RNA-seq**

A robust pipeline for detecting differential intron usage in both short-read (Illumina) and long-read (PacBio Kinnex) RNA-seq data using edgeR with compositional normalization.

## Overview

Diff-Splice-Finder identifies changes in splicing by testing **intron usage proportions within local splicing clusters**, rather than absolute intron counts. This approach:

- ✅ Works consistently across short and long-read technologies
- ✅ Separates splicing changes from expression changes
- ✅ Uses edgeR's robust statistical framework with shared offsets
- ✅ Tests each intron once with comprehensive information from both splice sites
- ✅ Requires minimal annotation (splice junctions define features)

## Key Concepts

### Compositional Analysis
Splicing is inherently compositional - it reflects **choices among introns** within a locus. We model this using shared cluster-total offsets in edgeR GLMs:

```
log(μ_i,s) = X_s × β_i + log(T_shared,s)
```

where `T_shared,s = max(T_donor,s, T_acceptor,s)` is the maximum of the donor and acceptor cluster totals for intron i in sample s. This ensures consistent normalization and prevents singleton cluster artifacts.

### How the Statistical Testing Works

#### The Core Problem

When a gene is alternatively spliced, different introns compete for usage - if one intron is used more, others in the same cluster must be used less (they sum to the total splicing events at that site). This is called a **compositional constraint**.

#### Why We Need Offsets

Consider a donor site where 3 introns compete:
- **Intron A**: 80 reads in condition 1, 40 reads in condition 2  
- **Intron B**: 15 reads in condition 1, 10 reads in condition 2
- **Intron C**: 5 reads in condition 1, 50 reads in condition 2
- **Cluster total**: 100 reads in condition 1, 100 reads in condition 2

Looking at raw counts, Intron A decreased by 50%. But the cluster total stayed at 100 - so this isn't about expression changes, it's about **switching** between introns.

The **offset** represents the log-transformed cluster total. By incorporating this as an offset in the statistical model, we tell the algorithm: "Don't treat these as independent events - they're parts of a whole."

#### How edgeR Works with Offsets

edgeR uses a **negative binomial generalized linear model** (GLM):

```
log(μ) = β₀ + β₁×Group + log(ClusterTotal)
         ↑              ↑
    baseline    group effect    ← offset (fixed)
```

**Key points:**

1. **The offset is fixed** - it's not estimated, it's given. This forces the model to compare **proportions within clusters** rather than raw counts.

2. **The test asks**: "Is the proportion of this intron (relative to cluster total) different between groups?"

3. **Log fold-change interpretation**: 
   - logFC = 2 means the intron is used **4× more** in group A vs B (as a proportion of the cluster)
   - logFC = -1 means it's used **50% less** (2× less)

#### The Statistical Test

For each intron, edgeR:
1. **Estimates dispersion** (biological variability between replicates)
2. **Fits the model** accounting for the offset
3. **Tests the null hypothesis**: "Group coefficient β₁ = 0" (no difference in proportional usage)
4. **Computes p-values** using a quasi-likelihood F-test
5. **Adjusts for multiple testing** → FDR (False Discovery Rate)

#### What Makes an Intron "Significant"?

You need **both**:
- **FDR < 0.05** (statistically reliable, accounting for testing thousands of introns)
- **|logFC| ≥ threshold** (biologically meaningful effect size)

An intron with FDR=0.001 but logFC=0.1 might be statistically significant but biologically uninteresting (only 7% change in proportion).

#### Why This Approach Works

Traditional differential expression tools (like analyzing total gene counts) would fail here because:
- They'd detect changes even when **total splicing stays the same** but switching occurs
- They'd miss changes when **expression increases** but one intron's proportion drops

The offset-based approach **isolates the splicing changes** from expression changes, giving you a clean answer to: "Did the splicing pattern change between conditions?"

### Clustering and Shared Offsets
- **Donor clusters**: Introns sharing the same 5' splice site (captures alternative acceptors)
- **Acceptor clusters**: Introns sharing the same 3' splice site (captures alternative donors)
- **Shared offsets**: For each intron, uses `max(donor_total, acceptor_total)` as the denominator

This approach:
- Prevents singleton cluster artifacts (e.g., novel acceptor paired with common donor)
- Ensures consistent offsets across all analyses
- Each intron tested once with comprehensive information from both splice sites

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
    
# Optional: Compress matrix to save space (gzipped files are supported)
gzip intron_counts.matrix
```

**Note**: The pipeline automatically handles gzipped input files (`.gz` extension).

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

That's it! Results are in `results/` directory.

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

1. **Cluster introns** by shared donor AND acceptor sites
   - Groups introns into local splicing clusters
   - Creates both donor_cluster and acceptor_cluster columns

2. **Annotate with genes** (if GTF provided)
   - Maps introns to gene names
   - Marks known vs novel intron status

3. **Compute shared offsets**
   - For each intron, calculates donor cluster total and acceptor cluster total
   - Uses `max(donor_total, acceptor_total)` as shared offset
   - Prevents singleton cluster artifacts

4. **Filter low-confidence features**
   - Remove non-canonical splice sites (GT-AG, GC-AG, AT-AC only)
   - Filter introns with insufficient read support
   - Require threshold met for BOTH donor and acceptor clusters

5. **Prepare edgeR inputs**
   - Creates count, offset, and annotation files
   - Log-transform shared offsets for GLM

6. **Run edgeR analysis**
   - Negative binomial GLM with QL framework
   - NO library size normalization (norm.factors = 1)
   - All normalization via shared cluster-total offsets
   - Tests intron usage proportions, not expression
   - Each intron tested once

7. **Compute PSI values**
   - PSI = intron_count / shared_cluster_total
   - Uses same denominators as edgeR for consistency
   - Calculates group means and delta PSI

8. **Add PSI and filter** (optional)
   - Merges PSI values with edgeR results
   - Optionally filters by minimum |delta_PSI|
   - Recalculates FDR on filtered set

### Running Individual Modules

For more control, run modules separately:

```bash
# 1. Cluster introns (both donor and acceptor)
python3 util/cluster_introns.py \
    --matrix intron_counts.matrix \
    --output_donor introns_clustered.tsv \
    --cluster_type both

# 2. Compute shared offsets
python3 util/compute_offsets.py \
    --matrix introns_clustered.tsv \
    --output_prefix shared_offsets \
    --shared_offsets

# 3. Filter (requires both cluster thresholds)
python3 util/filter_introns.py \
    --matrix introns_clustered.tsv \
    --output introns_filtered.tsv \
    --cluster_type both \
    --min_intron_count 10 \
    --min_cluster_count 20

# 4. Prepare edgeR inputs
python3 util/compute_offsets.py \
    --matrix introns_filtered.tsv \
    --output_prefix edgeR_input \
    --shared_offsets_file shared_offsets.raw_cluster_totals.tsv

# 5. Run edgeR
Rscript util/run_edgeR_analysis.R \
    --counts edgeR_input.counts.tsv \
    --offsets edgeR_input.offsets.tsv \
    --annotations edgeR_input.annotations.tsv \
    --samples sample_metadata.tsv \
    --output results \
    --contrast "TDP43-control"

# 6. Compute PSI
python3 util/compute_psi.py \
    --counts edgeR_input.counts.tsv \
    --annotations edgeR_input.annotations.tsv \
    --samples sample_metadata.tsv \
    --shared_cluster_totals shared_offsets.raw_cluster_totals.tsv \
    --output psi_values.tsv
```

## Output Files

### Key Results Files

**Intron-level results:**
- `edgeR_results.intron_results.tsv` - All tested introns with statistics from edgeR
- `edgeR_results.intron_results_with_psi.tsv` - Results with PSI values added (unfiltered)
- `edgeR_results.intron_results_with_psi.psi_filtered.tsv` - Filtered by delta PSI threshold (if specified)
- `edgeR_results.significant_introns.tsv` - Significant introns only (FDR < threshold)

**PSI values:**
- `psi.psi_values.tsv` - Per-sample PSI values, group means, and delta PSI

**Shared offsets:**
- `shared_offsets.raw_cluster_totals.tsv` - Raw cluster totals (max of donor/acceptor)
- `shared_offsets.log_offsets.tsv` - Log-transformed offsets used by edgeR

**Intermediate files:**
- `introns_clustered.tsv` - Clustered matrix with donor_cluster and acceptor_cluster columns
- `introns_filtered.tsv` - After filtering low-confidence features

**Diagnostics:**
- `edgeR_results.diagnostics.pdf` - BCV, dispersion, MA, volcano plots

### Result Interpretation

**Intron-level columns:**
- `intron_id`: Intron coordinates and splice sites (chr:start-end^splice_pair^flag)
- `donor_cluster`: Donor cluster ID (chr:donor_pos:strand)
- `acceptor_cluster`: Acceptor cluster ID (chr:acceptor_pos:strand)
- `gene_name`: Gene name (if GTF provided)
- `intron_status`: known/novel (if GTF provided)
- `logFC`: Log2 fold-change in **intron usage proportion** (NOT expression)
  - Positive = increased usage in first group of contrast
  - Negative = decreased usage in first group of contrast
- `logCPM`: Average log2 counts per million
- `F`: F-statistic from quasi-likelihood F-test
- `PValue`: P-value from quasi-likelihood F-test
- `FDR`: Benjamini-Hochberg adjusted p-value
- `*_mean_PSI`: Mean PSI in each group (if PSI computed)
- `delta_PSI`: Difference in mean PSI between groups (if PSI computed)

**Important notes:**
- Each intron is tested **once** with shared offsets from both splice sites
- LogFC represents change in proportional usage within competing alternatives
- PSI uses the same shared denominators as edgeR for consistency
- Shared offsets prevent singleton cluster artifacts

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

### With Control Groups
When you have specific control samples and want to compare all treatment groups against them:
```bash
./examples/run_with_control_groups.sh
```

Or specify directly:
```bash
python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/with_controls \
    --control_groups control
```

For multiple control types (e.g., control and wildtype), use comma-separated values:
```bash
python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata_with_controls.tsv \
    --output_dir results/with_controls \
    --control_groups control,wildtype
```

**Benefits of control-based comparisons:**
- Focuses on biologically relevant comparisons (treatment vs control)
- Reduces multiple testing burden (fewer comparisons = better FDR)
- More interpretable results for experimental designs with clear control groups

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
