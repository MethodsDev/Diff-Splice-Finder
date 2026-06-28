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

#### How to Interpret `logFC`

The reported `logFC` is estimated from the edgeR GLM after adding the shared
cluster-total offset:

```
log(expected intron count) = group effect + log(shared cluster total)
```

Because the offset is on the log scale, it is equivalent to a multiplier on the
count scale:

```
expected intron count = shared cluster total × exp(group effect)
```

You can think of the model as comparing:

```
log(intron count) - log(shared cluster total)
= log(intron count / shared cluster total)
```

So `logFC` is an **offset-adjusted intron usage log fold-change**, not a raw
intron-count fold-change and not a whole-gene expression fold-change. It asks how
the intron's proportional usage changes after accounting for the local splice
cluster abundance.

This is related to a PSI ratio, but not identical. A rough intuition is:

```
logFC ≈ log2(mean PSI in group A / mean PSI in group B)
```

However, the edgeR `logFC` is model-based: it uses counts, replicate structure,
dispersion estimates, and the GLM offset. A direct PSI ratio is a simple ratio of
summary PSI values and can be unstable when the denominator PSI is near zero.

`delta_PSI` is complementary:

```
delta_PSI = mean PSI in group A - mean PSI in group B
```

`delta_PSI` measures the absolute change in usage proportion, while `logFC`
measures the relative, model-estimated fold change in usage proportion. They are
two scales of the same shift, linked by `delta_PSI = PSI_control·(2^logFC − 1)` —
so you can convert between them only with the baseline (control) PSI, and neither
can be derived from the other alone. See
[docs/PSI_and_logFC.md](docs/PSI_and_logFC.md) for the derivation and worked examples.

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

For the gene-level DEXSeq-like plotting utility, also install:
```R
install.packages(c("dplyr", "tidyr", "readr", "stringr", "tibble", "ggplot2", "patchwork", "ggrepel"))
BiocManager::install("rtracklayer")
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

# Build count and site-depth offset matrices
python3 util/build_intron_count_matrix.py \
    --intron_files sample*.introns \
    --output_matrix intron_counts.matrix
    
# Optional: Compress matrix to save space (gzipped files are supported)
gzip intron_counts.matrix
```

When the `.introns` files contain `site_depth_offset`, the matrix builder also
writes `intron_counts.offsets.matrix`.

For complete site-depth offsets across all samples, run
`count_introns_from_bam.py` with a shared target intron list once that list is
known:

```bash
python3 util/count_introns_from_bam.py \
    --genome_fa reference.fa \
    --bam sample1.bam \
    --target_introns intron_counts.matrix > sample1.targeted.introns
```

**Note**: The pipeline automatically handles gzipped input files (`.gz` extension).

### 2. Create Sample Metadata
```tsv
sample_id	group
sample1	perturb
sample2	perturb
sample3	control
sample4	control
```

**Sample Subsetting:** The pipeline supports analyzing a subset of samples from your count matrix. Simply include only the samples you want to analyze in your metadata file. The pipeline will automatically filter the count and offset matrices to match. This is useful for:
- Comparing specific conditions from a larger dataset
- Excluding outlier samples
- Focusing on particular experimental contrasts

For example, if your count matrix has 50 samples but you only want to compare 10 of them, just list those 10 samples in your metadata file.

### 3. Run Differential Splicing Analysis
```bash
python3 run_diff_splice_analysis.py \
    --matrix intron_counts.matrix \
    --site_depth_offsets intron_counts.offsets.matrix \
    --samples sample_metadata.tsv \
    --output_dir results \
    --contrast "perturb,control"
```

That's it! Results are in `results/` directory.

### Alternative: Count BAMs During the Pipeline

Instead of providing count and offset matrices, `run_diff_splice_analysis.py`
can start from a BAM manifest:

```tsv
sample_type	replicate_id	bam_file
perturb	perturb_1	/path/to/perturb_1.bam
perturb	perturb_2	/path/to/perturb_2.bam
control	control_1	/path/to/control_1.bam
control	control_2	/path/to/control_2.bam
```

Then run:

```bash
python3 run_diff_splice_analysis.py \
    --samples samples.tsv \
    --genome_fa reference.fa \
    --output_dir results \
    --contrast "perturb,control"
```

In this mode the pipeline writes intermediate matrices under
`results/workdir/bam_inputs/` and uses site-depth offsets by default.

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

## Gene-Level Visualization

You can generate a gene-centric PDF view of intron results with transcript structure and per-group mean logCPM:

```bash
Rscript util/plot_intron_dexseq_like.R \
    --output_dir results \
    --gene RPL3 \
    --gtf annotation.gtf \
    --contrast perturb_vs_control \
    --output results/RPL3_perturb_vs_control_intron_DEXseq_like.pdf
```

Notes:
- The script expects the PSI-enhanced results file, typically `results/workdir/edgeR_results.intron_results_with_psi.tsv`.
- If your results contain multiple contrasts, pass `--contrast`. It accepts either `GroupA_vs_GroupB` or the pipeline contrast syntax `GroupA,GroupB1;GroupB2`.
- If `gene_name` is absent from the results, the script falls back to selecting introns by overlap with the target gene locus from the GTF.
- Add `--use_fitted_logcpm` to plot fitted rather than observed mean logCPM values.

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
   - Optional: use a precomputed splice-site depth matrix as the denominator for all introns, which supports singleton and retained-intron-like cases

4. **Filter low-confidence features**
   - Remove non-canonical splice sites (GT-AG, GC-AG, AT-AC only)
   - Filter introns with insufficient read support
   - Require threshold met for EITHER the donor or acceptor cluster when using shared offsets

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

8. **Add PSI and filter**
   - Merges PSI values with edgeR results
   - Filters by minimum |delta_PSI| (default: 0.05; set `--min_delta_psi 0` to disable)
   - Recalculates FDR on filtered set

### Running Individual Modules

For more control, run modules separately:

```bash
# 1. Cluster introns (both donor and acceptor)
python3 util/cluster_introns.py \
    --matrix intron_counts.matrix \
    --output_donor introns_clustered.tsv \
    --cluster_type both

# 2. Compute shared offsets only
python3 util/compute_offsets.py \
    --matrix introns_clustered.tsv \
    --output shared_offsets.tsv \
    --compute_offsets_only

# Optional: compute unstranded splice-site depth offsets from BAMs
# bams.tsv must contain columns: sample_id, bam
python3 util/compute_splice_site_depth_offsets.py \
    --matrix introns_clustered.tsv \
    --bam_list bams.tsv \
    --output site_depth_offsets.tsv \
    --window_radius 10

# Optional: compute offsets from a precomputed splice-site depth matrix
python3 util/compute_offsets.py \
    --matrix introns_clustered.tsv \
    --output shared_offsets.tsv \
    --compute_offsets_only \
    --site_depth_offsets site_depth_offsets.tsv \
    --offset_mode site_depth

# 3. Filter (requires donor OR acceptor support)
python3 util/filter_introns.py \
    --matrix introns_clustered.tsv \
    --output introns_filtered.tsv \
    --min_intron_count 10 \
    --min_cluster_count 20 \
    --min_cluster_samples 3

# 4. Prepare edgeR inputs
python3 util/compute_offsets.py \
    --matrix introns_filtered.tsv \
    --output_prefix edgeR_input \
    --shared_offsets shared_offsets.tsv

# 5. Run edgeR
Rscript util/run_edgeR_analysis.R \
    --counts edgeR_input.counts.tsv \
    --offsets edgeR_input.offsets.tsv \
    --annotations edgeR_input.annotations.tsv \
    --samples sample_metadata.tsv \
    --output results \
    --contrast "perturb,control"

# 6. Compute PSI from the prepared edgeR inputs
python3 util/compute_psi.py \
    --counts edgeR_input.counts.tsv \
    --annotations edgeR_input.annotations.tsv \
    --samples sample_metadata.tsv \
    --results results.intron_results.tsv \
    --output results_with_psi.tsv
```

The main pipeline can also compute site-depth offsets directly from BAMs:

```bash
python3 run_diff_splice_analysis.py \
    --matrix intron_counts.matrix \
    --samples sample_metadata.tsv \
    --output_dir results \
    --site_depth_bam_list bams.tsv
```

When either `--site_depth_bam_list` or `--site_depth_offsets` is provided and
`--offset_mode` is left at its default `auto`, the pipeline uses
`--offset_mode site_depth`.

## Output Files

### Key Results Files

**Intron-level results:**
- `edgeR_results.all.tsv` - All tested introns with edgeR statistics and PSI values
- `edgeR_results.significant_introns.tsv` - Final significant introns after the configured filters are applied. When `--min_delta_psi` is enabled, the pipeline first filters by `|delta_PSI|` only, recomputes `FDR` on that delta-PSI-filtered set, then reports rows passing the recomputed `FDR`, the configured `--fdr_threshold`, and `--min_logFC`.

**Intermediate files:**
- `workdir/introns_clustered.tsv` - Clustered matrix with donor_cluster and acceptor_cluster columns
- `workdir/introns_filtered.tsv` - After filtering low-confidence features
- `workdir/site_depth_offsets.tsv` - Site-depth denominator matrix, when computed from `--site_depth_bam_list`
- `workdir/shared_offsets.tsv` - Shared cluster totals used to build edgeR offsets and PSI denominators
- `workdir/shared_offsets.metadata.tsv` - Cluster sizes and offset source metadata
- `workdir/edgeR_input.counts.tsv` - Filtered count matrix used by edgeR
- `workdir/edgeR_input.offsets.tsv` - Log-transformed shared offsets used by edgeR
- `workdir/edgeR_input.annotations.tsv` - Intron annotations used by edgeR
- `workdir/edgeR_results.intron_results.tsv` - Raw edgeR output before PSI columns are added
- `workdir/psi.psi_values.tsv` - Per-sample PSI values, group means, and delta PSI
- `workdir/edgeR_results.intron_results_with_psi.psi_filtered.tsv` - Delta-PSI-filtered full result set used to generate the significant-only final file. In this file, `FDR_original` is the Benjamini-Hochberg FDR from the full tested set before delta-PSI filtering, while `FDR` is recalculated from `PValue` over only the rows passing `|delta_PSI| >= --min_delta_psi`. `logFC` is not used to choose the rows for this FDR recalculation.

**Diagnostics:**
- `workdir/edgeR_results.diagnostics.pdf` - BCV, dispersion, MA, volcano plots

### Result Interpretation

**Intron-level columns:**
- `intron_id`: Intron coordinates and splice sites (chr:start-end^splice_pair^flag)
- `donor_cluster`: Donor cluster ID (chr:donor_pos:strand)
- `acceptor_cluster`: Acceptor cluster ID (chr:acceptor_pos:strand)
- `donor_cluster_size`, `acceptor_cluster_size`: Number of introns observed in each splice-site cluster during offset calculation
- `both_splice_sites_singleton`: Whether the intron was singleton at both splice sites
- `offset_mode`: Denominator strategy, such as `site_depth` or `cluster_max`
- `offset_source`: `site_depth`, `cluster_max_donor`, `cluster_max_acceptor`, or `site_depth_fallback`
- `site_depth_fallback_used`: Whether site-depth was used only as a singleton fallback for this intron
- `gene_name`: Gene name (if GTF provided)
- `intron_status`: known/novel (if GTF provided)
- `logFC`: Model-based log2 fold-change in **offset-adjusted intron usage proportion** (not raw counts or gene expression)
  - Positive = increased usage in first group of contrast
  - Negative = decreased usage in first group of contrast
  - Roughly analogous to `log2(mean_PSI_group1 / mean_PSI_group2)`, but estimated by edgeR from counts with shared cluster-total offsets and dispersion modeling
- `logCPM`: Average log2 counts per million
- `F`: F-statistic from quasi-likelihood F-test
- `PValue`: P-value from quasi-likelihood F-test
- `FDR`: Benjamini-Hochberg adjusted p-value. In unfiltered result files, this is adjusted over all tested introns. In delta-PSI-filtered result files, this is recalculated over only the rows that pass the `--min_delta_psi` threshold; `logFC` is applied later when assigning final significance, not before FDR recalculation.
- `FDR_original`: Original Benjamini-Hochberg adjusted p-value from the full tested set before delta-PSI filtering; present in delta-PSI-filtered outputs for comparison with the recomputed `FDR`.
- `*_mean_PSI`: Mean PSI in each group (if PSI computed)
- `delta_PSI`: Difference in mean PSI between groups (if PSI computed)

**Important notes:**
- Each intron is tested **once** with shared offsets from both splice sites
- LogFC represents relative change in proportional usage within competing alternatives
- Delta PSI represents absolute change in proportional usage
- PSI uses the same shared denominators as edgeR for consistency
- Shared offsets prevent singleton cluster artifacts
- The final `edgeR_results.significant_introns.tsv` file contains only rows marked significant after the active filters are applied. With `--min_delta_psi` enabled, the FDR recalculation set is defined by `|delta_PSI| >= --min_delta_psi` alone; final reported rows then must pass the recomputed `FDR <= --fdr_threshold` and `|logFC| >= --min_logFC`.

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

## Testing

### Quick Integration Test
Run automated tests to verify the pipeline works correctly:
```bash
make test
```

This runs a quick test (~1-2 minutes) using a small dataset (550 introns). It validates:
- Clustering with both donor and acceptor
- Shared offset computation
- Filtering with dual cluster thresholds
- edgeR analysis
- PSI calculation with shared denominators
- All expected output files are created

### Visualization Test
Run the DEXSeq-like PDF plotting smoke test:
```bash
make test-viz
```

This exercises the plotting CLI against committed quick-test results plus a synthetic GTF fixture and checks that a non-empty PDF is written.

### Full Integration Test
Test all features including gene annotation and PSI filtering:
```bash
make test-full
```

### Testing From `testing/`
```bash
make -C testing test
make -C testing test_viz
make -C testing clean
```

The shell scripts in `testing/` are used by the local Makefile and are not the recommended public entrypoint.

See [testing/README.md](testing/README.md) for fixture details and expected outputs.

## Design Principles

This pipeline implements specific design choices detailed in [docs/AI_ONBOARDING.md](docs/AI_ONBOARDING.md) and [docs/SHARED_OFFSETS.md](docs/SHARED_OFFSETS.md):

1. **Splicing is compositional** - Test relative usage within clusters using shared offsets
2. **Shared offsets** - Uses max(donor_total, acceptor_total) to prevent singleton artifacts
3. **Intron-level analysis** - Each intron tested once with comprehensive information
4. **Counts remain counts** - No PSI/TPM transformation; normalization via GLM offsets
5. **Technology-agnostic** - Same framework for short and long reads
6. **Annotation-light** - Introns defined from observed junctions
7. **Robust statistics** - edgeR's QL framework with robust dispersion

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
