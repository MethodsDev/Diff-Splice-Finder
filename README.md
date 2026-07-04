# Diff-Splice-Finder

**Intron-focused differential splicing analysis for bulk RNA-seq**

A robust pipeline for detecting differential intron usage in bulk RNA-seq data using edgeR with selected local depth offsets.

## Overview

Diff-Splice-Finder identifies changes in splicing by testing **individual intron usage relative to local splice-site read depth**, rather than absolute intron counts. This approach:

- ✅ Works consistently across short and long-read technologies
- ✅ Separates splicing changes from expression changes
- ✅ Uses edgeR's robust statistical framework with local depth offsets
- ✅ Tests each intron once without donor/acceptor clustering
- ✅ Requires minimal annotation (splice junctions define features)

## Key Concepts

### Depth Offset Model

For each intron and sample, DSF uses a selected local read-depth denominator.
The default denominator is:

```
max_splice_plus_retained_depth =
    max(left_splice_depth + left_retained_depth,
        right_splice_depth + right_retained_depth)
```

Splice depths are sums of canonical intron counts sharing each splice-site
boundary. Retained depths are intron-interior coverage windows computed with
`samtools depth`, and by default begin 20 bases inside the intron.
An optional gene-scoped mode (`splice_vs_rest`) extends this with gene-total
junction evidence; see [Offset Modes](#offset-modes) below.
The edgeR GLM offset mode uses the selected denominator as a fixed log-offset:

```
log(μ_i,s) = X_s × β_i + log(D_i,s)
```

where `D_i,s` is the selected depth denominator for intron `i` in sample `s`.

### How the Statistical Testing Works

#### The Core Problem

Raw intron junction counts mix splicing usage with local expression and coverage. A junction can have more reads simply because the locus is more deeply sequenced in one sample. The depth offset asks whether the intron count changes **relative to local splice-site coverage**.

#### Why We Need Offsets

Consider one intron:
- **Condition 1**: 20 junction reads, 100 reads of selected local depth
- **Condition 2**: 40 junction reads, 400 reads of selected local depth

The raw junction count doubled, but usage relative to local depth fell from 20% to 10%. The offset lets edgeR model that proportional usage change instead of treating the raw count increase as higher intron usage.

The **offset** is the log-transformed selected denominator. It is supplied to the model, not estimated from the group coefficient.

#### How edgeR Works with Offsets

edgeR uses a **negative binomial generalized linear model** (GLM):

```
log(μ) = β₀ + β₁×Group + log(Denominator)
         ↑              ↑
    baseline    group effect    ← offset (fixed)
```

**Key points:**

1. **The offset is fixed** - it's not estimated, it's given. This makes the model compare intron counts relative to selected local depth rather than raw counts.

2. **The test asks**: "Is this intron's usage relative to local splice-site coverage different between groups?"

3. **Log fold-change interpretation**: 
   - logFC = 2 means the intron is used **4× more** in group A vs B after accounting for the selected denominator
   - logFC = -1 means it's used **50% less** (2× less)

#### How to Interpret `logFC`

The reported `logFC` is estimated from the edgeR GLM after adding the selected depth
offset:

```
log(expected intron count) = group effect + log(selected denominator)
```

Because the offset is on the log scale, it is equivalent to a multiplier on the
count scale:

```
expected intron count = selected denominator × exp(group effect)
```

You can think of the model as comparing:

```
log(intron count) - log(selected denominator)
= log(intron count / selected denominator)
```

So `logFC` is an **offset-adjusted intron usage log fold-change**, not a raw
intron-count fold-change and not a whole-gene expression fold-change. It asks how
the intron's usage changes after accounting for local splice-site coverage.

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
5. **Adjusts for multiple testing over the prefiltered test set** → FDR (False Discovery Rate)

#### What Makes an Intron "Significant"?

You need **both**:
- **FDR < 0.05** (statistically reliable, accounting for testing thousands of introns)
- **|logFC| ≥ threshold** (biologically meaningful effect size)
- The intron must have passed the pre-test `--min_delta_psi` filter if enabled

An intron with FDR=0.001 but logFC=0.1 might be statistically significant but biologically uninteresting (only 7% change in proportion).

#### Why This Approach Works

Traditional differential expression tools (like analyzing total gene counts) would fail here because:
- They'd detect changes even when **total splicing stays the same** but switching occurs
- They'd miss changes when **expression increases** but one intron's proportion drops

The offset-based approach **isolates intron usage changes** from coverage changes, giving you a clean answer to: "Did this intron's usage relative to local splice-site depth change between conditions?"

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

# Build count and depth denominator matrices
python3 util/build_intron_count_matrix.py \
    --intron_files sample*.introns \
    --output_matrix intron_counts.matrix
    
# Optional: Compress matrix to save space (gzipped files are supported)
gzip intron_counts.matrix
```

When the `.introns` files contain depth columns, the matrix builder also writes
column-specific denominator matrices such as
`intron_counts.max_adjacent_depth.matrix` and
`intron_counts.max_splice_plus_retained_depth.matrix`.
The `splice_plus_retained` analysis mode uses
`intron_counts.max_splice_plus_retained_depth.matrix`.

For complete depth denominators across all samples, run
`count_introns_from_bam.py` with a shared target intron list once that list is
known:

```bash
python3 util/count_introns_from_bam.py \
    --genome_fa reference.fa \
    --bam sample1.bam \
    --target_introns intron_counts.matrix > sample1.targeted.introns
```

For strand-specific depth denominators, add
`--site_depth_strand_mode F`, `R`, `FR`, or `RF`. The default is
`unstranded`.

Retained-depth windows skip the first `--retained_depth_inner_offset` intronic
bases at each splice boundary before measuring intron-body depth. The default
inner offset is `20`; the retained window width defaults to
`--site_depth_window_radius` unless `--retained_depth_window_radius` is set.

Stranded mode applies to the **depth denominator**. Junction counts are
discovered from split alignments without an explicit read-orientation filter;
for canonical splice motifs, the intron strand is inferred from the reference
dinucleotides. This is usually sufficient for the junction numerator, while the
stranded depth mode prevents antisense/local-coverage contamination in the
denominator. Use `unstranded` for unstranded libraries.

When the intron-counting WDL is run with a stranded mode, it also returns
`*.transcript_plus.bam(.bai)` and `*.transcript_minus.bam(.bai)` files.

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
    --offset_matrix intron_counts.max_splice_plus_retained_depth.matrix \
    --samples sample_metadata.tsv \
    --output_dir results \
    --contrast "perturb,control"
```

That's it! Results are in `results/` directory.

### Alternative: Count BAMs During the Pipeline

Instead of providing count and denominator matrices, `run_diff_splice_analysis.py`
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
    --contrast "perturb,control" \
    --site_depth_strand_mode RF
```

In this mode the pipeline writes intermediate matrices under
`results/workdir/bam_inputs/`. If `--site_depth_strand_mode` is omitted, depth
denominators are computed as unstranded.

### Offset Modes

The default mode uses splice-site boundary depth (splice counts + intron-body
retained depth):

```bash
--offset_mode splice_plus_retained   # default
```

A gene-scoped alternative adds gene-total junction evidence to the denominator:

```bash
--offset_mode splice_vs_rest \
--gtf annotation.gtf
```

`splice_vs_rest` computes `D = gene_total_junction_count + max(0, D^spr - Y_i)`,
combining all intron counts in the gene with the local spr remainder (competing
site splices and retained depth). Introns with no gene assignment fall back to
`splice_plus_retained` automatically.

Both modes use `max_splice_plus_retained_depth.matrix` as `--offset_matrix`;
the `splice_vs_rest` transformation is applied in-pipeline from the count matrix
and GTF.

### Statistical Modes

Two statistical engines are available via `--stat_mode`:

```bash
--stat_mode offset        # default: edgeR NB GLM with log-depth fixed offset
--stat_mode interaction   # DEXSeq-style stacked NB interaction model
```

**`offset` mode** (default): The log-transformed denominator is supplied as a
fixed edgeR offset. `logFC` approximates `log2(PSI_A / PSI_B)`.

**`interaction` mode**: Builds a doubled count matrix (focal + other
pseudo-samples per original sample) and fits
`~ sample_id + exon_type + group:exon_type` by LRT. The `sample_id` blocking
factor absorbs per-sample depth variation without a fixed offset. `logFC` is a
log2 focal/other odds-ratio change — numerically larger than the offset-mode
`logFC` for the same PSI shift when baseline PSI is high. `delta_PSI` is
unaffected by `--stat_mode` and remains the PSI-scale effect size.

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

1. **Load count and selected denominator matrices**
   - Count matrix rows are introns
   - Offset matrix has the same introns and samples
   - Offsets are raw denominators, not log-transformed

2. **Annotate introns** (if GTF provided)
   - Maps introns to gene names
   - Marks known vs novel intron status

3. **Filter low-confidence features before statistical testing**
   - Keep canonical splice sites only (GT-AG, GC-AG, AT-AC)
   - Filter introns with insufficient read support
   - Filter introns with insufficient selected-denominator support
   - Compute PSI and filter by minimum `|delta_PSI|` for the requested contrast

4. **Prepare model inputs**
   - Creates count, offset, and annotation files
   - Log-transform selected denominators for GLM

5. **Run statistical analysis** (controlled by `--stat_mode`)
   - `offset` (default): edgeR NB QL GLM with log-depth fixed offset
   - `interaction`: DEXSeq-style doubled count matrix, LRT on `group:exon_type` interaction
   - Tests intron usage proportions, not expression
   - Each intron tested once

6. **Add PSI values to results**
   - Merges PSI values with model results
   - Reports model FDR over the prefiltered tested intron set
   - Does not recompute FDR after statistical testing

### Running Individual Modules

The main pipeline is the recommended interface. Matrix mode consumes existing
count and offset matrices; BAMs are only used in BAM-manifest mode.

Useful lower-level steps are:

```bash
# Run edgeR on prepared matrices
Rscript util/run_edgeR_analysis.R \
    --counts edgeR_input.counts.tsv \
    --offsets edgeR_input.offsets.tsv \
    --annotations edgeR_input.annotations.tsv \
    --samples sample_metadata.tsv \
    --output results \
    --contrast "perturb,control"

# Compute PSI summaries from prepared matrices
python3 util/compute_psi.py \
    --counts edgeR_input.counts.tsv \
    --annotations edgeR_input.annotations.tsv \
    --samples sample_metadata.tsv \
    --results results.intron_results.tsv \
    --output results_with_psi.tsv
```

## Output Files

### Key Results Files

**Intron-level results:**
- `edgeR_results.all.tsv` - All tested introns with edgeR statistics and PSI values
- `edgeR_results.significant_introns.tsv` - Significant introns among the prefiltered tested set, using edgeR's FDR and the configured `--min_logFC`.

**Intermediate files:**
- `workdir/introns_filtered.tsv` - Counts and annotations after count, denominator-depth, and delta-PSI prefilters
- `workdir/site_depth_offsets.filtered.tsv` - Filtered raw denominators used for PSI
- `workdir/edgeR_input.counts.tsv` - Filtered count matrix used by edgeR
- `workdir/edgeR_input.offsets.tsv` - Log-transformed selected denominators used by edgeR
- `workdir/edgeR_input.annotations.tsv` - Intron annotations used by edgeR
- `workdir/edgeR_results.intron_results.tsv` - Raw edgeR output before PSI columns are added
- `workdir/psi.psi_values.tsv` - Per-sample PSI values, group means, and delta PSI

**Diagnostics:**
- `workdir/edgeR_results.diagnostics.pdf` - BCV, dispersion, MA, volcano plots

### Result Interpretation

**Intron-level columns:**
- `intron_id`: Intron coordinates and splice sites (chr:start-end^splice_pair^flag)
- `chr`, `start`, `end`, `strand`, `donor`, `acceptor`: Parsed intron coordinates
- `splice_pair`, `splice_flag`: Splice motif and canonical flag
- `offset_mode`: selected offset mode (`splice_plus_retained` or `splice_vs_rest`)
- `offset_source`: denominator source; `splice_vs_rest` introns without a gene assignment show `splice_plus_retained`
- `stat_mode`: statistical engine used (`offset` or `interaction`)
- `gene_name`: Gene name (if GTF provided)
- `intron_status`: known/novel (if GTF provided)
- `logFC`: Model-based log2 effect size.
  In **offset mode**: ≈ `log2(PSI_A / PSI_B)`, estimated from counts with depth offsets and dispersion modeling.
  In **interaction mode**: log2 focal/other odds-ratio change; larger in magnitude than offset logFC for the same PSI shift when baseline PSI is high.
  Positive = increased usage in the first group of the contrast.
- `logCPM`: Average log2 counts per million (focal pseudo-samples in interaction mode)
- `F` / `LR`: F-statistic (offset mode QL F-test) or likelihood-ratio statistic (interaction mode LRT)
- `PValue`: p-value from the selected statistical test
- `FDR`: Benjamini-Hochberg adjusted p-value over introns that passed pre-test filters
- `*_mean_PSI`: Mean PSI in each group (if PSI computed)
- `delta_PSI`: Difference in mean PSI between groups (if PSI computed)

**Important notes:**
- Each canonical intron is tested **once** with the selected denominator
- In **offset mode**, `logFC` represents relative change in usage after accounting for the selected depth denominator; it approximates `log2(PSI_A/PSI_B)`
- In **interaction mode**, `logFC` is a log2 odds ratio; use `delta_PSI` as the PSI-scale effect size regardless of stat mode
- `--min_delta_psi` is applied before statistical testing, so FDR is not recomputed after results are generated

## Parameter Tuning

### Default Parameters (balanced)
```bash
--min_intron_count 10 \
--min_intron_samples 2 \
--min_offset_depth 20 \
--min_offset_samples 3
```

### For Long-Read Data (more lenient)
```bash
--min_intron_count 5 \
--min_intron_samples 2 \
--min_offset_depth 10 \
--min_offset_samples 2
```

### For Conservative Analysis (fewer false positives)
```bash
--min_intron_count 20 \
--min_offset_depth 50 \
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

Run one explicit contrast per invocation:

```bash
python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --offset_matrix data/intron_counts.max_adjacent_depth.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/perturb_vs_control \
    --contrast perturb,control \
    --offset_mode exon_adjacent_depth
```

### Custom Parameters
```bash
python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --offset_matrix data/intron_counts.max_adjacent_depth.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/custom \
    --min_intron_count 20 \
    --min_offset_depth 30 \
    --min_logFC 0.5 \
    --fdr_threshold 0.05 \
    --contrast perturb,control \
    --offset_mode exon_adjacent_depth \
    --batch_col batch
```

## Testing

### Quick Integration Test
Run automated tests to verify the pipeline works correctly:
```bash
make test
```

This runs a quick test (~1-2 minutes) using a small dataset (550 introns). It validates:
- Offset matrix input handling
- Count, denominator-depth, and delta-PSI prefiltering
- statistical analysis
- PSI calculation with selected denominators
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
make -C testing test_modes
make -C testing test_viz
make -C testing clean
```

The shell scripts in `testing/` are used by the local Makefile and are not the recommended public entrypoint.

See [testing/README.md](testing/README.md) for fixture details and expected outputs.

## Design Principles

This pipeline implements specific design choices:

1. **Intron-level analysis** - Each intron is tested once
2. **Depth denominators** - Usage is modeled relative to selected local read depth
3. **Counts remain counts** - No PSI/TPM transformation before statistical testing; edgeR modes normalize with GLM offsets
4. **Pre-test effect-size filtering** - `--min_delta_psi` defines the tested intron set
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
