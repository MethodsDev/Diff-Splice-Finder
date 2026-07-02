# Testing Directory

## Quick Test Dataset

### Files
- **test_intron_counts.matrix**: Small test dataset with 550 introns
  - Subset of the full data/intron_counts.matrix (970K introns)
  - Contains introns with reasonable counts and variation
  - Includes some low-count introns for filter testing
  - Runs in ~1-2 minutes instead of hours

- **test_metadata_control.tsv**: Sample metadata for testing
  - 3 perturb samples
  - 3 control samples

- **test_site_depth_offsets.matrix**: Synthetic site-depth denominator matrix
  - Same shape as `test_intron_counts.matrix`
  - Used to exercise the streamlined site-depth offset pipeline

- **test_annotation.gtf**: GTF annotation for gene annotation testing
- **test_annotation.intron_cache.tsv**: Pre-computed intron cache

- **Makefile**: Preferred testing entrypoint
  - `make test` runs the BAM intron-count smoke test and Diff-Splice-Finder quick test
  - `make test_viz` runs the DEXSeq-like PDF smoke test
  - `make clean` removes local test outputs

- **run_quick_test.sh** / **run_plot_quick_test.sh**: Implementation scripts used by the Makefile targets
- **run_site_depth_strand_test.py**: Synthetic BAM checks for strand-specific depth and paired-end overlap handling
- **run_bam_introns_test.py**: Uses `data/alignments.b38.sorted.bam` to generate and validate a `.introns` file

### Usage

From the repository root:
```bash
make test
make test-viz
```

Or from the testing directory:
```bash
cd testing
make test
make test_viz
make clean
```

This will create `quick_test_output/` with results.
The `test` target also runs the synthetic site-depth strand checks and the
checked-in BAM intron-count smoke test before the pipeline smoke test.

### BAM Intron-Count Smoke Test

`run_bam_introns_test.py` exercises the per-sample BAM-to-`.introns` stage using
`data/alignments.b38.sorted.bam`. The test creates a temporary minimal `chr7`
FASTA with canonical `GT--AG` splice motifs at the BAM's discovered junctions,
runs `util/count_introns_from_bam.py`, and writes:

```text
bam_introns_test_output/alignments.b38.introns
```

The generated file has the same row format used before matrix construction:

```text
intron    splice_pair    splice_flag    count    site_depth_offset
```

The test validates that 49 introns are reported and that the known junction
`chr7:55155947-55156532` has count 342 with a positive site-depth offset.

### How the Quick Test Inputs Work

The quick test runs the pipeline in **matrix mode**. It does not generate
per-sample `.introns` files from BAMs. Instead, it starts from prebuilt intron x
sample matrices:

- `test_intron_counts.matrix` supplies the read counts.
- `test_site_depth_offsets.matrix` supplies the raw site-depth denominators.
- `test_metadata_control.tsv` defines the sample order and groups.

The six sample columns are:

```text
perturb_bc04, perturb_bc05, perturb_bc06
control_bc01, control_bc02, control_bc03
```

In this fixture, the site-depth offset matrix is intentionally synthetic: every
cell is exactly the count plus 20. For each intron/sample pair:

```text
offset = count + 20
PSI = count / offset
```

For example, if `perturb_bc04` has count 11 for an intron, its offset is 31 and
its PSI is `11 / 31 = 0.3548`.

The pipeline aligns the count and offset matrices to the samples listed in
`test_metadata_control.tsv`, filters introns, computes per-sample PSI values,
and only then runs edgeR. The main pre-edgeR intermediate files are:

- `workdir/introns_filtered.tsv` - filtered intron annotations plus sample counts
- `workdir/site_depth_offsets.filtered.tsv` - raw offsets for those filtered introns
- `workdir/psi.psi_values.tsv` - per-sample PSI, group summaries, and delta PSI
- `workdir/edgeR_input.counts.tsv` - count matrix passed to edgeR
- `workdir/edgeR_input.offsets.tsv` - `log(offset + 0.5)` matrix passed to edgeR
- `workdir/edgeR_input.annotations.tsv` - intron annotations passed to edgeR

### Where Per-Sample Outputs Fit In

Per-sample files are produced only in **BAM-manifest mode**, when
`run_diff_splice_analysis.py` is called without `--matrix` and with a sample
manifest plus `--genome_fa`. That mode runs two BAM-counting passes:

1. Discovery pass:
   `workdir/bam_inputs/discovery_introns/<sample>.introns`
2. Targeted pass:
   `workdir/bam_inputs/targeted_introns/<sample>.introns`

Each `.introns` file contains one sample's intron-level records:

```text
intron    splice_pair    splice_flag    count    site_depth_offset
```

Those per-sample files are then combined into:

- `workdir/bam_inputs/intron_counts.matrix`
- `workdir/bam_inputs/intron_counts.offsets.matrix`

The quick test skips this BAM-to-matrix stage by providing
`test_intron_counts.matrix` and `test_site_depth_offsets.matrix` directly.

### Test Dataset Statistics

- **Size**: 551 lines (550 introns + header) vs 970K in full dataset
- **File size**: 28KB vs ~100MB+ for full dataset
- **Samples**: 6 (perturb_bc04-06, control_bc01-03)
- **Chromosomes**: chr1-22, chrX, chrY
- **Count range**: 0-164
- **Mean counts**: ~28 reads per intron

### Expected Outputs

After running the quick test, you should see:
- `edgeR_results.all.tsv` - Main full result table
- `edgeR_results.significant_introns.tsv` - Significant hits from tested introns
- `workdir/introns_filtered.tsv` - Introns passing count, site-depth, and delta-PSI prefilters
- `workdir/site_depth_offsets.filtered.tsv` - Raw site-depth offsets for tested introns
- `workdir/psi.psi_values.tsv` - PSI values with site-depth denominators
- `workdir/edgeR_results.diagnostics.pdf` - QC plots
- `bam_introns_test_output/alignments.b38.introns` - BAM-derived per-sample intron file
- `plot_quick_test_output/TESTGENE_perturb_vs_control_intron_DEXseq_like.pdf` - Plotting utility output

### Validating Results

Key things to check:
1. `workdir/edgeR_input.offsets.tsv` contains log-transformed site-depth offsets
2. PSI values use the raw site-depth offsets as denominators
3. `--min_delta_psi` filtering happens before edgeR
4. Each intron is tested once
5. The synthetic strand test confirms F/R/FR/RF orientation handling and paired-mate overlap de-duplication
6. `bam_introns_test_output/alignments.b38.introns` contains the expected BAM-derived junctions

### Creating Custom Test Datasets

To create a different test subset:

```python
import pandas as pd

# Load full matrix
df = pd.read_csv("../data/intron_counts.matrix", sep="\t", index_col=0)

# Apply your filters (e.g., introns with good counts)
df = df[df.sum(axis=1) >= 100]

# Sample and save
subset = df.sample(n=500, random_state=42)
subset.to_csv("custom_test.matrix", sep="\t")
```

## Other Test Files

- **data/**: Test BAM files for intron counting validation
  - `alignments.b38.sorted.bam`
  - `alignments.b38.sorted.bam.bai`

## Cleaning Test Outputs

```bash
make clean-test
# or:
make -C testing clean
```

This removes all test output directories.
