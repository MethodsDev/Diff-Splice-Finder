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
  - `make test` runs the Diff-Splice-Finder quick test
  - `make test_viz` runs the DEXSeq-like PDF smoke test
  - `make clean` removes local test outputs

- **run_quick_test.sh** / **run_plot_quick_test.sh**: Implementation scripts used by the Makefile targets
- **run_site_depth_strand_test.py**: Synthetic BAM checks for strand-specific depth and paired-end overlap handling

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
The `test` target also runs the synthetic site-depth strand checks before the
pipeline smoke test.

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
- `plot_quick_test_output/TESTGENE_perturb_vs_control_intron_DEXseq_like.pdf` - Plotting utility output

### Validating Results

Key things to check:
1. `workdir/edgeR_input.offsets.tsv` contains log-transformed site-depth offsets
2. PSI values use the raw site-depth offsets as denominators
3. `--min_delta_psi` filtering happens before edgeR
4. Each intron is tested once
5. The synthetic strand test confirms F/R/FR/RF orientation handling and paired-mate overlap de-duplication

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
