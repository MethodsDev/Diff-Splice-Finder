# Testing Directory

## Quick Test Dataset

### Files
- **test_intron_counts.matrix**: Small test dataset with 550 introns
  - Subset of the full data/intron_counts.matrix (970K introns)
  - Contains introns with reasonable counts and variation
  - Includes some low-count introns for filter testing
  - Runs in ~1-2 minutes instead of hours

- **test_metadata_control.tsv**: Sample metadata for testing
  - 3 TDP43 samples
  - 3 control samples

- **test_annotation.gtf**: GTF annotation for gene annotation testing
- **test_annotation.intron_cache.tsv**: Pre-computed intron cache

- **run_quick_test.sh**: Quick test script
  - Tests the intron-level analysis with shared offsets
  - Uses control groups feature
  - Good for rapid testing of code changes

### Usage

From the repository root:
```bash
make test
```

Or from the testing directory:
```bash
cd testing
./run_quick_test.sh
```

This will create `quick_test_output/` with results.

### Running Full Integration Test

Test all features including gene annotation and PSI filtering:
```bash
make test-full
```

### Test Dataset Statistics

- **Size**: 551 lines (550 introns + header) vs 970K in full dataset
- **File size**: 28KB vs ~100MB+ for full dataset
- **Samples**: 6 (TDP43_bc04-06, control_bc01-03)
- **Chromosomes**: chr1-22, chrX, chrY
- **Count range**: 0-164
- **Mean counts**: ~28 reads per intron

### Expected Outputs

After running the quick test, you should see:
- `introns_clustered.tsv` - Both donor_cluster and acceptor_cluster columns
- `shared_offsets.raw_cluster_totals.tsv` - Max(donor, acceptor) offsets
- `introns_filtered.tsv` - Filtered by both cluster thresholds
- `edgeR_results.intron_results.tsv` - Main statistical results
- `psi.psi_values.tsv` - PSI values with shared denominators
- `edgeR_results.intron_results_with_psi.tsv` - Combined results
- `edgeR_results.diagnostics.pdf` - QC plots

### Validating Results

Key things to check:
1. Both `donor_cluster` and `acceptor_cluster` columns exist in clustered file
2. Shared offsets file contains max(donor, acceptor) for each intron
3. PSI values use shared cluster totals as denominators
4. Each intron tested once (not separately for donor and acceptor)
5. No singleton cluster artifacts (PSI should not be 1.0 for rare events)

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
```

This removes all test output directories.
