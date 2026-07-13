# Testing Directory

## Test Fixtures

### Files
- **mode_refactor_inputs/**: Tiny matrix-mode fixture used by the current quick
  test and offset-mode tests.

- **test_intron_counts.matrix**, **test_metadata_control.tsv**, and
  **test_site_depth_offsets.matrix**: Older 550-intron matrix fixture retained
  for compatibility and plotting smoke tests.

- **test_annotation.gtf**: GTF annotation for gene annotation testing
- **test_annotation.intron_cache.tsv**: Pre-computed intron cache

- **Makefile**: Preferred testing entrypoint
  - `make test` runs the BAM intron-count smoke test and Diff-Splice-Finder quick test
  - `make test_modes` runs the offset-mode refactor fixture test
  - `make test_viz` runs the DEXSeq-like PDF smoke test
  - `make clean` removes local test outputs

- **run_quick_test.sh** / **run_plot_quick_test.sh**: Implementation scripts used by the Makefile targets
- **run_mode_refactor_test.sh**: Runs both `--offset_mode` execution paths from a tiny committed fixture
- **run_site_depth_strand_test.py**: Synthetic BAM checks for strand-specific depth and paired-end overlap handling
- **run_bam_introns_test.py**: Uses `data/alignments.b38.sorted.bam` to generate and validate a `.introns` file

### Usage

From the repository root:
```bash
make test
make test-modes
make test-viz
```

Or from the testing directory:
```bash
cd testing
make test
make test_modes
make test_viz
make clean
```

This will create `quick_test_output/` with results.
The `test` target also runs the synthetic site-depth strand checks and the
checked-in BAM intron-count smoke test before the pipeline smoke test.

### Offset-Mode Refactor Fixture

`mode_refactor_inputs/` is a tiny committed fixture for exercising both
supported denominator modes without large reference files:

- `chrTiny.fa`: synthetic 1 kb reference sequence; intentionally not gzipped
- `source_bams/*.bam`: four tiny synthetic BAMs (`A1`, `A2`, `B1`, `B2`)
- `targeted_introns/*.introns`: per-sample intron files with all depth columns
- `intron_counts.matrix`: read-count matrix
- `intron_counts.max_splice_plus_retained_depth.matrix`: denominator for
  `splice_plus_retained`; it is also the starting denominator transformed by
  `splice_vs_rest`
- `chrTiny_fixture.gtf`: single-gene annotation so `splice_vs_rest` can compute
  gene-total junction counts instead of falling back

Run it with:

```bash
make test-modes
# or:
make -C testing test_modes
```

The runner writes outputs under `mode_refactor_inputs/mode_runs/`, which is
ignored by git.

The fixture runner validates both `splice_plus_retained` and `splice_vs_rest`
from the same count and `max_splice_plus_retained_depth` matrices.

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
intron    splice_pair    splice_flag    count    left_adjacent_depth    right_adjacent_depth    max_adjacent_depth    ...    max_splice_plus_retained_depth    site_depth_offset
```

`site_depth_offset` is retained as a compatibility alias for
`max_adjacent_depth`. The matrix builder writes additional denominator matrices
for numeric depth columns such as `max_adjacent_depth` and
`max_splice_plus_retained_depth`.

The test validates that 49 introns are reported and that the known junction
`chr7:55155947-55156532` has count 342 with a positive `max_adjacent_depth`
and compatibility `site_depth_offset`.

### How the Quick Test Inputs Work

The quick test runs the pipeline in **matrix mode**. It does not generate
per-sample `.introns` files from BAMs. Instead, it starts from prebuilt intron x
sample matrices:

- `mode_refactor_inputs/intron_counts.matrix` supplies the read counts.
- `mode_refactor_inputs/intron_counts.max_splice_plus_retained_depth.matrix`
  supplies the raw `splice_plus_retained` denominators.
- `mode_refactor_inputs/sample_metadata.tsv` defines the sample order and groups.

The four sample columns are:

```text
A1, A2, B1, B2
```

This fixture has three synthetic introns and exercises low-threshold filtering
so the full DSF path runs quickly. For each intron/sample pair:

```text
PSI = count / max_splice_plus_retained_depth
```

The pipeline aligns the count and offset matrices to the samples listed in
`sample_metadata.tsv`, filters introns, computes per-sample PSI values, and only
then runs statistical testing. The main pre-test intermediate files are:

- `workdir/introns_filtered.tsv` - filtered intron annotations plus sample counts
- `workdir/site_depth_offsets.filtered.tsv` - raw selected denominators for those filtered introns
- `workdir/psi.psi_values.tsv` - per-sample PSI, group summaries, and delta PSI
- `workdir/edgeR_input.counts.tsv` - count matrix passed to edgeR
- `workdir/edgeR_input.offsets.tsv` - `log(offset + 0.5)` matrix passed to edgeR
- `workdir/edgeR_input.annotations.tsv` - intron annotations passed to edgeR

### Where Per-Sample Outputs Fit In

Per-sample files are produced only in **BAM-manifest mode**, when
`DSF.py` is called without `--matrix` and with a sample
manifest plus `--genome_fa`. That mode runs two BAM-counting passes:

1. Discovery pass:
   `workdir/bam_inputs/discovery_introns/<sample>.introns`
2. Targeted pass:
   `workdir/bam_inputs/targeted_introns/<sample>.introns`

Each `.introns` file contains one sample's intron-level records:

```text
intron    splice_pair    splice_flag    count    left_adjacent_depth    right_adjacent_depth    max_adjacent_depth    ...    max_splice_plus_retained_depth    site_depth_offset
```

Those per-sample files are then combined into:

- `workdir/bam_inputs/intron_counts.matrix`
- `workdir/bam_inputs/intron_counts.offsets.matrix`
- `workdir/bam_inputs/intron_counts.max_adjacent_depth.matrix`
- `workdir/bam_inputs/intron_counts.max_splice_plus_retained_depth.matrix`

The quick test skips this BAM-to-matrix stage by providing the committed
`mode_refactor_inputs/` matrices directly.

### Legacy Matrix Fixture

The older 550-intron fixture remains available for ad hoc tests:

- **Size**: 551 lines (550 introns + header)
- **Samples**: 6 (perturb_bc04-06, control_bc01-03)
- **Files**: `test_intron_counts.matrix`, `test_site_depth_offsets.matrix`,
  `test_metadata_control.tsv`

`test_site_depth_offsets.matrix` is an older synthetic denominator matrix and is
not used by the current `run_quick_test.sh`.

### Expected Outputs

After running `make test`, you should see:
- `edgeR_results.all.tsv` - Main full result table
- `edgeR_results.significant_introns.tsv` - Significant hits from tested introns
- `workdir/introns_filtered.tsv` - Introns passing count, denominator-depth, and delta-PSI prefilters
- `workdir/site_depth_offsets.filtered.tsv` - Raw selected denominators for tested introns
- `workdir/psi.psi_values.tsv` - PSI values with selected denominators
- `workdir/edgeR_results.diagnostics.pdf` - QC plots
- `bam_introns_test_output/alignments.b38.introns` - BAM-derived per-sample intron file

After running `make test-viz`, you should see:

- `plot_quick_test_output/TESTGENE_perturb_vs_control_intron_DEXseq_like.pdf` - Plotting utility output

### Validating Results

Key things to check:
1. `workdir/edgeR_input.offsets.tsv` contains log-transformed selected denominators
2. PSI values use the raw selected denominators
3. `--min_delta_psi` filtering happens before statistical testing
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
