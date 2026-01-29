# Resume Functionality

## Overview

The differential splicing analysis pipeline now supports **automatic resume** functionality. If the pipeline is interrupted (crash, timeout, manual termination), simply rerun the same command and it will pick up where it left off.

## How It Works

### Checkpointing System

The pipeline checks for the existence of output files before running each step:

1. **Clustering** - Checks for `{cluster_type}_clustered.tsv`
2. **Filtering** - Checks for `{cluster_type}_filtered.tsv`
3. **Offset computation** - Checks for all three files:
   - `{cluster_type}_edgeR_input.counts.tsv`
   - `{cluster_type}_edgeR_input.offsets.tsv`
   - `{cluster_type}_edgeR_input.annotations.tsv`
4. **edgeR analysis** - Checks for `{cluster_type}_edgeR_results.intron_results.tsv`
5. **Cluster aggregation** - Checks for `{cluster_type}_aggregated.cluster_results.tsv`

If an output file exists and is non-empty (size > 0), that step is skipped.

### Logging

When resuming, you'll see log messages like:
```
=== Clustering introns by donor ===
SKIPPING - Output already exists: results/donor/donor_clustered.tsv
```

This confirms the pipeline is resuming correctly.

## Usage Examples

### Normal Run (with automatic resume)

```bash
# First run - processes everything
python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/analysis

# If interrupted, just rerun the same command
python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/analysis
# Will automatically skip completed steps
```

### Force Complete Rerun

If you need to regenerate all outputs (e.g., after changing parameters or fixing bugs):

```bash
python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir results/analysis \
    --force_rerun
```

## Benefits

1. **Time savings** - No need to recompute expensive steps like edgeR analysis
2. **Robustness** - Pipeline can recover from crashes, timeouts, or cluster job limits
3. **Development friendly** - Test individual steps without rerunning everything
4. **Resource efficiency** - Avoid wasting compute resources on redundant calculations

## Limitations and Considerations

### When Resume Works Well

- Pipeline crashed unexpectedly
- Timeout on cluster job
- Manual interruption (Ctrl-C)
- Testing/debugging specific steps

### When to Use --force_rerun

- **Changed input files** - If you modified the count matrix or sample metadata
- **Changed parameters** - If you altered filtering thresholds, clustering method, etc.
- **Updated code** - If you fixed a bug or updated the analysis scripts
- **Corrupted outputs** - If intermediate files are incomplete or corrupted

### File Validation

The resume system only checks:
- File exists
- File size > 0 bytes

It does **NOT** validate:
- File contents are correct
- File is complete/uncorrupted
- File was generated with current parameters

If you're unsure about file integrity, use `--force_rerun`.

## Manual Checkpoint Management

You can manually remove specific checkpoint files to force rerun of specific steps:

```bash
# Rerun only the edgeR analysis for donor clusters
rm results/donor/donor_edgeR_results.intron_results.tsv
rm results/donor/donor_aggregated.cluster_results.tsv

# Then rerun pipeline (will skip clustering/filtering, rerun edgeR)
python3 run_diff_splice_analysis.py [same arguments]
```

## Integration with Cluster/HPC Systems

This resume functionality is particularly useful on cluster systems with:
- **Time limits** - Job may timeout before completion
- **Preemption** - Jobs may be killed and restarted
- **Resource constraints** - May need to split analysis across multiple jobs

Example SLURM workflow:
```bash
#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=32G

# Run analysis - will resume if resubmitted
python3 run_diff_splice_analysis.py \
    --matrix $DATA/intron_counts.matrix \
    --samples $DATA/samples.tsv \
    --output_dir $OUTPUT

# If job times out, simply resubmit:
# sbatch run_pipeline.sh
# Will automatically resume from last completed step
```

## Best Practices

1. **Use consistent commands** - Always use the same arguments when resuming
2. **Check logs** - Verify which steps were skipped/run
3. **Validate outputs** - After completion, check key output files
4. **Use --force_rerun after changes** - Parameter or code changes require fresh run
5. **Keep backups** - Consider backing up results before --force_rerun

## Troubleshooting

### Pipeline skips step but output is incomplete

**Solution**: Remove the checkpoint file and rerun:
```bash
rm results/donor/donor_edgeR_results.intron_results.tsv
python3 run_diff_splice_analysis.py [args]
```

### Want to restart specific cluster type

**Solution**: Remove that cluster's directory:
```bash
rm -rf results/donor/
python3 run_diff_splice_analysis.py [args]
# Will rerun donor, skip acceptor
```

### Need to change parameters

**Solution**: Use --force_rerun or new output directory:
```bash
# Option 1: Force rerun
python3 run_diff_splice_analysis.py [args] --force_rerun

# Option 2: New output directory
python3 run_diff_splice_analysis.py [args] --output_dir results/run2
```

## Implementation Details

The resume functionality is implemented in [run_diff_splice_analysis.py](../run_diff_splice_analysis.py):

- `file_exists_and_valid()` - Checks if file exists and is non-empty
- `run_command()` - Optional `skip_if_exists` parameter
- `--force_rerun` - Command-line flag to disable resume

All pipeline functions accept a `force_rerun` parameter that propagates through the execution chain.
