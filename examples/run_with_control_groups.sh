#!/bin/bash

# Example: Compare all treatment groups against control samples
# This script demonstrates how to use the --control_groups parameter
# to perform control-based comparisons instead of all pairwise comparisons

./run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir output_with_controls \
    --group_col group \
    --control_groups control \
    --min_delta_psi 0.1 \
    --fdr_threshold 0.05 \
    --cpu 4

# Notes:
# - All non-control groups (e.g., TDP43) will be compared against the control group
# - If you have multiple control types, use: --control_groups control,wildtype
# - This reduces the number of comparisons and focuses on treatment vs control
# - Results will show: TDP43 vs control (instead of all pairwise combinations)
