#!/bin/bash

# Example: run one explicit treatment-vs-control contrast.

python3 run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --site_depth_offsets data/intron_counts.offsets.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir output_perturb_vs_control \
    --group_col group \
    --contrast perturb,control \
    --min_delta_psi 0.1 \
    --fdr_threshold 0.05
