#!/bin/bash

# Example script to run differential splicing analysis
# This demonstrates the full pipeline workflow

# Input files
INTRON_MATRIX="data/intron_counts.matrix"
OFFSET_MATRIX="data/intron_counts.max_adjacent_depth.matrix"
SAMPLE_METADATA="examples/sample_metadata.tsv"
OUTPUT_DIR="results/diff_splice_analysis"

# Run the full pipeline
# Note: Pipeline automatically resumes if interrupted
# Use --force_rerun to start from scratch
python3 run_diff_splice_analysis.py \
    --matrix ${INTRON_MATRIX} \
    --offset_matrix ${OFFSET_MATRIX} \
    --samples ${SAMPLE_METADATA} \
    --output_dir ${OUTPUT_DIR} \
    --offset_mode exon_adjacent_depth \
    --min_intron_count 10 \
    --min_intron_samples 2 \
    --min_offset_depth 20 \
    --min_offset_samples 3 \
    --group_col group \
    --contrast "perturb,control" \
    --fdr_threshold 0.05 \
    --min_logFC 0.5

echo "Analysis complete! Results in: ${OUTPUT_DIR}"
echo ""
echo "To rerun from scratch: add --force_rerun flag"
