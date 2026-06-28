#!/bin/bash

# Example with batch correction

INTRON_MATRIX="data/intron_counts.matrix"
SITE_DEPTH_OFFSETS="data/intron_counts.offsets.matrix"
SAMPLE_METADATA="examples/sample_metadata_with_batch.tsv"
OUTPUT_DIR="results/diff_splice_with_batch"

python3 run_diff_splice_analysis.py \
    --matrix ${INTRON_MATRIX} \
    --site_depth_offsets ${SITE_DEPTH_OFFSETS} \
    --samples ${SAMPLE_METADATA} \
    --output_dir ${OUTPUT_DIR} \
    --min_intron_count 10 \
    --min_intron_samples 2 \
    --min_offset_depth 20 \
    --min_offset_samples 3 \
    --group_col group \
    --batch_col batch \
    --contrast "perturb,control" \
    --fdr_threshold 0.05

echo "Analysis with batch correction complete!"
