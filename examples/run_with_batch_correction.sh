#!/bin/bash

# Example with batch correction

INTRON_MATRIX="data/intron_counts.matrix"
OFFSET_MATRIX="data/intron_counts.max_splice_plus_retained_depth.matrix"
SAMPLE_METADATA="examples/sample_metadata_with_batch.tsv"
OUTPUT_DIR="results/diff_splice_with_batch"

python3 DSF.py \
    --matrix ${INTRON_MATRIX} \
    --offset_matrix ${OFFSET_MATRIX} \
    --samples ${SAMPLE_METADATA} \
    --output_dir ${OUTPUT_DIR} \
    --offset_mode splice_plus_retained \
    --min_intron_count 10 \
    --min_intron_samples 2 \
    --min_offset_depth 20 \
    --min_offset_samples 3 \
    --group_col group \
    --batch_col batch \
    --contrast "perturb,control" \
    --fdr_threshold 0.05

echo "Analysis with batch correction complete!"
