#!/bin/bash

set -ex


rm -rf ./test_out_dir

python3 ~/GITHUB/MDL/Diff-Splice-Finder/run_diff_splice_analysis.py \
    --matrix test.intron_counts.matrix.gz \
    --samples sample_metadata.tsv \
    --output_dir test_out_dir \
    --min_intron_count 20 \
    --min_cluster_count 30 \
    --min_logFC 0.5 \
    --fdr_threshold 0.05 \
    --gtf $CTAT_GENOME_LIB/ref_annot.gtf \
    --min_delta_psi 0.1


