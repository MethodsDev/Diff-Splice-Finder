#!/bin/bash

# Quick test script using the small test dataset
# This runs in ~1-2 minutes instead of hours
# Tests the refactored intron-level analysis with shared offsets

set -e

echo "=== Quick Test with Small Dataset ==="
echo "Dataset: 550 introns (vs 970K in full dataset)"
echo "Mode: Intron-level analysis with shared offsets"
echo ""

# Run the analysis with the small test matrix
../run_diff_splice_analysis.py \
    --matrix test_intron_counts.matrix \
    --samples test_metadata_control.tsv \
    --output_dir quick_test_output \
    --control_groups control \
    --min_intron_count 5 \
    --min_intron_samples 2 \
    --min_cluster_count 10 \
    --min_cluster_samples 2 \
    --fdr_threshold 0.05 \
    --cpu 2

echo ""
echo "=== Test Complete! ==="
echo "Results in: testing/quick_test_output/"
echo ""
echo "Key output files:"
ls -lh quick_test_output/edgeR_results*.tsv 2>/dev/null || echo "  (check quick_test_output/ directory)"
echo ""
echo "Checking key files exist:"
test -f quick_test_output/introns_clustered.tsv && echo "  ✓ introns_clustered.tsv"
test -f quick_test_output/shared_offsets.raw_cluster_totals.tsv && echo "  ✓ shared_offsets.raw_cluster_totals.tsv"
test -f quick_test_output/edgeR_results.intron_results.tsv && echo "  ✓ edgeR_results.intron_results.tsv"
test -f quick_test_output/psi.psi_values.tsv && echo "  ✓ psi.psi_values.tsv"
test -f quick_test_output/edgeR_results.intron_results_with_psi.tsv && echo "  ✓ edgeR_results.intron_results_with_psi.tsv"
echo ""
echo "Test passed!"
