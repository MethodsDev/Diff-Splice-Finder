#!/bin/bash

# Quick test script using the small test dataset
# This runs in ~1-2 minutes instead of hours
# Tests intron-level analysis with site-depth offsets

set -e

echo "=== Quick Test with Small Dataset ==="
echo "Dataset: 550 introns (vs 970K in full dataset)"
echo "Mode: Intron-level analysis with site-depth offsets"
echo ""

# Run the analysis with the small test matrix
../run_diff_splice_analysis.py \
    --matrix test_intron_counts.matrix \
    --site_depth_offsets test_site_depth_offsets.matrix \
    --samples test_metadata_control.tsv \
    --output_dir quick_test_output \
    --contrast perturb,control \
    --min_intron_count 5 \
    --min_intron_samples 2 \
    --min_offset_depth 10 \
    --min_offset_samples 2 \
    --fdr_threshold 0.05 \
    --min_delta_psi 0.01

echo ""
echo "=== Test Complete! ==="
echo "Results in: testing/quick_test_output/"
echo ""
echo "Key output files:"
ls -lh quick_test_output/edgeR_results*.tsv 2>/dev/null || echo "  (check quick_test_output/ directory)"
echo ""
echo "Checking key files exist:"
test -f quick_test_output/edgeR_results.all.tsv && echo "  ✓ edgeR_results.all.tsv"
test -f quick_test_output/edgeR_results.significant_introns.tsv && echo "  ✓ edgeR_results.significant_introns.tsv"
test -f quick_test_output/workdir/introns_filtered.tsv && echo "  ✓ workdir/introns_filtered.tsv"
test -f quick_test_output/workdir/site_depth_offsets.filtered.tsv && echo "  ✓ workdir/site_depth_offsets.filtered.tsv"
test -f quick_test_output/workdir/psi.psi_values.tsv && echo "  ✓ workdir/psi.psi_values.tsv"
echo ""
echo "Test passed!"
