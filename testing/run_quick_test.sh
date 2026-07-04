#!/bin/bash

# Quick test script using the small test dataset
# This runs in ~1-2 minutes instead of hours
# Tests intron-level analysis with matrix-supplied depth denominators

set -e

echo "=== Quick Test with Small Dataset ==="
echo "Dataset: mode_refactor fixture introns"
echo "Mode: Intron-level analysis with splice-plus-retained depth offsets"
echo ""

# Run the analysis with the mode_refactor fixture matrices
../run_diff_splice_analysis.py \
    --matrix mode_refactor_inputs/intron_counts.matrix \
    --offset_matrix mode_refactor_inputs/intron_counts.max_splice_plus_retained_depth.matrix \
    --samples mode_refactor_inputs/sample_metadata.tsv \
    --output_dir quick_test_output \
    --contrast A,B \
    --offset_mode splice_plus_retained \
    --min_intron_count 1 \
    --min_intron_samples 1 \
    --min_offset_depth 1 \
    --min_offset_samples 1 \
    --fdr_threshold 1 \
    --min_delta_psi 0

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
