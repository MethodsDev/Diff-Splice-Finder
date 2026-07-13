#!/usr/bin/env bash

# Quick test script using the small test dataset
# This runs in ~1-2 minutes instead of hours
# Tests intron-level analysis with matrix-supplied depth denominators

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
INPUT_DIR="$HERE/mode_refactor_inputs"
OUT_DIR="$HERE/quick_test_output"

echo "=== Quick Test with Small Dataset ==="
echo "Dataset: mode_refactor fixture introns"
echo "Mode: Intron-level analysis with splice-plus-retained depth offsets"
echo ""

# Run the analysis with the mode_refactor fixture matrices
"$ROOT/DSF.py" \
    --matrix "$INPUT_DIR/intron_counts.matrix" \
    --offset_matrix "$INPUT_DIR/intron_counts.max_splice_plus_retained_depth.matrix" \
    --samples "$INPUT_DIR/sample_metadata.tsv" \
    --output_dir "$OUT_DIR" \
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
ls -lh "$OUT_DIR"/edgeR_results*.tsv 2>/dev/null || echo "  (check quick_test_output/ directory)"
echo ""
echo "Checking key files exist:"
test -f "$OUT_DIR/edgeR_results.all.tsv" && echo "  ✓ edgeR_results.all.tsv"
test -f "$OUT_DIR/edgeR_results.significant_introns.tsv" && echo "  ✓ edgeR_results.significant_introns.tsv"
test -f "$OUT_DIR/workdir/introns_filtered.tsv" && echo "  ✓ workdir/introns_filtered.tsv"
test -f "$OUT_DIR/workdir/site_depth_offsets.filtered.tsv" && echo "  ✓ workdir/site_depth_offsets.filtered.tsv"
test -f "$OUT_DIR/workdir/psi.psi_values.tsv" && echo "  ✓ workdir/psi.psi_values.tsv"
echo ""
echo "Test passed!"
