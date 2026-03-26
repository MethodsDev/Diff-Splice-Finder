#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
FIXTURE_DIR="${SCRIPT_DIR}/quick_test_output_local"
TEST_OUTPUT_DIR="${SCRIPT_DIR}/plot_quick_test_output"
TEST_PDF="${TEST_OUTPUT_DIR}/TESTGENE_perturb_vs_control_intron_DEXseq_like.pdf"

echo "=== Plot Quick Test ==="
echo "Fixture: testing/quick_test_output_local"
echo "Gene: TESTGENE"
echo ""

test -f "${FIXTURE_DIR}/workdir/edgeR_results.intron_results_with_psi.tsv"
test -f "${SCRIPT_DIR}/test_plot_gene.gtf"

rm -rf "${TEST_OUTPUT_DIR}"
mkdir -p "${TEST_OUTPUT_DIR}"

Rscript "${REPO_ROOT}/util/plot_intron_dexseq_like.R" \
    --output_dir "${FIXTURE_DIR}" \
    --gene TESTGENE \
    --gtf "${SCRIPT_DIR}/test_plot_gene.gtf" \
    --contrast perturb_vs_control \
    --output "${TEST_PDF}"

test -s "${TEST_PDF}"

echo ""
echo "Created:"
ls -lh "${TEST_PDF}"
echo ""
echo "Plot test passed!"
