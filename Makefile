.PHONY: test test-quick test-full clean-test help check-deps

# Default target
help:
	@echo "Diff-Splice-Finder Testing Targets:"
	@echo ""
	@echo "  make test         - Run quick integration test (~1-2 min)"
	@echo "  make test-quick   - Same as 'make test'"
	@echo "  make test-full    - Run full integration test with all features"
	@echo "  make clean-test   - Clean all test output directories"
	@echo "  make check-deps   - Check Python and R dependencies"
	@echo ""
	@echo "Test files are in testing/ directory"

# Check dependencies
check-deps:
	@echo "Checking Python dependencies..."
	@python3 -c "import pandas; import numpy; import scipy" 2>/dev/null && echo "  ✓ Python packages OK" || (echo "  ✗ Missing Python packages. Run: pip install -r requirements.txt" && exit 1)
	@echo "Checking R dependencies..."
	@Rscript -e "library(edgeR); library(optparse)" 2>/dev/null && echo "  ✓ R packages OK" || (echo "  ✗ Missing R packages. See README for installation" && exit 1)
	@echo "All dependencies satisfied!"

# Quick test using small dataset (~1-2 minutes)
test: check-deps
	@echo "Running quick integration test..."
	@cd testing && ./run_quick_test.sh
	@echo ""
	@echo "✓ Quick test passed!"

test-quick: test

# Full integration test with all features
test-full: check-deps
	@echo "Running full integration test with all features..."
	@cd testing && bash -c '\
		../run_diff_splice_analysis.py \
			--matrix test_intron_counts.matrix \
			--samples test_metadata_control.tsv \
			--output_dir full_test_output \
			--gtf test_annotation.gtf \
			--contrast "TDP43-control" \
			--min_intron_count 5 \
			--min_intron_samples 2 \
			--min_cluster_count 10 \
			--min_cluster_samples 2 \
			--min_delta_psi 0.1 \
			--fdr_threshold 0.05 \
			--cpu 2'
	@echo ""
	@echo "Validating full test output..."
	@test -f testing/full_test_output/introns_clustered.tsv || (echo "✗ Missing introns_clustered.tsv" && exit 1)
	@test -f testing/full_test_output/shared_offsets.raw_cluster_totals.tsv || (echo "✗ Missing shared_offsets" && exit 1)
	@test -f testing/full_test_output/edgeR_results.intron_results.tsv || (echo "✗ Missing edgeR results" && exit 1)
	@test -f testing/full_test_output/psi.psi_values.tsv || (echo "✗ Missing PSI values" && exit 1)
	@test -f testing/full_test_output/edgeR_results.intron_results_with_psi.tsv || (echo "✗ Missing unfiltered PSI results" && exit 1)
	@test -f testing/full_test_output/edgeR_results.intron_results_with_psi.psi_filtered.tsv || (echo "✗ Missing PSI filtered results" && exit 1)
	@test -f testing/full_test_output/edgeR_results.diagnostics.pdf || (echo "✗ Missing diagnostics PDF" && exit 1)
	@echo "  ✓ All expected output files present"
	@echo ""
	@echo "Checking for gene annotations..."
	@grep -q "gene_name" testing/full_test_output/edgeR_results.intron_results.tsv && echo "  ✓ Gene annotations present" || echo "  ⚠ Gene annotations not found (may be expected)"
	@echo ""
	@echo "Checking cluster columns..."
	@grep -q "donor_cluster" testing/full_test_output/introns_clustered.tsv && echo "  ✓ donor_cluster column present" || (echo "  ✗ Missing donor_cluster" && exit 1)
	@grep -q "acceptor_cluster" testing/full_test_output/introns_clustered.tsv && echo "  ✓ acceptor_cluster column present" || (echo "  ✗ Missing acceptor_cluster" && exit 1)
	@echo ""
	@echo "✓ Full integration test passed!"

# Clean test outputs
clean-test:
	@echo "Cleaning test output directories..."
	@rm -rf testing/quick_test_output
	@rm -rf testing/full_test_output
	@rm -rf testing/psi_demo/output
	@rm -rf testing/psi_filtered_test/output
	@rm -rf testing/psi_integration_test/output
	@rm -rf testing/gene_annotation_test/output
	@rm -rf testing/psi_column_test/output
	@echo "✓ Test outputs cleaned"
