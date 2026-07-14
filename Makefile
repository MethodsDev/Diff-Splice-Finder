.PHONY: test test-quick test-full test-viz test-modes test-intx clean-test help check-deps check-dexseq

# Default target
help:
	@echo "Diff-Splice-Finder Testing Targets:"
	@echo ""
	@echo "  make test         - Run quick integration test (~1-2 min)"
	@echo "  make test-modes   - Run offset-mode refactor fixture test"
	@echo "  make test-intx    - Run interaction-engine fixture tests"
	@echo "  make test-viz     - Run DEXSeq-like PDF visualization test"
	@echo "  make test-quick   - Same as 'make test'"
	@echo "  make test-full    - Run full integration test with all features"
	@echo "  make clean-test   - Clean all test output directories"
	@echo "  make check-deps   - Check Python and R dependencies"
	@echo "  make check-dexseq - Check optional DEXSeq dependency"
	@echo ""
	@echo "Test files are in testing/ directory"
	@echo "You can also run 'make -C testing help'"

# Check dependencies
check-deps:
	@echo "Checking Python dependencies..."
	@python3 -c "import pandas; import numpy; import scipy" 2>/dev/null && echo "  ✓ Python packages OK" || (echo "  ✗ Missing Python packages. Run: pip install -r requirements.txt" && exit 1)
	@echo "Checking R dependencies..."
	@Rscript -e "library(edgeR); library(optparse)" 2>/dev/null && echo "  ✓ R packages OK" || (echo "  ✗ Missing R packages. See README for installation" && exit 1)
	@echo "All dependencies satisfied!"

check-dexseq:
	@echo "Checking optional DEXSeq dependency..."
	@Rscript -e "library(DEXSeq)" 2>/dev/null && echo "  ✓ DEXSeq OK" || (echo "  ✗ Missing optional DEXSeq package. Run: BiocManager::install('DEXSeq')" && exit 1)

# Quick test using small dataset (~1-2 minutes)
test: check-deps
	@echo "Running quick integration test..."
	@$(MAKE) -C testing test
	@echo ""
	@echo "✓ Quick test passed!"

test-quick: test

test-viz: check-deps
	@echo "Running DEXSeq-like PDF visualization test..."
	@$(MAKE) -C testing test_viz
	@echo ""
	@echo "✓ Visualization test passed!"

test-modes: check-deps
	@echo "Running offset-mode refactor fixture test..."
	@$(MAKE) -C testing test_modes
	@echo ""
	@echo "✓ Offset-mode fixture test passed!"

test-intx: check-deps
	@echo "Running interaction-engine fixture tests..."
	@$(MAKE) -C testing test_intx
	@echo ""
	@echo "✓ Interaction-engine fixture tests passed!"

# Full integration test with all features
test-full: check-deps
	@echo "Running full integration test with all features..."
	@cd testing && bash -c '\
		../DSF.py \
			--matrix test_intron_counts.matrix \
			--site_depth_offsets test_site_depth_offsets.matrix \
			--samples test_metadata_control.tsv \
			--output_dir full_test_output \
			--gtf test_annotation.gtf \
			--contrast "perturb,control" \
			--min_intron_count 5 \
			--min_intron_samples 2 \
			--min_offset_depth 10 \
			--min_offset_samples 2 \
			--min_delta_psi 0.1 \
			--fdr_threshold 0.05'
	@echo ""
	@echo "Validating full test output..."
	@test -f testing/full_test_output/workdir/introns_filtered.tsv || (echo "✗ Missing workdir/introns_filtered.tsv" && exit 1)
	@test -f testing/full_test_output/workdir/site_depth_offsets.filtered.tsv || (echo "✗ Missing workdir/site_depth_offsets.filtered.tsv" && exit 1)
	@test -f testing/full_test_output/edgeR_results.all.tsv || (echo "✗ Missing full PSI results" && exit 1)
	@test -f testing/full_test_output/edgeR_results.significant_introns.tsv || (echo "✗ Missing significant results" && exit 1)
	@test -f testing/full_test_output/workdir/edgeR_results.intron_results.tsv || (echo "✗ Missing workdir edgeR results" && exit 1)
	@test -f testing/full_test_output/workdir/psi.psi_values.tsv || (echo "✗ Missing workdir PSI values" && exit 1)
	@test -f testing/full_test_output/workdir/edgeR_results.diagnostics.pdf || (echo "✗ Missing workdir diagnostics PDF" && exit 1)
	@echo "  ✓ All expected output files present"
	@echo ""
	@echo "Checking for gene annotations..."
	@grep -q "gene_name" testing/full_test_output/edgeR_results.all.tsv && echo "  ✓ Gene annotations present" || echo "  ⚠ Gene annotations not found (may be expected)"
	@echo ""
	@echo "Checking site-depth annotation columns..."
	@grep -q "offset_source" testing/full_test_output/workdir/edgeR_input.annotations.tsv && echo "  ✓ offset_source column present" || (echo "  ✗ Missing offset_source" && exit 1)
	@echo ""
	@echo "✓ Full integration test passed!"

# Clean test outputs
clean-test:
	@echo "Cleaning test output directories..."
	@$(MAKE) -C testing clean
	@rm -rf testing/full_test_output
	@rm -rf testing/psi_demo/output
	@rm -rf testing/psi_filtered_test/output
	@rm -rf testing/psi_integration_test/output
	@rm -rf testing/gene_annotation_test/output
	@rm -rf testing/psi_column_test/output
	@echo "✓ Test outputs cleaned"
