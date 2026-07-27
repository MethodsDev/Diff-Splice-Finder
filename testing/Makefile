.PHONY: help test test_modes test_intx test_viz test_wdl clean

help:
	@echo "Testing targets:"
	@echo ""
	@echo "  make test      - Run the Diff-Splice-Finder quick test"
	@echo "  make test_modes - Run the offset-mode refactor fixture test"
	@echo "  make test_intx - Run interaction-engine fixture tests"
	@echo "  make test_viz  - Run the DEXSeq-like PDF visualization test"
	@echo "  make test_wdl  - Run the WDL count_introns execution test (needs docker + miniwdl)"
	@echo "  make clean     - Remove test outputs via cleanMe.sh"

test:
	@python3 run_site_depth_strand_test.py
	@python3 run_bam_introns_test.py
	@./run_quick_test.sh

test_modes:
	@./run_mode_refactor_test.sh

test_intx:
	@./run_intx_engine_test.sh

test_wdl:
	@python3 run_wdl_count_introns_test.py

test_viz:
	@./run_plot_quick_test.sh

clean:
	@bash ./cleanMe.sh
