.PHONY: help test test_viz clean

help:
	@echo "Testing targets:"
	@echo ""
	@echo "  make test      - Run the Diff-Splice-Finder quick test"
	@echo "  make test_viz  - Run the DEXSeq-like PDF visualization test"
	@echo "  make clean     - Remove test outputs via cleanMe.sh"

test:
	@./run_quick_test.sh

test_viz:
	@./run_plot_quick_test.sh

clean:
	@bash ./cleanMe.sh
