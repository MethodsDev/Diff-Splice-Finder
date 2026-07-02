#!/bin/bash

set -euo pipefail

rm -rf ./quick_test_output
rm -rf ./plot_quick_test_output
rm -rf ./bam_introns_test_output
rm -rf ./mode_refactor_inputs/mode_runs
rm -f ./mode_refactor_inputs/*.intron_cache.tsv
