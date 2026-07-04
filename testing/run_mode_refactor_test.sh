#!/usr/bin/env bash

# Exercise the offset-mode execution paths from committed fixture matrices.
# The fixture uses a tiny synthetic chromosome, so the FASTA can be committed
# directly without gzip or GitHub large-file concerns.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
INPUT_DIR="$HERE/mode_refactor_inputs"
OUT_DIR="$INPUT_DIR/mode_runs"

FA="$INPUT_DIR/chrTiny.fa"

if [[ ! -s "$FA" ]]; then
  echo "Missing $FA" >&2
  exit 1
fi

if [[ ! -s "$FA.fai" ]]; then
  echo "[mode-refactor] indexing $FA"
  python3 - <<PY
import pysam
pysam.faidx("$FA")
PY
fi

run_mode() {
  local mode="$1"
  local offset_matrix="$2"
  shift 2

  local out="$OUT_DIR/$mode"
  rm -rf "$out"
  mkdir -p "$out"

  "$ROOT/run_diff_splice_analysis.py" \
    --matrix "$INPUT_DIR/intron_counts.matrix" \
    --offset_matrix "$offset_matrix" \
    --samples "$INPUT_DIR/sample_metadata.tsv" \
    --output_dir "$out" \
    --contrast A,B \
    --offset_mode "$mode" \
    "$@" \
    --min_intron_count 1 \
    --min_intron_samples 1 \
    --min_offset_depth 1 \
    --min_offset_samples 1 \
    --min_delta_psi 0 \
    --fdr_threshold 1 \
    --min_logFC 0

  test -s "$out/edgeR_results.all.tsv"
  test -s "$out/workdir/edgeR_results.intron_results_with_psi.tsv"
  python3 - <<PY
import pandas as pd
ann = pd.read_csv("$out/workdir/edgeR_input.annotations.tsv", sep="\t", index_col=0)
mode = "$mode"
observed = set(ann["offset_mode"].astype(str))
if observed != {mode}:
    raise SystemExit(f"{mode}: unexpected offset_mode values: {observed}")
print(f"[mode-refactor] {mode}: {len(ann)} introns")
PY
}

run_mode \
  splice_plus_retained \
  "$INPUT_DIR/intron_counts.max_splice_plus_retained_depth.matrix"

run_mode \
  splice_vs_rest \
  "$INPUT_DIR/intron_counts.max_splice_plus_retained_depth.matrix" \
  --gtf "$INPUT_DIR/chrTiny_fixture.gtf"

echo "[mode-refactor] all modes passed"
