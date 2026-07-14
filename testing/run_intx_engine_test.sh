#!/usr/bin/env bash

# Exercise interaction-mode engine selection from committed fixture matrices.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
INPUT_DIR="$HERE/mode_refactor_inputs"
OUT_DIR="$INPUT_DIR/intx_engine_runs"

run_engine() {
  local engine="$1"
  local out="$OUT_DIR/$engine"
  rm -rf "$out"
  mkdir -p "$out"

  "$ROOT/DSF.py" \
    --matrix "$INPUT_DIR/intron_counts.matrix" \
    --offset_matrix "$INPUT_DIR/intron_counts.max_splice_plus_retained_depth.matrix" \
    --samples "$INPUT_DIR/sample_metadata.tsv" \
    --output_dir "$out" \
    --contrast A,B \
    --offset_mode splice_plus_retained \
    --stat_mode interaction \
    --intx_engine "$engine" \
    --min_intron_count 1 \
    --min_intron_samples 1 \
    --min_offset_depth 1 \
    --min_offset_samples 1 \
    --min_delta_psi 0 \
    --fdr_threshold 1 \
    --min_logFC 0 \
    --force_rerun

  test -s "$out/edgeR_results.all.tsv"
  test -s "$out/workdir/edgeR_results.intron_results.tsv"
  python3 - <<PY
import pandas as pd
engine = "$engine"
df = pd.read_csv("$out/edgeR_results.all.tsv", sep="\t")
required = {"intron_id", "logFC", "PValue", "FDR", "stat_mode", "intx_engine", "delta_PSI"}
missing = required - set(df.columns)
if missing:
    raise SystemExit(f"{engine}: missing columns {sorted(missing)}")
if set(df["stat_mode"].astype(str)) != {"interaction"}:
    raise SystemExit(f"{engine}: unexpected stat_mode values {set(df['stat_mode'].astype(str))}")
if set(df["intx_engine"].astype(str)) != {engine}:
    raise SystemExit(f"{engine}: unexpected intx_engine values {set(df['intx_engine'].astype(str))}")
print(f"[intx-engine] {engine}: {len(df)} introns")
PY
}

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

run_engine edgeR

if Rscript -e 'quit(status = !requireNamespace("DEXSeq", quietly=TRUE))' >/dev/null 2>&1; then
  run_engine DEXSeq
else
  echo "[intx-engine] DEXSeq not installed; skipping DEXSeq engine fixture"
fi

echo "[intx-engine] engine tests passed"
