#!/usr/bin/env python3
"""Exercise count_introns_from_bam.py against the checked-in test BAM."""

import csv
import os
import subprocess
import sys
import tempfile
from collections import Counter

import pysam


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
BAM_PATH = os.path.join(SCRIPT_DIR, "data", "alignments.b38.sorted.bam")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "bam_introns_test_output")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "alignments.b38.introns")
EXPECTED_INTRON = "chr7:55155947-55156532"
EXPECTED_COUNT = 342
EXPECTED_INTRONS = 49
EXPECTED_DEPTH_COLUMNS = {
    "left_adjacent_depth",
    "right_adjacent_depth",
    "max_adjacent_depth",
    "left_splice_depth",
    "right_splice_depth",
    "left_retained_depth",
    "right_retained_depth",
    "left_splice_plus_retained_depth",
    "right_splice_plus_retained_depth",
    "max_splice_plus_retained_depth",
    "splice_plus_retained_depth_source",
    "site_depth_offset",
}


def discover_introns_and_required_fasta_len():
    intron_counts = Counter()
    max_coord = 0
    with pysam.AlignmentFile(BAM_PATH, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.mapping_quality < 60 or read.is_secondary:
                continue
            chrom = bam.get_reference_name(read.reference_id)
            for start0, end in bam.find_introns([read]):
                start = start0 + 1
                intron_counts[f"{chrom}:{start}-{end}"] += 1
                max_coord = max(max_coord, start, end)

    if EXPECTED_INTRON not in intron_counts:
        raise AssertionError(f"Expected intron missing from BAM discovery: {EXPECTED_INTRON}")
    if intron_counts[EXPECTED_INTRON] != EXPECTED_COUNT:
        raise AssertionError(
            f"Expected {EXPECTED_INTRON} count {EXPECTED_COUNT}, "
            f"observed {intron_counts[EXPECTED_INTRON]}"
        )
    if len(intron_counts) != EXPECTED_INTRONS:
        raise AssertionError(f"Expected {EXPECTED_INTRONS} introns, observed {len(intron_counts)}")

    return intron_counts, max_coord + 100


def write_minimal_chr7_fasta(path, intron_counts, length):
    seq = bytearray(b"A" * length)
    for intron in intron_counts:
        chrom, coords = intron.split(":", 1)
        if chrom != "chr7":
            continue
        start, end = (int(value) for value in coords.split("-", 1))
        seq[start - 1 : start + 1] = b"GT"
        seq[end - 2 : end] = b"AG"

    with open(path, "wt") as fasta:
        fasta.write(">chr7\n")
        for i in range(0, len(seq), 80):
            fasta.write(seq[i : i + 80].decode("ascii") + "\n")
    pysam.faidx(path)


def run_count_introns(genome_fa):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    cmd = [
        sys.executable,
        os.path.join(REPO_ROOT, "util", "count_introns_from_bam.py"),
        "--genome_fa",
        genome_fa,
        "--bam",
        BAM_PATH,
        "--site_depth_window_radius",
        "10",
        "--min_mapping_quality",
        "60",
    ]
    with open(OUTPUT_FILE, "wt") as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "count_introns_from_bam.py failed with exit code "
            f"{result.returncode}\n{result.stderr}"
        )


def validate_output():
    with open(OUTPUT_FILE, "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    if len(rows) != EXPECTED_INTRONS:
        raise AssertionError(f"Expected {EXPECTED_INTRONS} output rows, observed {len(rows)}")
    missing_columns = EXPECTED_DEPTH_COLUMNS - set(rows[0])
    if missing_columns:
        raise AssertionError(f"Missing expected depth columns: {sorted(missing_columns)}")

    row_by_intron = {row["intron"]: row for row in rows}
    row = row_by_intron.get(EXPECTED_INTRON)
    if row is None:
        raise AssertionError(f"Expected output intron missing: {EXPECTED_INTRON}")
    if int(row["count"]) != EXPECTED_COUNT:
        raise AssertionError(
            f"Expected {EXPECTED_INTRON} count {EXPECTED_COUNT}, observed {row['count']}"
        )
    if int(float(row["site_depth_offset"])) <= 0:
        raise AssertionError(f"Expected positive site_depth_offset for {EXPECTED_INTRON}")
    if int(float(row["max_adjacent_depth"])) != int(float(row["site_depth_offset"])):
        raise AssertionError("Expected site_depth_offset compatibility alias to equal max_adjacent_depth")
    if int(float(row["max_splice_plus_retained_depth"])) < int(row["count"]):
        raise AssertionError("Expected splice-plus-retained denominator to be at least the focal count")


def main():
    if not os.path.exists(BAM_PATH):
        raise FileNotFoundError(f"Missing test BAM: {BAM_PATH}")

    intron_counts, fasta_len = discover_introns_and_required_fasta_len()
    with tempfile.TemporaryDirectory(prefix="dsf_bam_introns_test_") as tmpdir:
        genome_fa = os.path.join(tmpdir, "chr7.fa")
        write_minimal_chr7_fasta(genome_fa, intron_counts, fasta_len)
        run_count_introns(genome_fa)

    validate_output()
    print(f"BAM intron-count test passed: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
