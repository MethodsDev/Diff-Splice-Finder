#!/usr/bin/env python3

"""
Depth-window calculations backed by samtools depth.

The intron counter uses these helpers to report reusable read-depth columns in
per-sample .introns files. The depth values are pileup summaries over targeted
windows:

  left/right adjacent depth: exon-side windows outside the intron
  left/right retained depth: intron-side windows inside the intron

All depth values are max depth over the requested window. Splice-depth values
are derived separately from intron junction counts.
"""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile
from collections import defaultdict

from site_depth import STRAND_MODES, write_strand_partitioned_bams

logger = logging.getLogger(__name__)


def parse_intron_coord(intron):
    intron = str(intron).split("^", 1)[0]
    chrom, coords = intron.split(":", 1)
    lend, rend = coords.split("-", 1)
    return chrom, int(lend), int(rend)


def _window(start, end):
    start = max(0, int(start))
    end = max(start, int(end))
    return start, end


def intron_depth_windows(intron, radius):
    """Return BED-style half-open windows for an intron."""
    chrom, lend, rend = parse_intron_coord(intron)

    # Coordinates:
    #   lend = first intronic base (1-based) -> 0-based lend - 1
    #   rend = last intronic base (1-based)  -> 0-based rend - 1
    # Exon-adjacent windows are outside the intron. Retained windows are inside.
    left_adjacent = _window(lend - 1 - radius, lend - 1)
    right_adjacent = _window(rend, rend + radius)
    left_retained = _window(lend - 1, min(rend, lend - 1 + radius))
    right_retained = _window(max(lend - 1, rend - radius), rend)

    return {
        "left_adjacent_depth": (chrom, *left_adjacent),
        "right_adjacent_depth": (chrom, *right_adjacent),
        "left_retained_depth": (chrom, *left_retained),
        "right_retained_depth": (chrom, *right_retained),
    }


def _write_bed(path, windows):
    with open(path, "wt") as bed:
        for chrom, start, end in sorted(set(windows)):
            if end > start:
                bed.write(f"{chrom}\t{start}\t{end}\n")


def _run_samtools_depth(bam_file, windows, min_mapping_quality, workdir):
    """Return {(chrom, 1-based-pos): depth} for all requested windows."""
    if not windows:
        return {}

    bed_file = os.path.join(workdir, "depth_windows.bed")
    _write_bed(bed_file, windows)

    cmd = [
        "samtools",
        "depth",
        "-aa",
        "-s",
        "-Q",
        str(int(min_mapping_quality)),
        "-G",
        "SUPPLEMENTARY",
        "-b",
        bed_file,
        bam_file,
    ]
    logger.info("Running samtools depth for targeted depth windows")
    result = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE)

    depth_by_pos = {}
    for line in result.stdout.splitlines():
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        chrom = fields[0]
        pos1 = int(fields[1])
        depth = int(fields[2])
        depth_by_pos[(chrom, pos1)] = depth
    return depth_by_pos


def _max_depth(depth_by_pos, window):
    chrom, start, end = window
    if end <= start:
        return 0
    # samtools depth reports 1-based positions; BED [start, end) maps to
    # 1-based positions start + 1 through end.
    return max((depth_by_pos.get((chrom, pos1), 0) for pos1 in range(start + 1, end + 1)), default=0)


def _compute_unstranded_depths(bam_file, introns, radius, min_mapping_quality, workdir):
    intron_windows = {intron: intron_depth_windows(intron, radius) for intron in introns}
    all_windows = [window for windows in intron_windows.values() for window in windows.values()]
    depth_by_pos = _run_samtools_depth(bam_file, all_windows, min_mapping_quality, workdir)

    out = {}
    for intron, windows in intron_windows.items():
        values = {name: _max_depth(depth_by_pos, window) for name, window in windows.items()}
        values["max_adjacent_depth"] = max(values["left_adjacent_depth"], values["right_adjacent_depth"])
        out[intron] = values
    return out


def compute_depth_columns(
    bam_file,
    intron_to_strand,
    radius=10,
    min_mapping_quality=60,
    strand_mode="unstranded",
):
    """
    Compute adjacent and retained depth columns for target introns.

    intron_to_strand maps raw intron IDs (chr:lend-rend) to '+', '-', or '?'.
    In stranded modes, canonical '+'/'-' introns are evaluated against the
    corresponding transcript-strand BAM. Unknown-strand introns receive zero
    strand-specific depths.
    """
    if strand_mode not in STRAND_MODES:
        raise ValueError(f"Unsupported strand mode: {strand_mode}")
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM not found: {bam_file}")

    introns = sorted(intron_to_strand)
    if not introns:
        return {}

    if strand_mode == "unstranded":
        with tempfile.TemporaryDirectory(prefix="dsf_depth_") as tmpdir:
            return _compute_unstranded_depths(
                bam_file,
                introns,
                radius,
                min_mapping_quality,
                tmpdir,
            )

    by_strand = defaultdict(list)
    for intron, strand in intron_to_strand.items():
        if strand in {"+", "-"}:
            by_strand[strand].append(intron)

    out = {
        intron: {
            "left_adjacent_depth": 0,
            "right_adjacent_depth": 0,
            "max_adjacent_depth": 0,
            "left_retained_depth": 0,
            "right_retained_depth": 0,
        }
        for intron in introns
    }

    with tempfile.TemporaryDirectory(prefix="dsf_depth_stranded_") as tmpdir:
        prefix = os.path.join(tmpdir, "strand_partitioned")
        plus_bam, minus_bam = write_strand_partitioned_bams(
            bam_file,
            prefix,
            strand_mode,
            min_mapping_quality=min_mapping_quality,
        )
        strand_bams = {"+": plus_bam, "-": minus_bam}
        for strand, strand_introns in by_strand.items():
            if not strand_introns:
                continue
            strand_depths = _compute_unstranded_depths(
                strand_bams[strand],
                strand_introns,
                radius,
                min_mapping_quality,
                tmpdir,
            )
            out.update(strand_depths)

    return out
