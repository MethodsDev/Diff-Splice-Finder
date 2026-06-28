#!/usr/bin/env python3

"""
Shared splice-site depth offset utilities.

Depth is computed over aligned reference blocks, not full template spans. For
paired-end reads, both mates are collapsed to a per-fragment union within each
depth interval so overlapping mates contribute at most one count per position.
"""

from collections import Counter, defaultdict
import logging
import os

import numpy as np
import pysam


logger = logging.getLogger(__name__)

STRAND_MODES = ("unstranded", "F", "R", "FR", "RF")


def make_window(pos_1based, radius):
    pos0 = int(pos_1based) - 1
    return max(0, pos0 - radius), pos0 + radius + 1


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        last = merged[-1]
        if start <= last[1]:
            last[1] = max(last[1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def read_passes_depth_filters(read, min_mapq):
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    if read.is_duplicate or read.is_qcfail:
        return False
    return read.mapping_quality >= min_mapq


def single_end_transcript_strand(read, strand_mode):
    if strand_mode == "F":
        return "-" if read.is_reverse else "+"
    if strand_mode == "R":
        return "+" if read.is_reverse else "-"
    raise ValueError(f"Unsupported single-end strand mode: {strand_mode}")


def paired_end_transcript_strand(read, strand_mode):
    if not read.is_read1 and not read.is_read2:
        return None

    if strand_mode == "FR":
        if read.is_read1:
            return "-" if read.is_reverse else "+"
        return "+" if read.is_reverse else "-"

    if strand_mode == "RF":
        if read.is_read1:
            return "+" if read.is_reverse else "-"
        return "-" if read.is_reverse else "+"

    raise ValueError(f"Unsupported paired-end strand mode: {strand_mode}")


def read_transcript_strand(read, strand_mode):
    if strand_mode == "unstranded":
        return None
    if strand_mode in {"F", "R"}:
        if read.is_paired:
            return None
        return single_end_transcript_strand(read, strand_mode)
    if strand_mode in {"FR", "RF"}:
        if not read.is_paired or not read.is_proper_pair or read.mate_is_unmapped:
            return None
        if read.next_reference_id != read.reference_id:
            return None
        return paired_end_transcript_strand(read, strand_mode)
    raise ValueError(f"Unsupported strand mode: {strand_mode}")


def read_is_compatible_with_site(read, site_strand, strand_mode):
    if strand_mode == "unstranded":
        return True
    if site_strand not in {"+", "-"}:
        return False
    return read_transcript_strand(read, strand_mode) == site_strand


def clipped_read_blocks(read, interval_start, interval_end):
    blocks = []
    for block_start, block_end in read.get_blocks():
        start = max(block_start, interval_start)
        end = min(block_end, interval_end)
        if end > start:
            blocks.append((start, end))
    return blocks


def fragment_key(read, strand_mode):
    if read.is_paired and strand_mode in {"unstranded", "FR", "RF"}:
        return ("pair", read.query_name)
    return ("read", read.query_name, read.flag, read.reference_start)


def add_blocks_to_group(groups, read, interval_start, interval_end, strand_mode):
    blocks = clipped_read_blocks(read, interval_start, interval_end)
    if not blocks:
        return

    if strand_mode in {"F", "R"}:
        if read.is_paired:
            return
        transcript_strand = read_transcript_strand(read, strand_mode)
        if transcript_strand is None:
            return
    elif strand_mode in {"FR", "RF"}:
        transcript_strand = read_transcript_strand(read, strand_mode)
        if transcript_strand is None:
            return
    else:
        transcript_strand = None

    key = fragment_key(read, strand_mode)
    group = groups[key]
    if transcript_strand is not None:
        if group["strand"] is None:
            group["strand"] = transcript_strand
        elif group["strand"] != transcript_strand:
            group["ambiguous"] = True
            return

    group["blocks"].extend(blocks)


def merged_blocks(blocks):
    return merge_intervals(blocks)


def interval_depth_by_strand(bam, chrom, start, end, strand_mode, min_mapq):
    """
    Return depth arrays for an interval keyed by site strand.
    """
    length = end - start
    if strand_mode == "unstranded":
        depths = {None: np.zeros(length, dtype=np.int32)}
    else:
        depths = {
            "+": np.zeros(length, dtype=np.int32),
            "-": np.zeros(length, dtype=np.int32),
        }

    groups = defaultdict(lambda: {"blocks": [], "strand": None, "ambiguous": False})
    skipped = Counter()

    for read in bam.fetch(chrom, start, end):
        if not read_passes_depth_filters(read, min_mapq):
            skipped["filtered"] += 1
            continue
        if strand_mode in {"F", "R"} and read.is_paired:
            skipped["paired_in_single_end_mode"] += 1
            continue
        if strand_mode in {"FR", "RF"}:
            if not read.is_paired:
                skipped["unpaired_in_paired_end_mode"] += 1
                continue
            if not read.is_proper_pair or read.mate_is_unmapped or read.next_reference_id != read.reference_id:
                skipped["unusable_pair"] += 1
                continue
        add_blocks_to_group(groups, read, start, end, strand_mode)

    for group in groups.values():
        if group["ambiguous"] or not group["blocks"]:
            skipped["ambiguous_fragment"] += 1
            continue
        depth_key = group["strand"] if strand_mode != "unstranded" else None
        if depth_key not in depths:
            skipped["unknown_strand"] += 1
            continue
        for block_start, block_end in merged_blocks(group["blocks"]):
            depths[depth_key][block_start - start : block_end - start] += 1

    return depths, skipped


def compute_site_depth_offsets(
    bam_filename,
    intron_sites_by_chrom,
    window_radius=10,
    min_mapping_quality=60,
    strand_mode="unstranded",
):
    """
    Compute max splice-boundary window depth for each intron.

    intron_sites_by_chrom maps chrom -> intron_id -> [(site_pos, site_strand), ...].
    """
    if strand_mode not in STRAND_MODES:
        raise ValueError(f"Unsupported strand mode: {strand_mode}")
    if not os.path.exists(bam_filename):
        raise FileNotFoundError(f"BAM not found: {bam_filename}")

    logger.info(f"-computing splice-site depth offsets using strand mode: {strand_mode}")

    bam = pysam.AlignmentFile(bam_filename, "rb")
    contigs = set(bam.references)
    offsets = {}
    total_skipped = Counter()

    for chrom, intron_to_sites in sorted(intron_sites_by_chrom.items()):
        if not intron_to_sites:
            continue
        if chrom not in contigs:
            logger.warning(f"Contig {chrom} is not present in BAM; site-depth offsets set to zero")
            for intron in intron_to_sites:
                offsets[intron] = 0
            continue

        sites = sorted({site for site_list in intron_to_sites.values() for site in site_list})
        site_windows = {site: make_window(site[0], window_radius) for site in sites}
        merged_depth_intervals = merge_intervals(site_windows.values())

        logger.info(
            f"  {chrom}: {len(intron_to_sites)} introns, {len(sites)} sites, "
            f"{len(merged_depth_intervals)} merged depth windows"
        )

        site_depths = {site: 0 for site in sites}
        site_idx = 0

        for start, end in merged_depth_intervals:
            interval_depths, skipped = interval_depth_by_strand(
                bam,
                chrom,
                start,
                end,
                strand_mode,
                min_mapping_quality,
            )
            total_skipped.update(skipped)

            while site_idx < len(sites) and site_windows[sites[site_idx]][1] <= start:
                site_idx += 1

            scan_idx = site_idx
            while scan_idx < len(sites):
                site = sites[scan_idx]
                win_start, win_end = site_windows[site]
                if win_start >= end:
                    break
                local_start = max(win_start, start) - start
                local_end = min(win_end, end) - start
                if local_end > local_start:
                    depth_key = None if strand_mode == "unstranded" else site[1]
                    depth_array = interval_depths.get(depth_key)
                    if depth_array is not None:
                        value = int(depth_array[local_start:local_end].max())
                        if value > site_depths[site]:
                            site_depths[site] = value
                scan_idx += 1

        for intron, intron_sites in intron_to_sites.items():
            offsets[intron] = max(site_depths.get(site, 0) for site in intron_sites)

    bam.close()

    if total_skipped:
        logger.info(
            "Skipped depth alignments/fragments: "
            + ", ".join(f"{key}={value}" for key, value in sorted(total_skipped.items()))
        )

    return offsets


def write_strand_partitioned_bams(bam_filename, output_prefix, strand_mode, min_mapping_quality=60):
    """
    Write transcript-plus and transcript-minus BAMs for stranded modes.
    """
    if strand_mode == "unstranded":
        return []
    if strand_mode not in STRAND_MODES:
        raise ValueError(f"Unsupported strand mode: {strand_mode}")

    plus_bam = f"{output_prefix}.transcript_plus.bam"
    minus_bam = f"{output_prefix}.transcript_minus.bam"
    skipped = Counter()
    written = Counter()

    with pysam.AlignmentFile(bam_filename, "rb") as bam:
        with pysam.AlignmentFile(plus_bam, "wb", template=bam) as plus_out, pysam.AlignmentFile(
            minus_bam, "wb", template=bam
        ) as minus_out:
            for read in bam.fetch(until_eof=True):
                if not read_passes_depth_filters(read, min_mapping_quality):
                    skipped["filtered"] += 1
                    continue
                transcript_strand = read_transcript_strand(read, strand_mode)
                if transcript_strand == "+":
                    plus_out.write(read)
                    written["+"] += 1
                elif transcript_strand == "-":
                    minus_out.write(read)
                    written["-"] += 1
                else:
                    skipped["unassigned"] += 1

    pysam.index(plus_bam)
    pysam.index(minus_bam)
    logger.info(
        "Wrote strand-partitioned BAMs: "
        f"{plus_bam} ({written['+']} alignments), "
        f"{minus_bam} ({written['-']} alignments)"
    )
    if skipped:
        logger.info(
            "Skipped strand-BAM alignments: "
            + ", ".join(f"{key}={value}" for key, value in sorted(skipped.items()))
        )
    return [plus_bam, minus_bam]
