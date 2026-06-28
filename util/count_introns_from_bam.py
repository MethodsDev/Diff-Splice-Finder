#!/usr/bin/env python3

import sys, os, re
import pysam
from collections import defaultdict
import logging
import argparse
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


MIN_MAPPING_QUALITY = 60
DEFAULT_SITE_DEPTH_WINDOW_RADIUS = 10

OK_SPLICES = (
    "GT--AG",
    "GC--AG",
    "AT--AC",  # forward strand
    "CT--AC",
    "CT--GC",
    "GT--AT",  # reverse strand
)


def main():

    parser = argparse.ArgumentParser(
        description="annotate introns",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--genome_fa",
        type=str,
        required=True,
        help="reference genome fasta file",
    )

    parser.add_argument("--bam", type=str, required=True, help="bam file")

    parser.add_argument(
        "--target_introns",
        type=str,
        default=None,
        help=(
            "Optional intron list or count matrix. When provided, report every "
            "target intron with its count and splice-site depth offset, even if "
            "the junction count is zero in this BAM."
        ),
    )

    parser.add_argument(
        "--site_depth_window_radius",
        type=int,
        default=DEFAULT_SITE_DEPTH_WINDOW_RADIUS,
        help="Bases on each side of each intron boundary to include for site-depth offsets",
    )

    parser.add_argument(
        "--min_mapping_quality",
        type=int,
        default=MIN_MAPPING_QUALITY,
        help="Minimum read mapping quality",
    )

    args = parser.parse_args()

    genome_fasta_filename = args.genome_fa
    bam_file = args.bam

    if not os.path.exists(genome_fasta_filename):
        exit(f"Error, genome fasta file not found: {genome_fasta_filename}")
    if not os.path.exists(bam_file):
        exit(f"Error, bam file not found: {bam_file}")

    fasta_reader = pysam.FastaFile(genome_fasta_filename)

    intron_counter = defaultdict(lambda: defaultdict(int))

    evaluate_introns_from_bam_file(
        bam_file,
        intron_counter,
        min_mapping_quality=args.min_mapping_quality,
    )

    target_introns = load_target_introns(args.target_introns) if args.target_introns else set()
    if target_introns:
        logger.info(f"Loaded {len(target_introns)} target introns")

    report_introns_by_chrom = defaultdict(set)
    for chrom, chrom_icounter in intron_counter.items():
        report_introns_by_chrom[chrom].update(chrom_icounter.keys())
    for intron in target_introns:
        chrom, _, _ = parse_intron_coord(intron)
        report_introns_by_chrom[chrom].add(intron)

    site_depth_offsets = compute_site_depth_offsets(
        bam_file,
        report_introns_by_chrom,
        window_radius=args.site_depth_window_radius,
        min_mapping_quality=args.min_mapping_quality,
    )

    ######
    logger.info("reporting intron annotations")

    # print header
    print(
        "\t".join(
            [
                "intron",
                "splice_pair",
                "splice_flag",
                "count",
                "site_depth_offset",
            ]
        )
    )

    for chrom, chrom_introns in report_introns_by_chrom.items():

        if "_" in chrom:  # only main chromosomes
            continue

        chrom_seq = fasta_reader.fetch(chrom)

        for intron in sorted(chrom_introns):

            chromval, lend, rend = parse_intron_coord(intron)
            count = intron_counter[chrom][intron]

            intron_key = intron

            left_dinuc = chrom_seq[lend - 1 : lend + 1]
            right_dinuc = chrom_seq[rend - 1 - 1 : rend]

            splice_tok = f"{left_dinuc}--{right_dinuc}"

            splice_flag = "OK" if splice_tok in OK_SPLICES else "NON"

            print(
                "\t".join(
                    [
                        intron_key,
                        splice_tok,
                        splice_flag,
                        str(count),
                        str(site_depth_offsets[intron]),
                    ]
                )
            )


def parse_intron_coord(intron):
    """
    Parse either raw intron coordinates or annotated intron IDs.
    """
    intron = intron.split("^", 1)[0]
    chromval, coords_val = intron.split(":")
    lend, rend = coords_val.split("-")
    return chromval, int(lend), int(rend)


def load_target_introns(target_introns_file):
    """
    Load target introns from a one-column list or a matrix whose first column is
    intron IDs.
    """
    target_introns = set()
    with open(target_introns_file, "rt") as fh:
        header = fh.readline()
        if not header:
            return target_introns

        first_field = header.rstrip("\n").split("\t", 1)[0]
        if first_field and first_field not in {"intron", "intron_id"} and ":" in first_field:
            target_introns.add(first_field.split("^", 1)[0])

        for line in fh:
            if not line.strip():
                continue
            intron = line.rstrip("\n").split("\t", 1)[0]
            if intron in {"intron", "intron_id"}:
                continue
            target_introns.add(intron.split("^", 1)[0])

    return target_introns


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


def read_passes_filters(read, min_mapping_quality):
    if read.is_unmapped or read.is_secondary:
        return False
    return read.mapping_quality >= min_mapping_quality


def compute_site_depth_offsets(
    bam_filename,
    report_introns_by_chrom,
    window_radius=DEFAULT_SITE_DEPTH_WINDOW_RADIUS,
    min_mapping_quality=MIN_MAPPING_QUALITY,
):
    """
    Compute max splice-boundary window depth for every reported intron.
    """
    logger.info("-computing splice-site depth offsets")

    bam_reader = pysam.AlignmentFile(bam_filename, "rb")
    contigs = set(bam_reader.references)
    offsets = {}

    for chrom, introns in sorted(report_introns_by_chrom.items()):
        if not introns:
            continue
        if chrom not in contigs:
            logger.warning(f"Contig {chrom} is not present in BAM; site-depth offsets set to zero")
            for intron in introns:
                offsets[intron] = 0
            continue

        sites = sorted({pos for intron in introns for pos in parse_intron_coord(intron)[1:]})
        site_windows = {site: make_window(site, window_radius) for site in sites}
        merged_intervals = merge_intervals(site_windows.values())

        logger.info(
            f"  {chrom}: {len(introns)} introns, {len(sites)} sites, "
            f"{len(merged_intervals)} merged depth windows"
        )

        site_depths = {site: 0 for site in sites}
        site_idx = 0

        for start, end in merged_intervals:
            coverage = bam_reader.count_coverage(
                chrom,
                start=start,
                stop=end,
                quality_threshold=0,
                read_callback=lambda read: read_passes_filters(read, min_mapping_quality),
            )
            depth = np.sum(np.vstack(coverage), axis=0)

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
                    site_depths[site] = max(site_depths[site], int(depth[local_start:local_end].max()))
                scan_idx += 1

        for intron in introns:
            _, lend, rend = parse_intron_coord(intron)
            offsets[intron] = max(site_depths.get(lend, 0), site_depths.get(rend, 0))

    bam_reader.close()
    return offsets


def evaluate_introns_from_bam_file(bam_filename, intron_counter, min_mapping_quality=MIN_MAPPING_QUALITY):

    logger.info("-searching bam file")

    bam_reader = pysam.AlignmentFile(bam_filename, "rb")
    for read in bam_reader.fetch():
        if read.is_unmapped:
            continue

        if read.mapping_quality < min_mapping_quality:
            continue

        if read.is_secondary:
            continue

        chrom = bam_reader.get_reference_name(read.reference_id)

        introns = bam_reader.find_introns([read])

        for coordpair in introns.keys():
            lend, rend = coordpair
            lend += 1

            intron_key = f"{chrom}:{lend}-{rend}"

            intron_counter[chrom][intron_key] += 1

    return


if __name__ == "__main__":
    main()
