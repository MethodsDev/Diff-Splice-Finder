#!/usr/bin/env python3

import sys, os, re
import pysam
from collections import defaultdict
import logging
import argparse

from site_depth import STRAND_MODES, compute_site_depth_offsets, write_strand_partitioned_bams

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

FORWARD_SPLICES = {"GT--AG", "GC--AG", "AT--AC"}
REVERSE_SPLICES = {"CT--AC", "CT--GC", "GT--AT"}


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

    parser.add_argument(
        "--site_depth_strand_mode",
        choices=STRAND_MODES,
        default="unstranded",
        help=(
            "Strand mode for site-depth offsets. F/R are single-end modes; "
            "FR/RF are paired-end modes describing read1/read2 orientations "
            "relative to the transcript."
        ),
    )

    parser.add_argument(
        "--strand_bam_prefix",
        type=str,
        default=None,
        help=(
            "Optional prefix for transcript_plus/transcript_minus BAM outputs. "
            "Only used when --site_depth_strand_mode is not unstranded."
        ),
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

    ######
    logger.info("reporting intron annotations")

    intron_annotations = {}
    intron_sites_by_chrom = defaultdict(dict)

    for chrom, chrom_introns in report_introns_by_chrom.items():

        if "_" in chrom:  # only main chromosomes
            continue

        chrom_seq = fasta_reader.fetch(chrom)

        for intron in sorted(chrom_introns):

            _, lend, rend = parse_intron_coord(intron)

            left_dinuc = chrom_seq[lend - 1 : lend + 1]
            right_dinuc = chrom_seq[rend - 1 - 1 : rend]

            splice_tok = f"{left_dinuc}--{right_dinuc}"

            splice_flag = "OK" if splice_tok in OK_SPLICES else "NON"
            strand = infer_splice_strand(splice_tok)
            intron_annotations[intron] = {
                "splice_pair": splice_tok,
                "splice_flag": splice_flag,
                "strand": strand,
            }
            intron_sites_by_chrom[chrom][intron] = [
                (lend, strand),
                (rend, strand),
            ]

    site_depth_offsets = compute_site_depth_offsets(
        bam_file,
        intron_sites_by_chrom,
        window_radius=args.site_depth_window_radius,
        min_mapping_quality=args.min_mapping_quality,
        strand_mode=args.site_depth_strand_mode,
    )

    if args.strand_bam_prefix and args.site_depth_strand_mode != "unstranded":
        write_strand_partitioned_bams(
            bam_file,
            args.strand_bam_prefix,
            args.site_depth_strand_mode,
            min_mapping_quality=args.min_mapping_quality,
        )

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

        for intron in sorted(chrom_introns):

            count = intron_counter[chrom][intron]
            annotation = intron_annotations[intron]

            print(
                "\t".join(
                    [
                        intron,
                        annotation["splice_pair"],
                        annotation["splice_flag"],
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


def infer_splice_strand(splice_pair):
    if splice_pair in FORWARD_SPLICES:
        return "+"
    if splice_pair in REVERSE_SPLICES:
        return "-"
    return "?"


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
