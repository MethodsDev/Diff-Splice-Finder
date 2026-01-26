#!/usr/bin/env python3

import sys, os, re
import pysam
from collections import defaultdict
import logging
import argparse

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


MIN_MAPPING_QUALITY = 60

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

    args = parser.parse_args()

    genome_fasta_filename = args.genome_fa
    bam_file = args.bam

    if not os.path.exists(genome_fasta_filename):
        exit(f"Error, genome fasta file not found: {genome_fasta_filename}")
    if not os.path.exists(bam_file):
        exit(f"Error, bam file not found: {bam_file}")

    fasta_reader = pysam.FastaFile(genome_fasta_filename)

    intron_counter = defaultdict(lambda: defaultdict(int))

    evaluate_introns_from_bam_file(bam_file, intron_counter)

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
            ]
        )
    )

    for chrom, chrom_icounter in intron_counter.items():

        if "_" in chrom:  # only main chromosomes
            continue

        chrom_seq = fasta_reader.fetch(chrom)

        for intron, count in sorted(chrom_icounter.items()):

            chromval, coords_val = intron.split(":")
            lend, rend = coords_val.split("-")

            lend = int(lend)
            rend = int(rend)

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
                    ]
                )
            )


def evaluate_introns_from_bam_file(bam_filename, intron_counter):

    logger.info("-searching bam file")

    bam_reader = pysam.AlignmentFile(bam_filename, "rb")
    for read in bam_reader.fetch():
        if read.mapping_quality < MIN_MAPPING_QUALITY:
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
