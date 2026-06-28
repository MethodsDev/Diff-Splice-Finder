#!/usr/bin/env python3

import sys, os, re
import csv
from collections import defaultdict
import argparse
import logging
import gzip


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="build sample-vs-intron count matrix from *.introns files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--output_matrix",
        type=str,
        required=True,
        help="output intron count matrix",
    )

    parser.add_argument(
        "--output_offset_matrix",
        type=str,
        default=None,
        help=(
            "output site-depth offset matrix. By default, writes next to "
            "--output_matrix when input intron files include site_depth_offset."
        ),
    )

    parser.add_argument(
        "--intron_files",
        type=str,
        nargs="+",
        required=True,
        help="space-delimited list of ${sample_name}.*.introns files to build into matrices (need at least 2)",
    )

    args = parser.parse_args()

    output_matrix = args.output_matrix
    output_offset_matrix = args.output_offset_matrix
    intron_files = args.intron_files

    if len(intron_files) < 2:
        print(
            "Error, need at least 2 introns files specified to build the matrices",
            file=sys.stderr,
        )

    intron_counts_matrix_data = defaultdict(lambda: defaultdict(int))
    intron_offsets_matrix_data = defaultdict(lambda: defaultdict(int))
    intron_ids = set()
    saw_offsets = False

    for intron_file in intron_files:
        sample_name = os.path.basename(intron_file)
        sample_name = sample_name.split(".")[0]

        logger.info("Parsing {}".format(intron_file))
        
        # Handle both plain and gzipped files
        if intron_file.endswith('.gz'):
            fh = gzip.open(intron_file, "rt")
        else:
            fh = open(intron_file, "rt")
        
        with fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:

                intron = row["intron"]
                splice_pair = row["splice_pair"]
                splice_flag = row["splice_flag"]

                count = int(row["count"])
                site_depth_offset = row.get("site_depth_offset")

                intron_token = "^".join([intron, splice_pair, splice_flag])

                intron_ids.add(intron_token)
                intron_counts_matrix_data[sample_name][intron_token] += count
                if site_depth_offset is not None and str(site_depth_offset).strip():
                    saw_offsets = True
                    intron_offsets_matrix_data[sample_name][intron_token] = max(
                        intron_offsets_matrix_data[sample_name][intron_token],
                        int(float(site_depth_offset)),
                    )

    logger.info("-writing output intron count matrix {}".format(output_matrix))

    with open(output_matrix, "wt") as counts_ofh:

        sample_names = list(intron_counts_matrix_data.keys())

        print("\t" + "\t".join(sample_names), file=counts_ofh)

        for intron_id in sorted(intron_ids):
            count_vals = [intron_id]

            for sample_name in sample_names:
                intron_id_sample_count = intron_counts_matrix_data[sample_name][
                    intron_id
                ]
                count_vals.append(str(intron_id_sample_count))

            print("\t".join(count_vals), file=counts_ofh)

    if saw_offsets:
        if not output_offset_matrix:
            if output_matrix.endswith(".matrix"):
                output_offset_matrix = output_matrix[:-7] + ".offsets.matrix"
            else:
                output_offset_matrix = output_matrix + ".offsets.tsv"

        logger.info("-writing output site-depth offset matrix {}".format(output_offset_matrix))
        missing_offsets = 0

        with open(output_offset_matrix, "wt") as offsets_ofh:
            print("\t" + "\t".join(sample_names), file=offsets_ofh)

            for intron_id in sorted(intron_ids):
                offset_vals = [intron_id]

                for sample_name in sample_names:
                    offset_val = intron_offsets_matrix_data[sample_name][intron_id]
                    if offset_val == 0 and intron_counts_matrix_data[sample_name][intron_id] == 0:
                        missing_offsets += 1
                    offset_vals.append(str(offset_val))

                print("\t".join(offset_vals), file=offsets_ofh)

        if missing_offsets:
            logger.warning(
                f"{missing_offsets} intron-sample combinations have zero site-depth offsets. "
                "For complete offsets, generate per-sample intron files with "
                "count_introns_from_bam.py --target_introns using a shared target intron list."
            )
    elif output_offset_matrix:
        logger.warning(
            "--output_offset_matrix was provided, but no input rows contained site_depth_offset; "
            "offset matrix was not written."
        )

    return


if __name__ == "__main__":
    main()
