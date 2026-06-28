#!/usr/bin/env python3

"""
Compute splice-site depth offsets for intron features from BAM files.

This utility creates an intron x sample matrix of local alignment-depth
denominators. For each intron and sample, it computes:

    max(donor_site_window_depth, acceptor_site_window_depth)

Depth is currently unstranded. The code keeps strand in the site identifiers and
has a single compatibility function so strand-aware library modes can be added
without changing the output format.
"""

import argparse
import logging
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

from cluster_introns import parse_intron_id


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


EXCLUDE_COLS = {
    "intron_info",
    "donor_cluster",
    "acceptor_cluster",
    "gene_name",
    "intron_status",
    "overlapping_genes",
}


def parse_bam_args(bam_args):
    """
    Parse repeated --bam sample_id=path arguments.
    """
    sample_to_bam = {}
    for bam_arg in bam_args or []:
        if "=" not in bam_arg:
            raise ValueError("--bam values must use sample_id=/path/to/sample.bam")
        sample_id, bam_path = bam_arg.split("=", 1)
        sample_id = sample_id.strip()
        bam_path = bam_path.strip()
        if not sample_id or not bam_path:
            raise ValueError("--bam values must use sample_id=/path/to/sample.bam")
        sample_to_bam[sample_id] = bam_path
    return sample_to_bam


def load_bam_list(bam_list_file):
    """
    Load a TSV with columns sample_id and bam.
    """
    if not bam_list_file:
        return {}

    bam_df = pd.read_csv(bam_list_file, sep="\t", comment="#")
    required_cols = {"sample_id", "bam"}
    missing = required_cols - set(bam_df.columns)
    if missing:
        raise ValueError(
            f"BAM list must contain columns: {', '.join(sorted(required_cols))}; "
            f"missing: {', '.join(sorted(missing))}"
        )

    sample_to_bam = {}
    for _, row in bam_df.iterrows():
        sample_id = str(row["sample_id"]).strip()
        bam_path = str(row["bam"]).strip()
        if sample_id and bam_path and sample_id.lower() != "nan" and bam_path.lower() != "nan":
            sample_to_bam[sample_id] = bam_path
    return sample_to_bam


def get_sample_columns(df):
    return [col for col in df.columns if col not in EXCLUDE_COLS]


def load_introns(matrix_file):
    logger.info(f"Loading intron matrix from {matrix_file}")
    df = pd.read_csv(matrix_file, sep="\t", index_col=0)
    sample_cols = get_sample_columns(df)
    logger.info(f"Loaded {len(df)} introns and {len(sample_cols)} matrix sample columns")

    intron_info = {}
    intron_sites = {}
    sites_by_chrom = defaultdict(set)

    for intron_id in df.index:
        info = parse_intron_id(intron_id)
        intron_info[intron_id] = info

        donor_site = (info["chr"], int(info["donor"]), info["strand"])
        acceptor_site = (info["chr"], int(info["acceptor"]), info["strand"])
        intron_sites[intron_id] = (donor_site, acceptor_site)
        sites_by_chrom[donor_site[0]].add(donor_site)
        sites_by_chrom[acceptor_site[0]].add(acceptor_site)

    total_sites = sum(len(sites) for sites in sites_by_chrom.values())
    logger.info(f"Identified {total_sites} unique splice sites across {len(sites_by_chrom)} contigs")

    return df.index.tolist(), sample_cols, intron_sites, sites_by_chrom


def make_window(pos_1based, radius):
    """
    Convert a 1-based splice-site position to a 0-based half-open window.
    """
    pos0 = int(pos_1based) - 1
    start = max(0, pos0 - radius)
    end = pos0 + radius + 1
    return start, end


def merge_intervals(intervals):
    if not intervals:
        return []

    intervals = sorted(intervals)
    merged = [list(intervals[0])]

    for start, end in intervals[1:]:
        last = merged[-1]
        if start <= last[1]:
            if end > last[1]:
                last[1] = end
        else:
            merged.append([start, end])

    return [(start, end) for start, end in merged]


def read_is_compatible_with_strand(read, site_strand, strandedness):
    """
    Return whether a read should contribute to a site on site_strand.

    Only unstranded mode is implemented for now. This function is intentionally
    isolated so future fr-firststrand/fr-secondstrand/direct-RNA modes can be
    added in one place.
    """
    if strandedness == "unstranded":
        return True

    raise NotImplementedError(f"Strandedness mode is not implemented yet: {strandedness}")


def make_read_callback(min_mapq, strandedness):
    def read_callback(read):
        if read.is_unmapped or read.is_secondary:
            return False
        if read.mapping_quality < min_mapq:
            return False
        return read_is_compatible_with_strand(read, None, strandedness)

    return read_callback


def compute_site_depths_for_bam(bam_path, sites_by_chrom, window_radius, min_mapq, strandedness):
    """
    Compute max depth in each splice-site window for one BAM.
    """
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"BAM not found: {bam_path}")

    logger.info(f"Opening BAM: {bam_path}")
    bam = pysam.AlignmentFile(bam_path, "rb")
    contigs = set(bam.references)
    read_callback = make_read_callback(min_mapq, strandedness)
    site_depths = {}

    for chrom in sorted(sites_by_chrom):
        sites = sorted(sites_by_chrom[chrom], key=lambda site: (site[1], site[2]))

        if chrom not in contigs:
            logger.warning(f"Contig {chrom} is not present in {bam_path}; assigning zero depth")
            for site in sites:
                site_depths[site] = 0
            continue

        site_windows = {
            site: make_window(site[1], window_radius)
            for site in sites
        }
        merged_intervals = merge_intervals(site_windows.values())
        logger.info(
            f"  {chrom}: {len(sites)} sites represented by {len(merged_intervals)} merged windows"
        )

        chrom_depths = {}
        site_idx = 0

        for start, end in merged_intervals:
            coverage = bam.count_coverage(
                chrom,
                start=start,
                stop=end,
                quality_threshold=0,
                read_callback=read_callback,
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
                    value = int(depth[local_start:local_end].max())
                    previous = chrom_depths.get(site, 0)
                    if value > previous:
                        chrom_depths[site] = value
                scan_idx += 1

        for site in sites:
            site_depths[site] = chrom_depths.get(site, 0)

    bam.close()
    return site_depths


def compute_intron_offsets(intron_ids, intron_sites, sample_to_bam, window_radius, min_mapq, strandedness):
    offsets = pd.DataFrame(index=intron_ids, columns=list(sample_to_bam), dtype=np.int64)

    sites_by_chrom = defaultdict(set)
    for donor_site, acceptor_site in intron_sites.values():
        sites_by_chrom[donor_site[0]].add(donor_site)
        sites_by_chrom[acceptor_site[0]].add(acceptor_site)

    for sample_id, bam_path in sample_to_bam.items():
        logger.info(f"Computing site-depth offsets for sample {sample_id}")
        site_depths = compute_site_depths_for_bam(
            bam_path=bam_path,
            sites_by_chrom=sites_by_chrom,
            window_radius=window_radius,
            min_mapq=min_mapq,
            strandedness=strandedness,
        )

        values = []
        for intron_id in intron_ids:
            donor_site, acceptor_site = intron_sites[intron_id]
            values.append(max(site_depths.get(donor_site, 0), site_depths.get(acceptor_site, 0)))
        offsets[sample_id] = values

        nonzero = offsets[sample_id][offsets[sample_id] > 0]
        if len(nonzero) > 0:
            logger.info(
                f"  {sample_id}: nonzero offsets for {len(nonzero)}/{len(offsets)} introns; "
                f"median={nonzero.median():.0f}, max={nonzero.max():.0f}"
            )
        else:
            logger.warning(f"  {sample_id}: all site-depth offsets are zero")

    return offsets


def main():
    parser = argparse.ArgumentParser(
        description="Compute unstranded splice-site depth offsets from BAM files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--matrix", required=True, help="Input intron matrix with intron IDs as rows")
    parser.add_argument("--output", required=True, help="Output intron x sample site-depth offset TSV")
    parser.add_argument(
        "--bam_list",
        default=None,
        help="TSV with sample_id and bam columns",
    )
    parser.add_argument(
        "--bam",
        action="append",
        default=[],
        help="Sample/BAM mapping as sample_id=/path/to/sample.bam. Can be repeated.",
    )
    parser.add_argument(
        "--window_radius",
        type=int,
        default=10,
        help="Bases on each side of the splice-site coordinate to include",
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=60,
        help="Minimum read mapping quality",
    )
    parser.add_argument(
        "--strandedness",
        choices=["unstranded"],
        default="unstranded",
        help="Library strandedness mode. Only unstranded is currently implemented.",
    )
    parser.add_argument(
        "--allow_extra_bams",
        action="store_true",
        help="Allow BAM samples not present as matrix columns; they will be ignored",
    )

    args = parser.parse_args()

    intron_ids, matrix_samples, intron_sites, _ = load_introns(args.matrix)

    sample_to_bam = load_bam_list(args.bam_list)
    sample_to_bam.update(parse_bam_args(args.bam))

    if not sample_to_bam:
        parser.error("Provide at least one BAM using --bam_list or --bam")

    missing_bams = [sample for sample in matrix_samples if sample not in sample_to_bam]
    if missing_bams:
        raise ValueError(
            "Missing BAM paths for matrix samples: " + ", ".join(missing_bams[:10]) +
            ("..." if len(missing_bams) > 10 else "")
        )

    extra_bams = [sample for sample in sample_to_bam if sample not in matrix_samples]
    if extra_bams and not args.allow_extra_bams:
        raise ValueError(
            "BAM samples not found in matrix columns: " + ", ".join(extra_bams[:10]) +
            ". Use --allow_extra_bams to ignore them."
        )

    sample_to_bam = {sample: sample_to_bam[sample] for sample in matrix_samples}

    offsets = compute_intron_offsets(
        intron_ids=intron_ids,
        intron_sites=intron_sites,
        sample_to_bam=sample_to_bam,
        window_radius=args.window_radius,
        min_mapq=args.min_mapq,
        strandedness=args.strandedness,
    )

    logger.info(f"Writing site-depth offsets to {args.output}")
    offsets.to_csv(args.output, sep="\t", na_rep="NA")
    logger.info("Site-depth offset computation complete")


if __name__ == "__main__":
    main()
