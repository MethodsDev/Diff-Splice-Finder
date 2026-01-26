#!/usr/bin/env python3

"""
Cluster introns by shared donor or acceptor splice sites.

Given an intron count matrix, groups introns that share either:
- The same donor (5' splice site): chr:donor_position
- The same acceptor (3' splice site): chr:acceptor_position

This creates the denominators for compositional splicing analysis.
"""

import sys
import os
import argparse
import logging
import pandas as pd
from collections import defaultdict

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def parse_intron_id(intron_id):
    """
    Parse intron ID into components.
    
    Format: chr:start-end^splice_pair^flag
    Example: chr20:3864789-3865682^GT--AG^OK
    
    Returns:
        dict with keys: chr, start, end, strand, splice_pair, flag, donor, acceptor
    """
    try:
        coords, splice_pair, flag = intron_id.split("^")
        chr_pos, end = coords.split("-")
        chrom, start = chr_pos.split(":")
        
        start = int(start)
        end = int(end)
        
        # Determine strand from splice dinucleotides
        # Forward strand: GT--AG, GC--AG, AT--AC
        # Reverse strand: CT--AC, CT--GC, GT--AT
        forward_splices = {"GT--AG", "GC--AG", "AT--AC"}
        reverse_splices = {"CT--AC", "CT--GC", "GT--AT"}
        
        if splice_pair in forward_splices:
            strand = "+"
            donor = start
            acceptor = end
        elif splice_pair in reverse_splices:
            strand = "-"
            donor = end
            acceptor = start
        else:
            # Unknown/non-canonical - still assign but mark
            strand = "?"
            donor = start
            acceptor = end
        
        return {
            "chr": chrom,
            "start": start,
            "end": end,
            "strand": strand,
            "splice_pair": splice_pair,
            "flag": flag,
            "donor": donor,
            "acceptor": acceptor,
            "intron_id": intron_id,
        }
    except Exception as e:
        logger.error(f"Failed to parse intron ID: {intron_id}, error: {e}")
        raise


def cluster_by_donor(introns_df):
    """
    Cluster introns by shared donor (5' splice site).
    
    Args:
        introns_df: DataFrame with intron annotations and counts
        
    Returns:
        DataFrame with added 'donor_cluster' column
    """
    logger.info("Clustering introns by donor site...")
    
    donor_clusters = []
    for idx, row in introns_df.iterrows():
        intron_info = row["intron_info"]
        # Cluster ID: chr:donor_position:strand
        cluster_id = f"{intron_info['chr']}:{intron_info['donor']}:{intron_info['strand']}"
        donor_clusters.append(cluster_id)
    
    introns_df["donor_cluster"] = donor_clusters
    
    cluster_sizes = introns_df["donor_cluster"].value_counts()
    logger.info(f"Created {len(cluster_sizes)} donor clusters")
    logger.info(f"Cluster size distribution: min={cluster_sizes.min()}, "
                f"median={cluster_sizes.median():.1f}, max={cluster_sizes.max()}")
    
    return introns_df


def cluster_by_acceptor(introns_df):
    """
    Cluster introns by shared acceptor (3' splice site).
    
    Args:
        introns_df: DataFrame with intron annotations and counts
        
    Returns:
        DataFrame with added 'acceptor_cluster' column
    """
    logger.info("Clustering introns by acceptor site...")
    
    acceptor_clusters = []
    for idx, row in introns_df.iterrows():
        intron_info = row["intron_info"]
        # Cluster ID: chr:acceptor_position:strand
        cluster_id = f"{intron_info['chr']}:{intron_info['acceptor']}:{intron_info['strand']}"
        acceptor_clusters.append(cluster_id)
    
    introns_df["acceptor_cluster"] = acceptor_clusters
    
    cluster_sizes = introns_df["acceptor_cluster"].value_counts()
    logger.info(f"Created {len(cluster_sizes)} acceptor clusters")
    logger.info(f"Cluster size distribution: min={cluster_sizes.min()}, "
                f"median={cluster_sizes.median():.1f}, max={cluster_sizes.max()}")
    
    return introns_df


def load_and_annotate_introns(matrix_file):
    """
    Load intron count matrix and parse intron annotations.
    
    Args:
        matrix_file: Path to intron count matrix (tab-delimited)
        
    Returns:
        DataFrame with intron_info column containing parsed annotations
    """
    logger.info(f"Loading intron count matrix from {matrix_file}")
    
    df = pd.read_csv(matrix_file, sep="\t", index_col=0)
    logger.info(f"Loaded {len(df)} introns across {len(df.columns)} samples")
    
    # Parse intron IDs
    logger.info("Parsing intron annotations...")
    intron_infos = []
    for intron_id in df.index:
        intron_infos.append(parse_intron_id(intron_id))
    
    df["intron_info"] = intron_infos
    
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Cluster introns by shared donor or acceptor splice sites",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    parser.add_argument(
        "--matrix",
        type=str,
        required=True,
        help="Input intron count matrix (output from build_intron_count_matrix.py)",
    )
    
    parser.add_argument(
        "--output_donor",
        type=str,
        required=False,
        help="Output matrix with donor cluster annotations",
    )
    
    parser.add_argument(
        "--output_acceptor",
        type=str,
        required=False,
        help="Output matrix with acceptor cluster annotations",
    )
    
    parser.add_argument(
        "--cluster_type",
        type=str,
        choices=["donor", "acceptor", "both"],
        default="both",
        help="Type of clustering to perform",
    )
    
    args = parser.parse_args()
    
    # Load and annotate introns
    introns_df = load_and_annotate_introns(args.matrix)
    
    # Perform clustering
    if args.cluster_type in ["donor", "both"]:
        introns_df = cluster_by_donor(introns_df)
        if args.output_donor:
            logger.info(f"Writing donor-clustered matrix to {args.output_donor}")
            # Save with cluster annotation
            output_df = introns_df.drop(columns=["intron_info"])
            output_df.to_csv(args.output_donor, sep="\t")
    
    if args.cluster_type in ["acceptor", "both"]:
        introns_df = cluster_by_acceptor(introns_df)
        if args.output_acceptor:
            logger.info(f"Writing acceptor-clustered matrix to {args.output_acceptor}")
            # Save with cluster annotation
            output_df = introns_df.drop(columns=["intron_info"])
            output_df.to_csv(args.output_acceptor, sep="\t")
    
    logger.info("Clustering complete!")
    
    return introns_df


if __name__ == "__main__":
    main()
