#!/usr/bin/env python3

"""
Filter introns and clusters based on count thresholds and splice site quality.

Implements filtering strategy from the design doc:
1. Filter out non-canonical splice sites (keep only 'OK' flagged introns)
2. Filter introns with low total counts across samples
3. Filter introns with too few non-zero samples
4. Filter clusters with insufficient total reads
"""

import sys
import os
import argparse
import logging
import pandas as pd
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def filter_canonical_splices(introns_df, keep_noncanonical=False):
    """
    Filter introns to keep only canonical splice sites (OK flag).
    
    Args:
        introns_df: DataFrame with intron_info column
        keep_noncanonical: If True, keep non-canonical introns
        
    Returns:
        Filtered DataFrame
    """
    if keep_noncanonical:
        logger.info("Keeping all introns including non-canonical splice sites")
        return introns_df
    
    n_before = len(introns_df)
    
    canonical_mask = introns_df["intron_info"].apply(lambda x: x["flag"] == "OK")
    introns_df = introns_df[canonical_mask].copy()
    
    n_after = len(introns_df)
    n_removed = n_before - n_after
    
    logger.info(f"Canonical splice filter: kept {n_after} introns, removed {n_removed} "
                f"({100.0 * n_removed / n_before:.1f}%)")
    
    return introns_df


def filter_low_count_introns(introns_df, sample_cols, min_total_count=10, min_samples_nonzero=2):
    """
    Filter introns with insufficient read support.
    
    Args:
        introns_df: DataFrame with count data
        sample_cols: List of sample column names
        min_total_count: Minimum total reads across all samples
        min_samples_nonzero: Minimum number of samples with non-zero counts
        
    Returns:
        Filtered DataFrame
    """
    n_before = len(introns_df)
    
    # Calculate total counts and non-zero sample counts
    count_data = introns_df[sample_cols]
    total_counts = count_data.sum(axis=1)
    nonzero_samples = (count_data > 0).sum(axis=1)
    
    # Apply filters
    count_mask = total_counts >= min_total_count
    sample_mask = nonzero_samples >= min_samples_nonzero
    combined_mask = count_mask & sample_mask
    
    introns_df = introns_df[combined_mask].copy()
    
    n_after = len(introns_df)
    n_removed = n_before - n_after
    
    logger.info(f"Intron count filter (>={min_total_count} total, >={min_samples_nonzero} samples): "
                f"kept {n_after} introns, removed {n_removed} ({100.0 * n_removed / n_before:.1f}%)")
    
    # Log some stats
    if n_after > 0:
        logger.info(f"  Remaining introns - total count range: "
                    f"{introns_df[sample_cols].sum(axis=1).min():.0f} to "
                    f"{introns_df[sample_cols].sum(axis=1).max():.0f}")
    
    return introns_df


def filter_low_count_clusters(introns_df, sample_cols, cluster_col, 
                               min_cluster_count=20, min_cluster_samples=3):
    """
    Filter clusters with insufficient total read support.
    
    Args:
        introns_df: DataFrame with count data and cluster assignments
        sample_cols: List of sample column names
        cluster_col: Name of cluster column ('donor_cluster' or 'acceptor_cluster')
        min_cluster_count: Minimum total reads in cluster per sample
        min_cluster_samples: Minimum number of samples meeting count threshold
        
    Returns:
        Filtered DataFrame with low-count clusters removed
    """
    if cluster_col not in introns_df.columns:
        logger.warning(f"Cluster column '{cluster_col}' not found, skipping cluster filter")
        return introns_df
    
    n_before = len(introns_df)
    n_clusters_before = introns_df[cluster_col].nunique()
    
    # Calculate cluster totals per sample
    cluster_totals = introns_df.groupby(cluster_col)[sample_cols].sum()
    
    # Count how many samples meet the threshold for each cluster
    samples_passing = (cluster_totals >= min_cluster_count).sum(axis=1)
    
    # Keep clusters that pass in enough samples
    passing_clusters = samples_passing[samples_passing >= min_cluster_samples].index
    
    introns_df = introns_df[introns_df[cluster_col].isin(passing_clusters)].copy()
    
    n_after = len(introns_df)
    n_clusters_after = introns_df[cluster_col].nunique()
    n_removed = n_before - n_after
    n_clusters_removed = n_clusters_before - n_clusters_after
    
    logger.info(f"Cluster filter (>={min_cluster_count} reads in >={min_cluster_samples} samples): "
                f"kept {n_clusters_after} clusters ({n_after} introns), "
                f"removed {n_clusters_removed} clusters ({n_removed} introns, "
                f"{100.0 * n_removed / n_before:.1f}%)")
    
    return introns_df


def get_sample_columns(df):
    """
    Identify sample columns (exclude intron_info, cluster, gene_name, intron_status columns).
    
    Args:
        df: DataFrame with intron data
        
    Returns:
        List of sample column names
    """
    exclude_cols = {"intron_info", "donor_cluster", "acceptor_cluster", "gene_name", "intron_status", "overlapping_genes"}
    sample_cols = [col for col in df.columns if col not in exclude_cols]
    return sample_cols


def main():
    parser = argparse.ArgumentParser(
        description="Filter introns and clusters based on count thresholds",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    parser.add_argument(
        "--matrix",
        type=str,
        required=True,
        help="Input clustered intron matrix (from cluster_introns.py)",
    )
    
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output filtered matrix",
    )
    
    parser.add_argument(
        "--cluster_type",
        type=str,
        choices=["donor", "acceptor", "both"],
        default="both",
        help="Type of clustering used - 'both' requires thresholds for both donor and acceptor clusters",
    )
    
    parser.add_argument(
        "--min_intron_count",
        type=int,
        default=10,
        help="Minimum total count across all samples for an intron",
    )
    
    parser.add_argument(
        "--min_intron_samples",
        type=int,
        default=2,
        help="Minimum number of samples with non-zero counts for an intron",
    )
    
    parser.add_argument(
        "--min_cluster_count",
        type=int,
        default=20,
        help="Minimum total reads per sample for a cluster",
    )
    
    parser.add_argument(
        "--min_cluster_samples",
        type=int,
        default=3,
        help="Minimum number of samples meeting cluster count threshold",
    )
    
    parser.add_argument(
        "--keep_noncanonical",
        action="store_true",
        help="Keep non-canonical splice sites (by default only 'OK' sites are kept)",
    )
    
    args = parser.parse_args()
    
    # Load matrix
    logger.info(f"Loading clustered matrix from {args.matrix}")
    df = pd.read_csv(args.matrix, sep="\t", index_col=0)
    
    # Get sample columns
    sample_cols = get_sample_columns(df)
    logger.info(f"Found {len(sample_cols)} sample columns: {sample_cols}")
    
    # Parse intron info if not already present
    if "intron_info" not in df.columns:
        from cluster_introns import parse_intron_id
        logger.info("Parsing intron annotations...")
        df["intron_info"] = [parse_intron_id(idx) for idx in df.index]
    
    logger.info(f"Starting with {len(df)} introns")
    
    # Apply filters in order
    df = filter_canonical_splices(df, keep_noncanonical=args.keep_noncanonical)
    
    df = filter_low_count_introns(
        df, sample_cols,
        min_total_count=args.min_intron_count,
        min_samples_nonzero=args.min_intron_samples
    )
    
    # Apply cluster filters
    if args.cluster_type == "both":
        # Require both donor and acceptor clusters to pass thresholds
        logger.info("Filtering by BOTH donor and acceptor cluster thresholds...")
        df = filter_low_count_clusters(
            df, sample_cols, "donor_cluster",
            min_cluster_count=args.min_cluster_count,
            min_cluster_samples=args.min_cluster_samples
        )
        df = filter_low_count_clusters(
            df, sample_cols, "acceptor_cluster",
            min_cluster_count=args.min_cluster_count,
            min_cluster_samples=args.min_cluster_samples
        )
        logger.info(f"Final: {len(df)} introns in {df['donor_cluster'].nunique()} donor clusters "
                   f"and {df['acceptor_cluster'].nunique()} acceptor clusters")
    else:
        # Legacy: filter by single cluster type
        cluster_col = f"{args.cluster_type}_cluster"
        df = filter_low_count_clusters(
            df, sample_cols, cluster_col,
            min_cluster_count=args.min_cluster_count,
            min_cluster_samples=args.min_cluster_samples
        )
        logger.info(f"Final: {len(df)} introns in {df[cluster_col].nunique()} clusters")
    
    # Save filtered matrix
    logger.info(f"Writing filtered matrix to {args.output}")
    # Drop intron_info but keep gene_name and intron_status if present
    cols_to_drop = ["intron_info"]
    cols_to_drop = [c for c in cols_to_drop if c in df.columns]
    output_df = df.drop(columns=cols_to_drop)
    output_df.to_csv(args.output, sep="\t", na_rep='NA')
    
    logger.info("Filtering complete!")


if __name__ == "__main__":
    main()
