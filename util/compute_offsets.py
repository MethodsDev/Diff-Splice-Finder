#!/usr/bin/env python3

"""
Compute cluster-total offsets for compositional splicing analysis.

For each sample, calculates the total reads supporting all introns within
each cluster. These totals are used as offsets in the edgeR GLM to normalize
for cluster-level expression differences, allowing us to test intron usage
proportions rather than absolute counts.

log(μ_i,s) = X_s * β_i + log(T_C,s)

where T_C,s is the cluster total offset for cluster C in sample s.
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


def compute_cluster_offsets(introns_df, sample_cols, cluster_col):
    """
    Compute cluster-total offsets for each sample.
    
    For each cluster and sample, sums all intron counts within that cluster.
    These totals become the offsets for the GLM.
    
    Args:
        introns_df: DataFrame with intron counts and cluster assignments
        sample_cols: List of sample column names
        cluster_col: Name of cluster column ('donor_cluster' or 'acceptor_cluster')
        
    Returns:
        DataFrame with same index as introns_df, columns are samples,
        values are cluster totals for each intron's cluster
    """
    logger.info(f"Computing cluster-total offsets using {cluster_col}...")
    
    # Calculate cluster totals per sample
    cluster_totals = introns_df.groupby(cluster_col)[sample_cols].sum()
    
    logger.info(f"Calculated totals for {len(cluster_totals)} clusters")
    
    # Log offset distribution statistics
    all_offsets = cluster_totals.values.flatten()
    all_offsets = all_offsets[all_offsets > 0]  # exclude zeros
    
    if len(all_offsets) > 0:
        logger.info(f"Offset distribution: min={all_offsets.min():.0f}, "
                    f"median={np.median(all_offsets):.0f}, "
                    f"mean={all_offsets.mean():.0f}, "
                    f"max={all_offsets.max():.0f}")
        
        # Check for very small offsets (potential issues)
        small_offsets = all_offsets[all_offsets < 10]
        if len(small_offsets) > 0:
            pct_small = 100.0 * len(small_offsets) / len(all_offsets)
            logger.warning(f"Warning: {len(small_offsets)} cluster-sample combinations "
                          f"({pct_small:.1f}%) have offsets < 10. Consider stricter filtering.")
    
    # Map cluster totals back to each intron
    offsets_df = pd.DataFrame(index=introns_df.index, columns=sample_cols, dtype=float)
    
    for cluster_id in cluster_totals.index:
        cluster_mask = introns_df[cluster_col] == cluster_id
        cluster_introns = introns_df.index[cluster_mask]
        
        for sample in sample_cols:
            offsets_df.loc[cluster_introns, sample] = cluster_totals.loc[cluster_id, sample]
    
    # Check for any missing offsets
    missing = offsets_df.isna().sum().sum()
    if missing > 0:
        logger.warning(f"Warning: {missing} missing offset values found")
    
    return offsets_df


def prepare_edgeR_input(introns_df, offsets_df, sample_cols, cluster_col, output_prefix):
    """
    Prepare input files for edgeR analysis.
    
    Creates:
    1. Count matrix (introns x samples)
    2. Offset matrix (introns x samples, log-transformed)
    3. Intron annotation file (intron metadata including cluster assignments)
    
    Args:
        introns_df: DataFrame with intron counts and annotations
        offsets_df: DataFrame with cluster-total offsets
        sample_cols: List of sample column names
        cluster_col: Name of cluster column
        output_prefix: Prefix for output files
    """
    logger.info("Preparing edgeR input files...")
    
    # 1. Count matrix
    count_matrix = introns_df[sample_cols]
    count_file = f"{output_prefix}.counts.tsv"
    count_matrix.to_csv(count_file, sep="\t")
    logger.info(f"Wrote count matrix: {count_file}")
    
    # 2. Offset matrix (log-transformed, handling zeros)
    # Add small pseudocount to avoid log(0)
    log_offsets = np.log(offsets_df + 0.5)
    offset_file = f"{output_prefix}.offsets.tsv"
    log_offsets.to_csv(offset_file, sep="\t")
    logger.info(f"Wrote offset matrix (log-transformed): {offset_file}")
    
    # 3. Intron annotations
    annotation_cols = [cluster_col]
    if "intron_info" in introns_df.columns:
        # Expand intron_info dict into separate columns
        info_df = pd.DataFrame(introns_df["intron_info"].tolist(), index=introns_df.index)
        annotation_df = pd.concat([info_df, introns_df[[cluster_col]]], axis=1)
    else:
        annotation_df = introns_df[[cluster_col]].copy()
        # Parse intron ID for basic info
        annotation_df["intron_id"] = annotation_df.index
    
    annotation_file = f"{output_prefix}.annotations.tsv"
    annotation_df.to_csv(annotation_file, sep="\t")
    logger.info(f"Wrote intron annotations: {annotation_file}")
    
    logger.info("edgeR input preparation complete!")
    
    return {
        "counts": count_file,
        "offsets": offset_file,
        "annotations": annotation_file,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compute cluster-total offsets for edgeR analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    parser.add_argument(
        "--matrix",
        type=str,
        required=True,
        help="Input filtered intron matrix (from filter_introns.py)",
    )
    
    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="Output prefix for edgeR input files (.counts.tsv, .offsets.tsv, .annotations.tsv)",
    )
    
    parser.add_argument(
        "--cluster_type",
        type=str,
        choices=["donor", "acceptor"],
        required=True,
        help="Type of clustering used",
    )
    
    args = parser.parse_args()
    
    # Load matrix
    logger.info(f"Loading filtered matrix from {args.matrix}")
    df = pd.read_csv(args.matrix, sep="\t", index_col=0)
    
    cluster_col = f"{args.cluster_type}_cluster"
    
    # Get sample columns
    exclude_cols = {"intron_info", "donor_cluster", "acceptor_cluster"}
    sample_cols = [col for col in df.columns if col not in exclude_cols]
    logger.info(f"Found {len(sample_cols)} samples")
    
    # Parse intron info if needed
    if "intron_info" not in df.columns:
        from cluster_introns import parse_intron_id
        logger.info("Parsing intron annotations...")
        df["intron_info"] = [parse_intron_id(idx) for idx in df.index]
    
    logger.info(f"Processing {len(df)} introns in {df[cluster_col].nunique()} clusters")
    
    # Compute offsets
    offsets_df = compute_cluster_offsets(df, sample_cols, cluster_col)
    
    # Prepare edgeR input files
    output_files = prepare_edgeR_input(df, offsets_df, sample_cols, cluster_col, args.output_prefix)
    
    logger.info("Offset computation complete!")
    logger.info(f"Output files:")
    for key, path in output_files.items():
        logger.info(f"  {key}: {path}")


if __name__ == "__main__":
    main()
