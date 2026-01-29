#!/usr/bin/env python3

"""
Compute PSI (Percent Spliced In) values for introns.

PSI represents the proportion of reads supporting a particular intron
relative to all introns in its cluster (shared donor or acceptor site).

PSI_i,s = count_i,s / sum(count_j,s for all j in cluster_i)

where:
- PSI_i,s is the PSI for intron i in sample s
- count_i,s is the read count for intron i in sample s
- cluster_i is the set of all introns sharing a splice site with intron i
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


def compute_psi_values(counts_df, annotations_df, sample_metadata, cluster_col='donor_cluster', group_col='group', contrast=None, shared_cluster_totals=None):
    """
    Compute PSI values for each intron across samples.
    
    Args:
        counts_df: DataFrame with intron counts (introns x samples)
        annotations_df: DataFrame with cluster assignments
        sample_metadata: DataFrame with sample information including group assignments
        cluster_col: Column name for cluster assignment
        group_col: Column name for sample groups in metadata
        contrast: Optional contrast string (e.g., 'GroupA-GroupB') to compute delta PSI
        shared_cluster_totals: Optional DataFrame with shared offsets (max of donor/acceptor cluster totals)
                               to use as consistent denominators for PSI calculation
        
    Returns:
        DataFrame with PSI values and summary statistics
    """
    # Get sample columns (exclude metadata columns)
    metadata_cols = {'gene_name', 'intron_status', 'overlapping_genes', cluster_col}
    sample_cols = [col for col in counts_df.columns if col not in metadata_cols]
    
    # Match samples to groups
    sample_to_group = dict(zip(sample_metadata['sample_id'], sample_metadata[group_col]))
    
    # Determine denominators for PSI calculation
    if shared_cluster_totals is not None:
        logger.info(f"Computing PSI values using shared offsets (max of donor and acceptor cluster totals)...")
        # Use shared cluster totals (already max of donor/acceptor) for all introns
        # Subset to match introns in counts_df
        cluster_totals_df = shared_cluster_totals.loc[counts_df.index, sample_cols]
    else:
        logger.info(f"Computing PSI values using {cluster_col}...")
        # Fall back to cluster-specific totals (legacy mode)
        counts_with_cluster = counts_df.copy()
        counts_with_cluster[cluster_col] = annotations_df.loc[counts_df.index, cluster_col]
        
        # Compute cluster totals for each sample
        cluster_totals_grouped = counts_with_cluster.groupby(cluster_col)[sample_cols].sum()
        
        # Map back to each intron
        cluster_totals_df = pd.DataFrame(index=counts_df.index, columns=sample_cols, dtype=float)
        for sample in sample_cols:
            intron_cluster_totals = counts_with_cluster[cluster_col].map(
                cluster_totals_grouped[sample].to_dict()
            )
            cluster_totals_df[sample] = intron_cluster_totals
    
    # Compute PSI for each intron in each sample
    psi_df = pd.DataFrame(index=counts_df.index)
    
    for sample in sample_cols:
        # PSI = intron_count / cluster_total
        # Only compute PSI where cluster_total > 0 to avoid division by zero
        # Set PSI to 0 where cluster has no reads
        psi_df[f'{sample}_PSI'] = np.where(
            cluster_totals_df[sample] > 0,
            counts_df[sample] / cluster_totals_df[sample],
            0.0
        )
    
    # Compute group-level statistics
    groups = sample_metadata[group_col].unique()
    
    for group in groups:
        group_samples = sample_metadata[sample_metadata[group_col] == group]['sample_id'].tolist()
        group_psi_cols = [f'{s}_PSI' for s in group_samples if f'{s}_PSI' in psi_df.columns]
        
        if group_psi_cols:
            psi_df[f'{group}_mean_PSI'] = psi_df[group_psi_cols].mean(axis=1)
            psi_df[f'{group}_median_PSI'] = psi_df[group_psi_cols].median(axis=1)
            psi_df[f'{group}_std_PSI'] = psi_df[group_psi_cols].std(axis=1)
    
    # Compute delta PSI based on contrast if provided
    if contrast:
        # Parse contrast (format: 'GroupA-GroupB')
        parts = contrast.split('-')
        if len(parts) == 2:
            group1, group2 = parts[0].strip(), parts[1].strip()
            if f'{group1}_mean_PSI' in psi_df.columns and f'{group2}_mean_PSI' in psi_df.columns:
                # delta_PSI = group1 - group2 (matching logFC direction)
                psi_df['delta_PSI'] = psi_df[f'{group1}_mean_PSI'] - psi_df[f'{group2}_mean_PSI']
                logger.info(f"Computed delta_PSI for contrast: {contrast} ({group1} - {group2})")
            else:
                logger.warning(f"Could not compute delta_PSI: groups {group1} or {group2} not found")
    elif len(groups) == 2:
        # If no contrast specified but exactly 2 groups, compute delta PSI
        group1, group2 = sorted(groups)
        psi_df['delta_PSI'] = psi_df[f'{group1}_mean_PSI'] - psi_df[f'{group2}_mean_PSI']
        logger.info(f"Computed delta_PSI: {group1} - {group2}")
    
    logger.info(f"Computed PSI for {len(psi_df)} introns across {len(sample_cols)} samples")
    
    # Summary statistics
    if 'delta_PSI' in psi_df.columns:
        high_delta = (psi_df['delta_PSI'].abs() > 0.1).sum()
        logger.info(f"Introns with |ΔPSI| > 0.1: {high_delta} ({100*high_delta/len(psi_df):.1f}%)")
    
    return psi_df


def add_psi_to_results(results_df, psi_df, keep_individual_samples=False):
    """
    Add PSI columns to edgeR results DataFrame.
    
    Args:
        results_df: DataFrame with edgeR results
        psi_df: DataFrame with PSI values
        keep_individual_samples: If True, keep per-sample PSI columns
        
    Returns:
        Enhanced DataFrame with PSI information
    """
    logger.info("Adding PSI values to results...")
    
    # Determine which PSI columns to keep
    if keep_individual_samples:
        psi_cols = list(psi_df.columns)
    else:
        # Keep only summary columns (mean, median, std, delta)
        psi_cols = [col for col in psi_df.columns 
                   if 'mean_PSI' in col or 'median_PSI' in col or 
                      'std_PSI' in col or 'delta_PSI' in col]
    
    # Merge PSI data with results
    enhanced_df = results_df.copy()
    
    for col in psi_cols:
        if col in psi_df.columns:
            enhanced_df[col] = psi_df.loc[results_df.index, col]
    
    logger.info(f"Added {len(psi_cols)} PSI columns to results")
    
    return enhanced_df


def main():
    parser = argparse.ArgumentParser(
        description="Compute PSI values for differential splicing results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    parser.add_argument(
        "--counts",
        type=str,
        required=True,
        help="Input count matrix (introns x samples)",
    )
    
    parser.add_argument(
        "--annotations",
        type=str,
        required=True,
        help="Intron annotations with cluster assignments",
    )
    
    parser.add_argument(
        "--samples",
        type=str,
        required=True,
        help="Sample metadata file",
    )
    
    parser.add_argument(
        "--results",
        type=str,
        required=True,
        help="edgeR results file to enhance",
    )
    
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output file with PSI-enhanced results",
    )
    
    parser.add_argument(
        "--cluster_col",
        type=str,
        default="donor_cluster",
        help="Column name for cluster assignments",
    )
    
    parser.add_argument(
        "--group_col",
        type=str,
        default="group",
        help="Column name for sample groups",
    )
    
    parser.add_argument(
        "--keep_individual_samples",
        action="store_true",
        help="Keep per-sample PSI values in output",
    )
    
    args = parser.parse_args()
    
    # Load data
    logger.info(f"Loading counts from {args.counts}")
    counts_df = pd.read_csv(args.counts, sep="\t", index_col=0)
    
    logger.info(f"Loading annotations from {args.annotations}")
    annotations_df = pd.read_csv(args.annotations, sep="\t", index_col=0)
    
    logger.info(f"Loading sample metadata from {args.samples}")
    samples_df = pd.read_csv(args.samples, sep="\t", comment='#')
    
    logger.info(f"Loading results from {args.results}")
    results_df = pd.read_csv(args.results, sep="\t")
    results_df = results_df.set_index('intron_id')
    
    # Compute PSI
    psi_df = compute_psi_values(
        counts_df, 
        annotations_df, 
        samples_df,
        cluster_col=args.cluster_col,
        group_col=args.group_col,
        contrast=None  # Will be extracted from results if needed
    )
    
    # Add PSI to results
    enhanced_results = add_psi_to_results(
        results_df, 
        psi_df,
        keep_individual_samples=args.keep_individual_samples
    )
    
    # Write output
    logger.info(f"Writing PSI-enhanced results to {args.output}")
    enhanced_results.reset_index().to_csv(args.output, sep="\t", index=False, na_rep='NA')
    
    logger.info("PSI computation complete!")


if __name__ == "__main__":
    main()
