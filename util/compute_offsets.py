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
    
    IMPORTANT: Always uses max(donor_cluster_total, acceptor_cluster_total) for
    ALL analyses. This ensures:
    1. The same intron has identical offsets in donor and acceptor analyses
    2. Novel splice sites paired with common sites get appropriate offsets
    3. The offset reflects the expression potential at BOTH ends of the intron
    
    Args:
        introns_df: DataFrame with intron counts and cluster assignments
        sample_cols: List of sample column names
        cluster_col: Name of cluster column ('donor_cluster' or 'acceptor_cluster')
                     or None to compute from both donor and acceptor clusters
        
    Returns:
        DataFrame with same index as introns_df, columns are samples,
        values are max(donor_total, acceptor_total) for each intron
    """
    logger.info(f"Computing cluster-total offsets (will use max of donor and acceptor)...")
    
    # Check if we have both cluster types
    has_donor = 'donor_cluster' in introns_df.columns
    has_acceptor = 'acceptor_cluster' in introns_df.columns
    
    # If cluster_col is None, we want to compute from both donor and acceptor
    if cluster_col is None:
        if not (has_donor and has_acceptor):
            raise ValueError("cluster_col=None requires both donor_cluster and acceptor_cluster columns in dataframe")
        # Force use of both clusters
        use_both = True
    else:
        # Specific cluster type requested, but still use both if available
        use_both = has_donor and has_acceptor
    
    if not use_both:
        # Fall back to single cluster type
        if cluster_col is None:
            raise ValueError("No cluster columns found in dataframe")
        logger.warning(f"Missing donor_cluster or acceptor_cluster column - using only {cluster_col}")
        cluster_totals = introns_df.groupby(cluster_col)[sample_cols].sum()
        logger.info(f"Calculated totals for {len(cluster_totals)} {cluster_col}s")
        
        # Map cluster totals back to each intron
        offsets_df = pd.DataFrame(index=introns_df.index, columns=sample_cols, dtype=float)
        for cluster_id in cluster_totals.index:
            cluster_mask = introns_df[cluster_col] == cluster_id
            cluster_introns = introns_df.index[cluster_mask]
            for sample in sample_cols:
                offsets_df.loc[cluster_introns, sample] = cluster_totals.loc[cluster_id, sample]
        return offsets_df
    
    # Calculate both donor and acceptor cluster totals
    logger.info(f"Computing donor_cluster totals...")
    donor_cluster_totals = introns_df.groupby('donor_cluster')[sample_cols].sum()
    logger.info(f"  Found {len(donor_cluster_totals)} donor clusters")
    
    logger.info(f"Computing acceptor_cluster totals...")
    acceptor_cluster_totals = introns_df.groupby('acceptor_cluster')[sample_cols].sum()
    logger.info(f"  Found {len(acceptor_cluster_totals)} acceptor clusters")
    
    # Map donor cluster totals to each intron (vectorized - much faster!)
    logger.info(f"Mapping donor cluster totals to introns...")
    donor_offsets_df = introns_df[['donor_cluster']].join(donor_cluster_totals, on='donor_cluster')[sample_cols]
    
    # Map acceptor cluster totals to each intron (vectorized - much faster!)
    logger.info(f"Mapping acceptor cluster totals to introns...")
    acceptor_offsets_df = introns_df[['acceptor_cluster']].join(acceptor_cluster_totals, on='acceptor_cluster')[sample_cols]
    
    # Log individual offset distributions
    donor_vals = donor_offsets_df.values.flatten()
    donor_vals = donor_vals[donor_vals > 0]
    if len(donor_vals) > 0:
        logger.info(f"Donor cluster offsets: min={donor_vals.min():.0f}, "
                    f"median={np.median(donor_vals):.0f}, "
                    f"mean={donor_vals.mean():.0f}, "
                    f"max={donor_vals.max():.0f}")
    
    acceptor_vals = acceptor_offsets_df.values.flatten()
    acceptor_vals = acceptor_vals[acceptor_vals > 0]
    if len(acceptor_vals) > 0:
        logger.info(f"Acceptor cluster offsets: min={acceptor_vals.min():.0f}, "
                    f"median={np.median(acceptor_vals):.0f}, "
                    f"mean={acceptor_vals.mean():.0f}, "
                    f"max={acceptor_vals.max():.0f}")
    
    # Use max(donor, acceptor) for all introns
    # This ensures the same intron has identical offsets in both analyses
    logger.info(f"Using max(donor_total, acceptor_total) for all offsets...")
    offsets_df = donor_offsets_df.combine(acceptor_offsets_df, np.maximum, fill_value=0)
    
    # Log final offset distribution
    final_vals = offsets_df.values.flatten()
    final_vals = final_vals[final_vals > 0]
    if len(final_vals) > 0:
        logger.info(f"Final offsets (max): min={final_vals.min():.0f}, "
                    f"median={np.median(final_vals):.0f}, "
                    f"mean={final_vals.mean():.0f}, "
                    f"max={final_vals.max():.0f}")
        
        # Count how many times each source was used
        used_donor = (donor_offsets_df >= acceptor_offsets_df).sum().sum()
        used_acceptor = (acceptor_offsets_df > donor_offsets_df).sum().sum()
        total = used_donor + used_acceptor
        logger.info(f"Offset source: donor={used_donor}/{total} ({100*used_donor/total:.1f}%), "
                   f"acceptor={used_acceptor}/{total} ({100*used_acceptor/total:.1f}%)")
        
        # Check for very small offsets
        small_offsets = final_vals[final_vals < 10]
        if len(small_offsets) > 0:
            pct_small = 100.0 * len(small_offsets) / len(final_vals)
            logger.warning(f"Warning: {len(small_offsets)} intron-sample combinations "
                          f"({pct_small:.1f}%) have offsets < 10. Consider stricter filtering.")
    
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
    count_matrix.to_csv(count_file, sep="\t", na_rep='NA')
    logger.info(f"Wrote count matrix: {count_file}")
    
    # 2. Offset matrix (log-transformed, handling zeros)
    # Filter offsets to match sample_cols (in case samples were filtered)
    offsets_filtered = offsets_df[sample_cols]
    # Add small pseudocount to avoid log(0)
    log_offsets = np.log(offsets_filtered + 0.5)
    offset_file = f"{output_prefix}.offsets.tsv"
    log_offsets.to_csv(offset_file, sep="\t", na_rep='NA')
    logger.info(f"Wrote offset matrix (log-transformed): {offset_file}")
    
    # 3. Intron annotations
    annotation_cols = [cluster_col]
    
    # Add both cluster columns if available (for intron-level analysis)
    if 'donor_cluster' in introns_df.columns and cluster_col != 'donor_cluster':
        annotation_cols.append('donor_cluster')
    if 'acceptor_cluster' in introns_df.columns and cluster_col != 'acceptor_cluster':
        annotation_cols.append('acceptor_cluster')
    
    # Add gene_name and intron_status if available
    metadata_cols = ['gene_name', 'intron_status', 'overlapping_genes']
    for col in metadata_cols:
        if col in introns_df.columns:
            annotation_cols.append(col)
    
    if "intron_info" in introns_df.columns:
        # Expand intron_info dict into separate columns
        info_df = pd.DataFrame(introns_df["intron_info"].tolist(), index=introns_df.index)
        annotation_df = pd.concat([info_df, introns_df[annotation_cols]], axis=1)
    else:
        annotation_df = introns_df[annotation_cols].copy()
        # Parse intron ID for basic info
        annotation_df["intron_id"] = annotation_df.index
    
    annotation_file = f"{output_prefix}.annotations.tsv"
    annotation_df.to_csv(annotation_file, sep="\t", na_rep='NA')
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
        help="Input intron matrix (filtered or full clustered)",
    )
    
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output file for offsets only (used with --compute_offsets_only)",
    )
    
    parser.add_argument(
        "--output_prefix",
        type=str,
        default=None,
        help="Output prefix for edgeR input files (.counts.tsv, .offsets.tsv, .annotations.tsv)",
    )
    
    parser.add_argument(
        "--compute_offsets_only",
        action="store_true",
        help="Only compute and save offsets (no counts or annotations). Requires --output.",
    )
    
    parser.add_argument(
        "--shared_offsets",
        type=str,
        default=None,
        help="Use pre-computed shared offsets file instead of computing new offsets",
    )
    
    parser.add_argument(
        "--samples",
        type=str,
        default=None,
        help="Sample metadata file - if provided, filter matrices to only include these samples",
    )
    
    parser.add_argument(
        "--cluster_type",
        type=str,
        choices=["donor", "acceptor"],
        default=None,
        help="Type of clustering used (required unless --compute_offsets_only)",
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.compute_offsets_only:
        if not args.output:
            parser.error("--compute_offsets_only requires --output")
    else:
        if not args.output_prefix:
            parser.error("--output_prefix required unless using --compute_offsets_only")
        # cluster_type only needed if not using shared offsets
        if not args.shared_offsets and not args.cluster_type:
            parser.error("--cluster_type required unless using --shared_offsets")
    
    # Load matrix
    logger.info(f"Loading matrix from {args.matrix}")
    df = pd.read_csv(args.matrix, sep="\t", index_col=0)
    logger.info(f"Loaded {len(df)} introns")
    
    # Get sample columns
    exclude_cols = {"intron_info", "donor_cluster", "acceptor_cluster", "gene_name", "intron_status", "overlapping_genes"}
    sample_cols = [col for col in df.columns if col not in exclude_cols]
    logger.info(f"Found {len(sample_cols)} samples")
    
    # Filter samples if metadata file provided
    if args.samples:
        logger.info(f"Loading sample metadata from {args.samples}")
        samples_df = pd.read_csv(args.samples, sep="\t", comment="#")
        
        if 'sample_id' not in samples_df.columns:
            raise ValueError("Sample metadata file must have a 'sample_id' column")
        
        # Convert sample_id to string and filter out NaN/empty values
        metadata_samples = set(
            str(x) for x in samples_df['sample_id'].tolist() 
            if pd.notna(x) and str(x).strip()
        )
        available_samples = set(sample_cols)
        
        # Find samples in metadata that exist in the matrix
        samples_to_keep = metadata_samples & available_samples
        samples_missing = metadata_samples - available_samples
        samples_excluded = available_samples - metadata_samples
        
        if not samples_to_keep:
            raise ValueError("No samples from metadata file found in count matrix")
        
        if samples_missing:
            missing_list = sorted(samples_missing)[:5]
            logger.warning(f"Warning: {len(samples_missing)} samples in metadata not found in matrix: {', '.join(missing_list)}{'...' if len(samples_missing) > 5 else ''}")
        
        if samples_excluded:
            excluded_list = sorted(samples_excluded)[:5]
            logger.info(f"Filtering out {len(samples_excluded)} samples not in metadata: {', '.join(excluded_list)}{'...' if len(samples_excluded) > 5 else ''}")
        
        # Filter to only the samples to keep
        sample_cols = [s for s in sample_cols if s in samples_to_keep]
        logger.info(f"Keeping {len(sample_cols)} samples that are in both metadata and matrix")

    
    # MODE 1: Compute offsets only (for shared offsets file)
    if args.compute_offsets_only:
        logger.info("Computing shared offsets from full clustered matrix...")
        
        # Parse intron info if needed
        if "intron_info" not in df.columns:
            from cluster_introns import parse_intron_id
            logger.info("Parsing intron annotations...")
            df["intron_info"] = [parse_intron_id(idx) for idx in df.index]
        
        # Compute offsets (uses both donor and acceptor clusters)
        offsets_df = compute_cluster_offsets(df, sample_cols, None)  # None means use both
        
        # Save offsets
        logger.info(f"Writing shared offsets to {args.output}")
        offsets_df.to_csv(args.output, sep="\t", na_rep='NA')
        logger.info("Shared offset computation complete!")
        return
    
    # MODE 2: Prepare edgeR inputs using shared offsets
    if args.shared_offsets:
        logger.info(f"Using pre-computed shared offsets from {args.shared_offsets}")
        offsets_full = pd.read_csv(args.shared_offsets, sep="\t", index_col=0)
        
        # Subset offsets to match filtered introns
        offsets_df = offsets_full.loc[df.index]
        logger.info(f"Subset offsets to {len(offsets_df)} introns")
        
        # Parse intron info if needed (for prepare_edgeR_input)
        if "intron_info" not in df.columns:
            from cluster_introns import parse_intron_id
            logger.info("Parsing intron annotations...")
            df["intron_info"] = [parse_intron_id(idx) for idx in df.index]
    else:
        # MODE 3: Legacy mode - compute offsets from filtered matrix (deprecated)
        logger.warning("Computing offsets from filtered matrix - consider using shared offsets for consistency")
        
        cluster_col = f"{args.cluster_type}_cluster"
        
        # Parse intron info if needed
        if "intron_info" not in df.columns:
            from cluster_introns import parse_intron_id
            logger.info("Parsing intron annotations...")
            df["intron_info"] = [parse_intron_id(idx) for idx in df.index]
        
        logger.info(f"Processing {len(df)} introns in {df[cluster_col].nunique()} clusters")
        offsets_df = compute_cluster_offsets(df, sample_cols, cluster_col)
    
    # For edgeR input preparation, determine cluster_col to use for annotations
    # Use whichever cluster columns are available
    if args.cluster_type:
        cluster_col = f"{args.cluster_type}_cluster"
    elif 'donor_cluster' in df.columns and 'acceptor_cluster' in df.columns:
        logger.info("Both donor_cluster and acceptor_cluster columns found - using donor_cluster for annotations")
        cluster_col = 'donor_cluster'
    elif 'donor_cluster' in df.columns:
        cluster_col = 'donor_cluster'
    elif 'acceptor_cluster' in df.columns:
        cluster_col = 'acceptor_cluster'
    else:
        raise ValueError("No cluster columns found in dataframe")
    
    # Parse intron info if needed (for prepare_edgeR_input)
    if "intron_info" not in df.columns:
        from cluster_introns import parse_intron_id
        logger.info("Parsing intron annotations...")
        df["intron_info"] = [parse_intron_id(idx) for idx in df.index]
    
    # Prepare edgeR input files
    output_files = prepare_edgeR_input(df, offsets_df, sample_cols, cluster_col, args.output_prefix)
    
    logger.info("Offset computation complete!")
    logger.info(f"Output files:")
    for key, path in output_files.items():
        logger.info(f"  {key}: {path}")


if __name__ == "__main__":
    main()
