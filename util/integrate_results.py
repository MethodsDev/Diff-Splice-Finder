#!/usr/bin/env python3

"""
Integrate results from donor and acceptor clustering analyses.

Creates unified summary files combining both analyses, allowing users to:
- See all tested introns across both clustering strategies
- Identify introns significant in donor, acceptor, or both
- Compare effect sizes between clustering methods
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


def load_results(donor_file, acceptor_file):
    """
    Load donor and acceptor results files.
    
    Args:
        donor_file: Path to donor intron results
        acceptor_file: Path to acceptor intron results
        
    Returns:
        Tuple of (donor_df, acceptor_df)
    """
    logger.info(f"Loading donor results from {donor_file}")
    donor_df = pd.read_csv(donor_file, sep="\t")
    
    logger.info(f"Loading acceptor results from {acceptor_file}")
    acceptor_df = pd.read_csv(acceptor_file, sep="\t")
    
    logger.info(f"Loaded {len(donor_df)} donor introns, {len(acceptor_df)} acceptor introns")
    
    return donor_df, acceptor_df


def integrate_intron_results(donor_df, acceptor_df):
    """
    Integrate donor and acceptor intron-level results.
    
    Args:
        donor_df: DataFrame with donor results
        acceptor_df: DataFrame with acceptor results
        
    Returns:
        Integrated DataFrame with columns indicating significance in each analysis
    """
    logger.info("Integrating intron-level results...")
    
    # Add clustering type indicator
    donor_df = donor_df.copy()
    acceptor_df = acceptor_df.copy()
    
    # Rename cluster-specific columns
    donor_df = donor_df.rename(columns={
        'logFC': 'donor_logFC',
        'logCPM': 'donor_logCPM',
        'F': 'donor_F',
        'PValue': 'donor_PValue',
        'FDR': 'donor_FDR',
        'cluster': 'donor_cluster',
        'significant': 'donor_significant'
    })
    
    acceptor_df = acceptor_df.rename(columns={
        'logFC': 'acceptor_logFC',
        'logCPM': 'acceptor_logCPM',
        'F': 'acceptor_F',
        'PValue': 'acceptor_PValue',
        'FDR': 'acceptor_FDR',
        'cluster': 'acceptor_cluster',
        'significant': 'acceptor_significant'
    })
    
    # Merge on intron_id (outer join to get all introns)
    integrated = pd.merge(
        donor_df,
        acceptor_df,
        on='intron_id',
        how='outer',
        suffixes=('', '_dup')
    )
    
    # Create summary columns
    integrated['tested_in'] = np.where(
        integrated['donor_PValue'].notna() & integrated['acceptor_PValue'].notna(),
        'both',
        np.where(integrated['donor_PValue'].notna(), 'donor_only', 'acceptor_only')
    )
    
    # Significance status
    donor_sig = integrated['donor_significant'].fillna(False)
    acceptor_sig = integrated['acceptor_significant'].fillna(False)
    
    integrated['significant_in'] = np.where(
        donor_sig & acceptor_sig,
        'both',
        np.where(
            donor_sig,
            'donor_only',
            np.where(acceptor_sig, 'acceptor_only', 'neither')
        )
    )
    
    # For introns significant in both, check direction consistency
    both_sig = (donor_sig & acceptor_sig)
    if both_sig.any():
        same_direction = np.sign(integrated.loc[both_sig, 'donor_logFC']) == \
                        np.sign(integrated.loc[both_sig, 'acceptor_logFC'])
        integrated.loc[both_sig, 'direction_consistent'] = same_direction
    else:
        integrated['direction_consistent'] = np.nan
    
    # Best (most significant) analysis
    integrated['best_FDR'] = integrated[['donor_FDR', 'acceptor_FDR']].min(axis=1)
    integrated['best_analysis'] = np.where(
        integrated['donor_FDR'] <= integrated['acceptor_FDR'],
        'donor',
        'acceptor'
    )
    
    # For best analysis, get the corresponding logFC
    integrated['best_logFC'] = np.where(
        integrated['best_analysis'] == 'donor',
        integrated['donor_logFC'],
        integrated['acceptor_logFC']
    )
    
    # Reorder columns for readability
    priority_cols = [
        'intron_id',
        'tested_in',
        'significant_in',
        'best_analysis',
        'best_FDR',
        'best_logFC',
        'direction_consistent',
        'donor_cluster',
        'donor_logFC',
        'donor_PValue',
        'donor_FDR',
        'donor_significant',
        'acceptor_cluster',
        'acceptor_logFC',
        'acceptor_PValue',
        'acceptor_FDR',
        'acceptor_significant',
    ]
    
    # Keep priority columns first, then any remaining columns
    other_cols = [col for col in integrated.columns if col not in priority_cols]
    integrated = integrated[priority_cols + other_cols]
    
    # Sort by best FDR
    integrated = integrated.sort_values('best_FDR')
    
    logger.info(f"Integrated {len(integrated)} unique introns")
    
    # Summary statistics
    logger.info(f"Tested in both analyses: {(integrated['tested_in'] == 'both').sum()}")
    logger.info(f"Tested in donor only: {(integrated['tested_in'] == 'donor_only').sum()}")
    logger.info(f"Tested in acceptor only: {(integrated['tested_in'] == 'acceptor_only').sum()}")
    
    sig_counts = integrated['significant_in'].value_counts()
    logger.info(f"\nSignificance summary:")
    for status, count in sig_counts.items():
        logger.info(f"  {status}: {count}")
    
    # Direction consistency for those significant in both
    both_sig_count = (integrated['significant_in'] == 'both').sum()
    if both_sig_count > 0:
        both_sig_mask = integrated['significant_in'] == 'both'
        consistent_count = (integrated.loc[both_sig_mask, 'direction_consistent'] == True).sum()
        opposite_count = (integrated.loc[both_sig_mask, 'direction_consistent'] == False).sum()
        logger.info(f"\nOf {both_sig_count} introns significant in both:")
        logger.info(f"  Consistent direction: {consistent_count}")
        logger.info(f"  Opposite direction: {opposite_count}")
    
    return integrated


def create_summary_table(integrated_df):
    """
    Create a high-level summary table of integration results.
    
    Args:
        integrated_df: Integrated results DataFrame
        
    Returns:
        Summary DataFrame
    """
    summary = {
        'total_introns': len(integrated_df),
        'tested_both': (integrated_df['tested_in'] == 'both').sum(),
        'tested_donor_only': (integrated_df['tested_in'] == 'donor_only').sum(),
        'tested_acceptor_only': (integrated_df['tested_in'] == 'acceptor_only').sum(),
        'significant_both': (integrated_df['significant_in'] == 'both').sum(),
        'significant_donor_only': (integrated_df['significant_in'] == 'donor_only').sum(),
        'significant_acceptor_only': (integrated_df['significant_in'] == 'acceptor_only').sum(),
        'significant_neither': (integrated_df['significant_in'] == 'neither').sum(),
    }
    
    # Direction consistency for those significant in both
    both_sig = integrated_df[integrated_df['significant_in'] == 'both']
    if len(both_sig) > 0:
        # Count True/False explicitly to avoid issues with NaN and boolean negation
        consistent = (both_sig['direction_consistent'] == True).sum()
        opposite = (both_sig['direction_consistent'] == False).sum()
        summary['both_sig_consistent_direction'] = int(consistent)
        summary['both_sig_opposite_direction'] = int(opposite)
    else:
        summary['both_sig_consistent_direction'] = 0
        summary['both_sig_opposite_direction'] = 0
    
    return pd.DataFrame([summary])


def main():
    parser = argparse.ArgumentParser(
        description="Integrate donor and acceptor differential splicing results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    parser.add_argument(
        "--donor_results",
        type=str,
        required=True,
        help="Donor intron results file (.intron_results.tsv)",
    )
    
    parser.add_argument(
        "--acceptor_results",
        type=str,
        required=True,
        help="Acceptor intron results file (.intron_results.tsv)",
    )
    
    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="Output prefix for integrated results",
    )
    
    parser.add_argument(
        "--fdr_threshold",
        type=float,
        default=0.05,
        help="FDR threshold for filtering significant results",
    )
    
    args = parser.parse_args()
    
    # Load results
    donor_df, acceptor_df = load_results(args.donor_results, args.acceptor_results)
    
    # Integrate
    integrated = integrate_intron_results(donor_df, acceptor_df)
    
    # Write all integrated results
    all_results_file = f"{args.output_prefix}.integrated_results.tsv"
    logger.info(f"\nWriting all integrated results to {all_results_file}")
    integrated.to_csv(all_results_file, sep="\t", index=False)
    
    # Filter to significant introns (in at least one analysis)
    significant = integrated[integrated['significant_in'] != 'neither']
    
    if len(significant) > 0:
        sig_results_file = f"{args.output_prefix}.significant_integrated.tsv"
        logger.info(f"Writing {len(significant)} significant introns to {sig_results_file}")
        significant.to_csv(sig_results_file, sep="\t", index=False)
        
        # Create summary table
        summary = create_summary_table(integrated)
        summary_file = f"{args.output_prefix}.integration_summary.tsv"
        logger.info(f"Writing integration summary to {summary_file}")
        summary.to_csv(summary_file, sep="\t", index=False)
    else:
        logger.warning("No significant introns found in either analysis")
    
    logger.info("\nIntegration complete!")


if __name__ == "__main__":
    main()
