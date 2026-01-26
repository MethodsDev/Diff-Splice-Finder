#!/usr/bin/env python3

"""
Aggregate intron-level results to cluster-level splicing events.

Takes individual intron p-values and aggregates them per cluster using
combined p-value methods (e.g., Fisher's method, Cauchy combination).

This creates interpretable "splicing event" level results similar to
LeafCutter's cluster analysis, where each cluster represents a local
alternative splicing decision.
"""

import sys
import os
import argparse
import logging
import pandas as pd
import numpy as np
from scipy import stats

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def fisher_combine_pvalues(pvalues):
    """
    Combine p-values using Fisher's method.
    
    Args:
        pvalues: Array of p-values
        
    Returns:
        Combined p-value
    """
    # Remove any NaN or invalid p-values
    pvalues = pvalues[~np.isnan(pvalues)]
    pvalues = pvalues[(pvalues > 0) & (pvalues <= 1)]
    
    if len(pvalues) == 0:
        return np.nan
    
    # Fisher's method: -2 * sum(log(p)) ~ chi-squared with 2k df
    chi_stat = -2 * np.sum(np.log(pvalues))
    df = 2 * len(pvalues)
    combined_p = stats.chi2.sf(chi_stat, df)
    
    return combined_p


def cauchy_combine_pvalues(pvalues, weights=None):
    """
    Combine p-values using Cauchy combination test (ACAT).
    
    More robust than Fisher's method for very small p-values.
    Reference: Liu & Xie (2020) JASA
    
    Args:
        pvalues: Array of p-values
        weights: Optional weights for each p-value
        
    Returns:
        Combined p-value
    """
    # Remove any NaN or invalid p-values
    valid = ~np.isnan(pvalues) & (pvalues > 0) & (pvalues <= 1)
    pvalues = pvalues[valid]
    
    if weights is not None:
        weights = weights[valid]
    
    if len(pvalues) == 0:
        return np.nan
    
    if weights is None:
        weights = np.ones(len(pvalues)) / len(pvalues)
    else:
        weights = weights / np.sum(weights)
    
    # Transform p-values to Cauchy statistics
    # tan((0.5 - p) * pi) is the Cauchy transform
    cauchy_stats = np.tan((0.5 - pvalues) * np.pi)
    
    # Weighted sum
    test_stat = np.sum(weights * cauchy_stats)
    
    # Combined p-value
    combined_p = 0.5 - np.arctan(test_stat) / np.pi
    
    return max(0, min(1, combined_p))  # ensure valid p-value


def aggregate_cluster_results(intron_results, method="cauchy"):
    """
    Aggregate intron-level results to cluster-level.
    
    Args:
        intron_results: DataFrame with columns [intron_id, cluster, PValue, logFC, FDR, significant]
        method: P-value combination method ('fisher' or 'cauchy')
        
    Returns:
        DataFrame with cluster-level results
    """
    logger.info(f"Aggregating intron results to clusters using {method} method...")
    
    if "cluster" not in intron_results.columns:
        logger.error("Input results missing 'cluster' column")
        raise ValueError("Results must contain cluster assignments")
    
    # Group by cluster
    cluster_groups = intron_results.groupby("cluster")
    
    cluster_results = []
    
    for cluster_id, cluster_data in cluster_groups:
        n_introns = len(cluster_data)
        
        # Combine p-values
        pvalues = cluster_data["PValue"].values
        
        if method == "fisher":
            cluster_pval = fisher_combine_pvalues(pvalues)
        elif method == "cauchy":
            # Can weight by absolute logFC if desired
            cluster_pval = cauchy_combine_pvalues(pvalues)
        else:
            raise ValueError(f"Unknown method: {method}")
        
        # Summary statistics
        sig_introns = cluster_data["significant"].sum() if "significant" in cluster_data.columns else 0
        
        # Direction: dominant direction of significant introns
        sig_data = cluster_data[cluster_data["significant"]] if "significant" in cluster_data.columns else cluster_data
        
        if len(sig_data) > 0:
            n_up = (sig_data["logFC"] > 0).sum()
            n_down = (sig_data["logFC"] < 0).sum()
            
            if n_up > n_down:
                dominant_direction = "increased"
            elif n_down > n_up:
                dominant_direction = "decreased"
            else:
                dominant_direction = "mixed"
            
            max_abs_logfc = sig_data["logFC"].abs().max()
            mean_logfc = sig_data["logFC"].mean()
        else:
            dominant_direction = "none"
            max_abs_logfc = 0
            mean_logfc = 0
            n_up = 0
            n_down = 0
        
        cluster_results.append({
            "cluster": cluster_id,
            "n_introns": n_introns,
            "n_significant_introns": sig_introns,
            "cluster_pvalue": cluster_pval,
            "dominant_direction": dominant_direction,
            "n_increased": n_up,
            "n_decreased": n_down,
            "max_abs_logFC": max_abs_logfc,
            "mean_logFC": mean_logfc,
        })
    
    cluster_df = pd.DataFrame(cluster_results)
    
    # Calculate cluster-level FDR
    from statsmodels.stats.multitest import multipletests
    
    valid_pvals = ~cluster_df["cluster_pvalue"].isna()
    cluster_df["cluster_FDR"] = np.nan
    
    if valid_pvals.sum() > 0:
        _, cluster_fdr, _, _ = multipletests(
            cluster_df.loc[valid_pvals, "cluster_pvalue"],
            method="fdr_bh"
        )
        cluster_df.loc[valid_pvals, "cluster_FDR"] = cluster_fdr
    
    # Sort by p-value
    cluster_df = cluster_df.sort_values("cluster_pvalue")
    
    logger.info(f"Aggregated {len(cluster_df)} clusters from {len(intron_results)} introns")
    
    return cluster_df


def annotate_cluster_introns(intron_results, cluster_results):
    """
    Create detailed annotation of introns within significant clusters.
    
    Args:
        intron_results: DataFrame with intron-level results
        cluster_results: DataFrame with cluster-level results
        
    Returns:
        DataFrame with introns from significant clusters
    """
    # Get significant clusters
    sig_clusters = cluster_results[cluster_results["cluster_FDR"] < 0.05]["cluster"]
    
    if len(sig_clusters) == 0:
        logger.warning("No significant clusters found")
        return pd.DataFrame()
    
    # Filter introns to significant clusters
    sig_introns = intron_results[intron_results["cluster"].isin(sig_clusters)].copy()
    
    # Merge with cluster-level info
    sig_introns = sig_introns.merge(
        cluster_results[["cluster", "cluster_pvalue", "cluster_FDR", "dominant_direction"]],
        on="cluster",
        how="left"
    )
    
    # Sort by cluster FDR, then by intron significance
    sig_introns = sig_introns.sort_values(["cluster_FDR", "PValue"])
    
    logger.info(f"Found {len(sig_introns)} introns in {len(sig_clusters)} significant clusters")
    
    return sig_introns


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate intron-level results to cluster-level splicing events",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    parser.add_argument(
        "--intron_results",
        type=str,
        required=True,
        help="Input intron-level results from edgeR (.intron_results.tsv)",
    )
    
    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="Output prefix for cluster-level results",
    )
    
    parser.add_argument(
        "--method",
        type=str,
        choices=["fisher", "cauchy"],
        default="cauchy",
        help="P-value combination method",
    )
    
    parser.add_argument(
        "--cluster_fdr",
        type=float,
        default=0.05,
        help="FDR threshold for cluster-level significance",
    )
    
    args = parser.parse_args()
    
    # Load intron results
    logger.info(f"Loading intron results from {args.intron_results}")
    intron_results = pd.read_csv(args.intron_results, sep="\t")
    
    logger.info(f"Loaded {len(intron_results)} intron results")
    
    # Check for required columns
    required_cols = ["PValue", "logFC"]
    missing = [col for col in required_cols if col not in intron_results.columns]
    if missing:
        logger.error(f"Missing required columns: {missing}")
        sys.exit(1)
    
    # Aggregate to clusters
    cluster_results = aggregate_cluster_results(intron_results, method=args.method)
    
    # Write cluster-level results
    cluster_file = f"{args.output_prefix}.cluster_results.tsv"
    logger.info(f"Writing cluster results to {cluster_file}")
    cluster_results.to_csv(cluster_file, sep="\t", index=False)
    
    # Summary
    n_sig_clusters = (cluster_results["cluster_FDR"] < args.cluster_fdr).sum()
    logger.info(f"\n=== Cluster-Level Summary ===")
    logger.info(f"Total clusters: {len(cluster_results)}")
    logger.info(f"Significant clusters (FDR < {args.cluster_fdr}): {n_sig_clusters}")
    
    if n_sig_clusters > 0:
        sig_summary = cluster_results[cluster_results["cluster_FDR"] < args.cluster_fdr]
        logger.info(f"  Increased usage: {(sig_summary['dominant_direction'] == 'increased').sum()}")
        logger.info(f"  Decreased usage: {(sig_summary['dominant_direction'] == 'decreased').sum()}")
        logger.info(f"  Mixed: {(sig_summary['dominant_direction'] == 'mixed').sum()}")
        
        # Annotate introns in significant clusters
        sig_introns = annotate_cluster_introns(intron_results, cluster_results)
        
        if len(sig_introns) > 0:
            sig_introns_file = f"{args.output_prefix}.significant_cluster_introns.tsv"
            logger.info(f"Writing introns from significant clusters to {sig_introns_file}")
            sig_introns.to_csv(sig_introns_file, sep="\t", index=False)
    
    logger.info("\nAggregation complete!")


if __name__ == "__main__":
    main()
