#!/usr/bin/env python3

"""
Main pipeline orchestrator for differential splicing analysis.

Coordinates the full workflow:
1. Load intron count matrix
2. Cluster introns by donor/acceptor sites
3. Filter introns and clusters
4. Compute cluster-total offsets
5. Run edgeR analysis
6. Aggregate results to cluster level

Runs analysis for both donor and acceptor clustering by default.
"""

import sys
import os
import argparse
import logging
import subprocess
import pandas as pd
from multiprocessing import Pool
from itertools import combinations
from functools import partial

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def file_exists_and_valid(filepath):
    """
    Check if a file exists and is non-empty.
    
    Args:
        filepath: Path to check
        
    Returns:
        True if file exists and has size > 0
    """
    return os.path.exists(filepath) and os.path.getsize(filepath) > 0


def run_command(cmd, description, skip_if_exists=None):
    """
    Run a shell command and handle errors.
    
    Args:
        cmd: Command string or list
        description: Description for logging
        skip_if_exists: Optional path to check; if exists, skip command
        
    Returns:
        Result object or None if skipped
    """
    # Check if we should skip this step
    if skip_if_exists and file_exists_and_valid(skip_if_exists):
        logger.info(f"=== {description} ===")
        logger.info(f"SKIPPING - Output already exists: {skip_if_exists}")
        return None
    
    logger.info(f"=== {description} ===")
    logger.info(f"Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    
    # Run command with real-time output to terminal
    process = subprocess.Popen(
        cmd,
        shell=isinstance(cmd, str),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1
    )
    
    # Stream output in real-time
    output_lines = []
    for line in process.stdout:
        print(line, end='')  # Print to terminal in real-time
        output_lines.append(line)
    
    process.wait()
    
    if process.returncode != 0:
        logger.error(f"Command failed with return code {process.returncode}")
        raise RuntimeError(f"Failed: {description}")
    
    # Create a result object similar to subprocess.run
    class Result:
        def __init__(self, returncode, stdout):
            self.returncode = returncode
            self.stdout = stdout
            self.stderr = ""
    
    return Result(process.returncode, ''.join(output_lines))


def cluster_and_filter(matrix_file, cluster_type, output_dir, filter_params, force_rerun=False):
    """
    Cluster introns and apply filtering.
    
    Args:
        matrix_file: Input intron count matrix
        cluster_type: 'donor' or 'acceptor'
        output_dir: Output directory
        filter_params: Dict with filtering parameters
        force_rerun: If True, rerun even if outputs exist
        
    Returns:
        Path to filtered matrix file
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    
    # Step 1: Cluster introns
    clustered_file = os.path.join(output_dir, f"{cluster_type}_clustered.tsv")
    
    cmd = [
        "python3",
        os.path.join(util_dir, "cluster_introns.py"),
        "--matrix", matrix_file,
        f"--output_{cluster_type}", clustered_file,
        "--cluster_type", cluster_type,
    ]
    
    run_command(
        cmd, 
        f"Clustering introns by {cluster_type}",
        skip_if_exists=None if force_rerun else clustered_file
    )
    
    # Step 2: Filter
    filtered_file = os.path.join(output_dir, f"{cluster_type}_filtered.tsv")
    
    cmd = [
        "python3",
        os.path.join(util_dir, "filter_introns.py"),
        "--matrix", clustered_file,
        "--output", filtered_file,
        "--cluster_type", cluster_type,
        "--min_intron_count", str(filter_params["min_intron_count"]),
        "--min_intron_samples", str(filter_params["min_intron_samples"]),
        "--min_cluster_count", str(filter_params["min_cluster_count"]),
        "--min_cluster_samples", str(filter_params["min_cluster_samples"]),
    ]
    
    if filter_params.get("keep_noncanonical", False):
        cmd.append("--keep_noncanonical")
    
    run_command(
        cmd, 
        f"Filtering {cluster_type} clusters",
        skip_if_exists=None if force_rerun else filtered_file
    )
    return filtered_file


def compute_offsets_and_prepare(filtered_file, cluster_type, output_dir, force_rerun=False):
    """
    Compute cluster-total offsets and prepare edgeR inputs.
    
    Args:
        filtered_file: Filtered intron matrix
        cluster_type: 'donor' or 'acceptor'
        output_dir: Output directory
        force_rerun: If True, rerun even if outputs exist
        
    Returns:
        Dict with paths to edgeR input files
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    
    output_prefix = os.path.join(output_dir, f"{cluster_type}_edgeR_input")
    
    # Check if all output files exist
    output_files = {
        "counts": f"{output_prefix}.counts.tsv",
        "offsets": f"{output_prefix}.offsets.tsv",
        "annotations": f"{output_prefix}.annotations.tsv",
    }
    
    all_exist = all(file_exists_and_valid(f) for f in output_files.values())
    
    if not force_rerun and all_exist:
        logger.info(f"=== Computing offsets for {cluster_type} clusters ===")
        logger.info(f"SKIPPING - All output files already exist")
        return output_files
    
    cmd = [
        "python3",
        os.path.join(util_dir, "compute_offsets.py"),
        "--matrix", filtered_file,
        "--output_prefix", output_prefix,
        "--cluster_type", cluster_type,
    ]
    
    run_command(cmd, f"Computing offsets for {cluster_type} clusters")
    
    return output_files


def get_groups_from_samples(samples_file, group_col):
    """
    Extract unique groups from sample metadata file.
    
    Args:
        samples_file: Path to sample metadata file
        group_col: Name of the group column
        
    Returns:
        List of unique group names
    """
    df = pd.read_csv(samples_file, sep="\t", comment="#")
    if group_col not in df.columns:
        raise ValueError(f"Group column '{group_col}' not found in {samples_file}")
    
    # Remove NaN values and convert to string
    groups = df[group_col].dropna().astype(str).unique()
    groups = sorted(groups)
    return groups


def run_single_contrast(contrast, edgeR_inputs, samples_file, cluster_type, output_dir, edgeR_params, force_rerun):
    """
    Run edgeR analysis for a single contrast.
    
    Args:
        contrast: Contrast string (e.g., "GroupA-GroupB")
        edgeR_inputs: Dict with paths to counts, offsets, annotations
        samples_file: Sample metadata file
        cluster_type: 'donor' or 'acceptor'
        output_dir: Output directory
        edgeR_params: Dict with edgeR parameters
        force_rerun: If True, rerun even if outputs exist
        
    Returns:
        Tuple of (contrast, intron_results_file)
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    
    # Create contrast-specific output prefix
    contrast_safe = contrast.replace("-", "_vs_")
    output_prefix = os.path.join(output_dir, f"{cluster_type}_{contrast_safe}_edgeR_results")
    intron_results_file = f"{output_prefix}.intron_results.tsv"
    
    # Check if results already exist
    if not force_rerun and file_exists_and_valid(intron_results_file):
        logger.info(f"SKIPPING contrast {contrast} - Results already exist: {intron_results_file}")
        return (contrast, intron_results_file)
    
    cmd = [
        "Rscript",
        os.path.join(util_dir, "run_edgeR_analysis.R"),
        "--counts", edgeR_inputs["counts"],
        "--offsets", edgeR_inputs["offsets"],
        "--annotations", edgeR_inputs["annotations"],
        "--samples", samples_file,
        "--output", output_prefix,
        "--group_col", edgeR_params["group_col"],
        "--contrast", contrast,
    ]
    
    if edgeR_params.get("batch_col"):
        cmd.extend(["--batch_col", edgeR_params["batch_col"]])
    
    if edgeR_params.get("fdr_threshold"):
        cmd.extend(["--fdr_threshold", str(edgeR_params["fdr_threshold"])])
    
    if edgeR_params.get("min_logFC"):
        cmd.extend(["--min_logFC", str(edgeR_params["min_logFC"])])
    
    logger.info(f"Running contrast: {contrast}")
    run_command(cmd, f"edgeR analysis for {cluster_type} - {contrast}")
    
    return (contrast, intron_results_file)


def run_edgeR(edgeR_inputs, samples_file, cluster_type, output_dir, edgeR_params, force_rerun=False, cpu=1):
    """
    Run edgeR differential analysis.
    
    Args:
        edgeR_inputs: Dict with paths to counts, offsets, annotations
        samples_file: Sample metadata file
        cluster_type: 'donor' or 'acceptor'
        output_dir: Output directory
        edgeR_params: Dict with edgeR parameters
        force_rerun: If True, rerun even if outputs exist
        cpu: Number of CPUs for parallel contrast execution
        
    Returns:
        Path to intron results file
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    
    output_prefix = os.path.join(output_dir, f"{cluster_type}_edgeR_results")
    intron_results_file = f"{output_prefix}.intron_results.tsv"
    
    # Check if combined results already exist
    if not force_rerun and file_exists_and_valid(intron_results_file):
        logger.info(f"=== Running edgeR analysis for {cluster_type} clusters ===")
        logger.info(f"SKIPPING - Results already exist: {intron_results_file}")
        return intron_results_file
    
    # Determine contrasts to run
    if edgeR_params.get("contrast"):
        # Single contrast specified
        logger.info(f"=== Running edgeR analysis for {cluster_type} clusters ===")
        logger.info(f"Single contrast: {edgeR_params['contrast']}")
        
        cmd = [
            "Rscript",
            os.path.join(util_dir, "run_edgeR_analysis.R"),
            "--counts", edgeR_inputs["counts"],
            "--offsets", edgeR_inputs["offsets"],
            "--annotations", edgeR_inputs["annotations"],
            "--samples", samples_file,
            "--output", output_prefix,
            "--group_col", edgeR_params["group_col"],
            "--contrast", edgeR_params["contrast"],
        ]
        
        if edgeR_params.get("batch_col"):
            cmd.extend(["--batch_col", edgeR_params["batch_col"]])
        
        if edgeR_params.get("fdr_threshold"):
            cmd.extend(["--fdr_threshold", str(edgeR_params["fdr_threshold"])])
        
        if edgeR_params.get("min_logFC"):
            cmd.extend(["--min_logFC", str(edgeR_params["min_logFC"])])
        
        run_command(cmd, f"Running edgeR analysis for {cluster_type} clusters")
    else:
        # Multiple contrasts - all pairwise comparisons
        groups = get_groups_from_samples(samples_file, edgeR_params["group_col"])
        contrasts = [f"{g1}-{g2}" for g1, g2 in combinations(groups, 2)]
        
        logger.info(f"=== Running edgeR analysis for {cluster_type} clusters ===")
        logger.info(f"All pairwise comparisons: {len(contrasts)} contrasts among {len(groups)} groups")
        logger.info(f"Groups: {', '.join(groups)}")
        logger.info(f"Using {cpu} CPU(s)")
        
        # Run contrasts in parallel
        run_contrast_partial = partial(
            run_single_contrast,
            edgeR_inputs=edgeR_inputs,
            samples_file=samples_file,
            cluster_type=cluster_type,
            output_dir=output_dir,
            edgeR_params=edgeR_params,
            force_rerun=force_rerun
        )
        
        if cpu > 1:
            logger.info(f"Running {len(contrasts)} contrasts in parallel with {cpu} workers...")
            with Pool(processes=cpu) as pool:
                results = pool.map(run_contrast_partial, contrasts)
        else:
            logger.info(f"Running {len(contrasts)} contrasts sequentially...")
            results = [run_contrast_partial(c) for c in contrasts]
        
        # Combine all results
        logger.info("Combining results from all contrasts...")
        all_dfs = []
        for contrast, result_file in results:
            if file_exists_and_valid(result_file):
                df = pd.read_csv(result_file, sep="\t")
                all_dfs.append(df)
            else:
                logger.warning(f"Results file not found for contrast {contrast}: {result_file}")
        
        if all_dfs:
            combined = pd.concat(all_dfs, ignore_index=True)
            combined.to_csv(intron_results_file, sep="\t", index=False)
            logger.info(f"Combined results written to: {intron_results_file}")
            logger.info(f"Total rows: {len(combined)}")
        else:
            raise RuntimeError("No contrast results were generated")
    
    return intron_results_file


def compute_psi_for_results(edgeR_inputs, samples_file, cluster_type, output_dir, edgeR_params, force_rerun=False):
    """
    Compute PSI values after edgeR analysis.
    
    Args:
        edgeR_inputs: Dict with paths to counts, offsets, annotations
        samples_file: Sample metadata file
        cluster_type: 'donor' or 'acceptor'
        output_dir: Output directory
        edgeR_params: Dict with edgeR parameters (for contrast info)
        force_rerun: If True, rerun even if outputs exist
        
    Returns:
        Path to PSI file
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    
    # Import compute_psi utilities
    sys.path.insert(0, util_dir)
    from compute_psi import compute_psi_values
    
    output_prefix = os.path.join(output_dir, f"{cluster_type}_psi")
    psi_file = f"{output_prefix}.psi_values.tsv"
    
    # Check if PSI file already exists
    if not force_rerun and file_exists_and_valid(psi_file):
        logger.info(f"=== Computing PSI for {cluster_type} clusters ===")
        logger.info(f"SKIPPING - PSI file already exists: {psi_file}")
        return psi_file
    
    logger.info(f"=== Computing PSI for {cluster_type} clusters ===")
    
    try:
        # Load required data
        counts_df = pd.read_csv(edgeR_inputs["counts"], sep="\t", index_col=0)
        annotations_df = pd.read_csv(edgeR_inputs["annotations"], sep="\t", index_col=0)
        samples_df = pd.read_csv(samples_file, sep="\t", comment='#')
        
        # Determine cluster column name
        cluster_col = f"{cluster_type}_cluster"
        
        # Get contrast from parameters
        contrast = edgeR_params.get("contrast")
        
        # Compute PSI values
        psi_df = compute_psi_values(
            counts_df,
            annotations_df,
            samples_df,
            cluster_col=cluster_col,
            group_col=edgeR_params["group_col"],
            contrast=contrast
        )
        
        # Write PSI values
        logger.info(f"Writing PSI values to {psi_file}")
        psi_df.to_csv(psi_file, sep="\t")
        
        logger.info("PSI computation complete!")
        return psi_file
        
    except Exception as e:
        logger.error(f"Error computing PSI: {e}")
        logger.warning("PSI computation failed")
        return None


def add_psi_and_filter(intron_results_file, psi_file, output_dir, cluster_type, min_delta_psi=None, force_rerun=False):
    """
    Add PSI values to edgeR results and optionally filter by delta PSI with FDR recalculation.
    
    Args:
        intron_results_file: Path to intron results from edgeR
        psi_file: Path to PSI values file
        output_dir: Output directory
        cluster_type: 'donor' or 'acceptor'
        min_delta_psi: Minimum absolute delta PSI to include (with FDR recalculation)
        force_rerun: If True, rerun even if outputs exist
        
    Returns:
        Path to PSI-enhanced results file (potentially filtered with recalculated FDR)
    """
    if not psi_file or not file_exists_and_valid(psi_file):
        logger.warning("No PSI file available, skipping PSI annotation")
        return intron_results_file
    
    output_prefix = os.path.join(output_dir, f"{cluster_type}_edgeR_results")
    
    if min_delta_psi:
        psi_enhanced_file = f"{output_prefix}.intron_results_with_psi.psi_filtered.tsv"
    else:
        psi_enhanced_file = f"{output_prefix}.intron_results_with_psi.tsv"
    
    # Check if PSI-enhanced results already exist
    if not force_rerun and file_exists_and_valid(psi_enhanced_file):
        logger.info(f"=== Adding PSI to {cluster_type} results ===")
        logger.info(f"SKIPPING - PSI-enhanced results already exist: {psi_enhanced_file}")
        return psi_enhanced_file
    
    logger.info(f"=== Adding PSI to {cluster_type} results ===")
    
    try:
        # Load results and PSI
        results_df = pd.read_csv(intron_results_file, sep="\t")
        psi_df = pd.read_csv(psi_file, sep="\t", index_col=0)
        
        # Keep only summary PSI columns (not per-sample)
        psi_summary_cols = [col for col in psi_df.columns 
                           if 'mean_PSI' in col or 'median_PSI' in col or 
                              'std_PSI' in col or col == 'delta_PSI']
        
        # Merge PSI data with results
        results_with_index = results_df.set_index('intron_id')
        
        for col in psi_summary_cols:
            if col in psi_df.columns:
                results_with_index[col] = psi_df.loc[results_with_index.index, col]
        
        results_with_psi = results_with_index.reset_index()
        
        # Apply delta PSI filtering and recalculate FDR if threshold specified
        if min_delta_psi and 'delta_PSI' in results_with_psi.columns:
            logger.info(f"Filtering results by |delta_PSI| >= {min_delta_psi} and recalculating FDR")
            
            # Filter by absolute delta PSI
            abs_delta_psi = results_with_psi['delta_PSI'].abs()
            pass_filter = abs_delta_psi >= min_delta_psi
            
            introns_before = len(results_with_psi)
            introns_after = pass_filter.sum()
            
            logger.info(f"Introns before PSI filter: {introns_before}")
            logger.info(f"Introns after PSI filter: {introns_after} ({100*introns_after/introns_before:.1f}%)")
            logger.info(f"Removed {introns_before - introns_after} introns with |delta_PSI| < {min_delta_psi}")
            
            if introns_after == 0:
                logger.warning(f"No introns pass delta_PSI >= {min_delta_psi} filter")
                # Still write the file with all results
                results_with_psi.to_csv(psi_enhanced_file, sep="\t", index=False)
                return psi_enhanced_file
            
            # Filter results
            filtered_results = results_with_psi[pass_filter].copy()
            
            # Recalculate FDR on filtered set using Benjamini-Hochberg
            if 'PValue' in filtered_results.columns:
                from scipy.stats import false_discovery_control
                
                # Get p-values, handling any NaN values
                pvalues = filtered_results['PValue'].values
                valid_pvalues = ~pd.isna(pvalues)
                
                if valid_pvalues.sum() > 0:
                    # Recalculate FDR on filtered subset
                    new_fdr = np.full(len(pvalues), np.nan)
                    new_fdr[valid_pvalues] = false_discovery_control(pvalues[valid_pvalues], method='bh')
                    
                    # Store original FDR for reference
                    filtered_results['FDR_original'] = filtered_results['FDR']
                    filtered_results['FDR'] = new_fdr
                    
                    # Recalculate significance based on new FDR
                    fdr_threshold = 0.05  # Could make this configurable
                    filtered_results['significant'] = filtered_results['FDR'] <= fdr_threshold
                    
                    sig_before = (filtered_results['FDR_original'] <= fdr_threshold).sum()
                    sig_after = (filtered_results['FDR'] <= fdr_threshold).sum()
                    
                    logger.info(f"Significant introns before FDR recalculation: {sig_before}")
                    logger.info(f"Significant introns after FDR recalculation: {sig_after}")
                    logger.info(f"Gained {sig_after - sig_before} significant introns from reduced multiple testing burden")
            
            results_with_psi = filtered_results
        
        # Write PSI-enhanced results
        logger.info(f"Writing PSI-enhanced results to {psi_enhanced_file}")
        results_with_psi.to_csv(psi_enhanced_file, sep="\t", index=False)
        
        logger.info(f"Added {len(psi_summary_cols)} PSI columns to results")
        return psi_enhanced_file
        
    except Exception as e:
        logger.error(f"Error adding PSI to results: {e}")
        logger.warning("Failed to add PSI, will use original results")
        return intron_results_file


def aggregate_results(intron_results_file, cluster_type, output_dir, agg_params, force_rerun=False):
    """
    Aggregate intron-level results to cluster level.
    
    Args:
        intron_results_file: Path to intron results from edgeR
        cluster_type: 'donor' or 'acceptor'
        output_dir: Output directory
        agg_params: Dict with aggregation parameters
        force_rerun: If True, rerun even if outputs exist
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    
    output_prefix = os.path.join(output_dir, f"{cluster_type}_aggregated")
    cluster_results_file = f"{output_prefix}.cluster_results.tsv"
    
    # Ensure parameters are properly set with defaults
    method = agg_params.get("method") or "cauchy"
    cluster_fdr = agg_params.get("cluster_fdr")
    if cluster_fdr is None:
        cluster_fdr = 0.05
    
    cmd = [
        "python3",
        os.path.join(util_dir, "aggregate_clusters.py"),
        "--intron_results", intron_results_file,
        "--output_prefix", output_prefix,
        "--method", method,
        "--cluster_fdr", str(cluster_fdr),
    ]
    
    run_command(
        cmd, 
        f"Aggregating {cluster_type} cluster results",
        skip_if_exists=None if force_rerun else cluster_results_file
    )


def integrate_donor_acceptor_results(output_dir, cluster_types, edgeR_params, gtf=None, output_prefix="integrated", force_rerun=False):
    """
    Integrate results from donor and acceptor analyses.
    
    Args:
        output_dir: Output directory
        cluster_types: List of cluster types that were run
        edgeR_params: Dict with edgeR parameters
        gtf: Optional GTF file for gene annotation
        output_prefix: Prefix for integrated output files
        force_rerun: If True, rerun even if outputs exist
    """
    # Only integrate if both donor and acceptor were run
    if 'donor' not in cluster_types or 'acceptor' not in cluster_types:
        logger.info("Skipping integration - both donor and acceptor analyses required")
        return
    
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    
    # Try to use PSI-filtered results if available, otherwise PSI-enhanced, then fall back to original
    donor_results = os.path.join(output_dir, "donor", "donor_edgeR_results.intron_results_with_psi.psi_filtered.tsv")
    if not file_exists_and_valid(donor_results):
        donor_results = os.path.join(output_dir, "donor", "donor_edgeR_results.intron_results_with_psi.tsv")
    if not file_exists_and_valid(donor_results):
        donor_results = os.path.join(output_dir, "donor", "donor_edgeR_results.intron_results.tsv")
    
    acceptor_results = os.path.join(output_dir, "acceptor", "acceptor_edgeR_results.intron_results_with_psi.psi_filtered.tsv")
    if not file_exists_and_valid(acceptor_results):
        acceptor_results = os.path.join(output_dir, "acceptor", "acceptor_edgeR_results.intron_results_with_psi.tsv")
    if not file_exists_and_valid(acceptor_results):
        acceptor_results = os.path.join(output_dir, "acceptor", "acceptor_edgeR_results.intron_results.tsv")
    
    # Check if both input files exist
    if not (file_exists_and_valid(donor_results) and file_exists_and_valid(acceptor_results)):
        logger.warning("Cannot integrate results - missing donor or acceptor intron results")
        return
    
    full_output_prefix = os.path.join(output_dir, output_prefix)
    integrated_file = f"{full_output_prefix}.integrated_results.tsv"
    
    cmd = [
        "python3",
        os.path.join(util_dir, "integrate_results.py"),
        "--donor_results", donor_results,
        "--acceptor_results", acceptor_results,
        "--output_prefix", full_output_prefix,
        "--fdr_threshold", str(edgeR_params.get("fdr_threshold", 0.05)),
    ]
    
    if gtf:
        cmd.extend(["--gtf", gtf])
    
    run_command(
        cmd,
        "Integrating donor and acceptor results",
        skip_if_exists=None if force_rerun else integrated_file
    )


def main():
    parser = argparse.ArgumentParser(
        description="Run full differential splicing analysis pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    # Input files
    parser.add_argument(
        "--matrix",
        type=str,
        required=True,
        help="Input intron count matrix (supports .tsv or .tsv.gz)",
    )
    
    parser.add_argument(
        "--samples",
        type=str,
        required=True,
        help="Sample metadata file (TSV with columns: sample_id, group, [batch])",
    )
    
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory for results",
    )
    
    parser.add_argument(
        "--gtf",
        type=str,
        default=None,
        help="GTF annotation file for gene annotation and known/novel intron status (optional)",
    )
    
    parser.add_argument(
        "--output-prefix",
        type=str,
        default="integrated",
        help="Prefix for integrated output files (default: integrated)",
    )
    
    # Clustering options
    parser.add_argument(
        "--cluster_types",
        type=str,
        nargs="+",
        default=["donor", "acceptor"],
        choices=["donor", "acceptor"],
        help="Types of clustering to perform",
    )
    
    # Filtering parameters
    parser.add_argument(
        "--min_intron_count",
        type=int,
        default=10,
        help="Minimum total count for an intron",
    )
    
    parser.add_argument(
        "--min_intron_samples",
        type=int,
        default=2,
        help="Minimum samples with non-zero counts for an intron",
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
        help="Minimum samples meeting cluster count threshold",
    )
    
    parser.add_argument(
        "--keep_noncanonical",
        action="store_true",
        help="Keep non-canonical splice sites",
    )
    
    # edgeR parameters
    parser.add_argument(
        "--group_col",
        type=str,
        default="group",
        help="Column name for sample groups in metadata file",
    )
    
    parser.add_argument(
        "--batch_col",
        type=str,
        default=None,
        help="Optional column name for batch effects",
    )
    
    parser.add_argument(
        "--contrast",
        type=str,
        default=None,
        help="Contrast to test (e.g., 'TDP43-control')",
    )
    
    parser.add_argument(
        "--fdr_threshold",
        type=float,
        default=0.05,
        help="FDR threshold for significance",
    )
    
    parser.add_argument(
        "--min_logFC",
        type=float,
        default=0.0,
        help="Minimum absolute log2FC for significance",
    )
    
    # PSI filtering parameters
    parser.add_argument(
        "--min_delta_psi",
        type=float,
        default=None,
        help="Minimum absolute delta PSI to include in final results. FDR will be recalculated on the filtered set (reduces multiple testing burden). Example: 0.1 for 10%% change",
    )
    
    # Aggregation parameters
    parser.add_argument(
        "--agg_method",
        type=str,
        default="cauchy",
        choices=["fisher", "cauchy"],
        help="P-value combination method for cluster aggregation",
    )
    
    parser.add_argument(
        "--cluster_fdr",
        type=float,
        default=0.05,
        help="FDR threshold for cluster-level significance",
    )
    
    # Pipeline control
    parser.add_argument(
        "--cpu",
        type=int,
        default=1,
        help="Number of CPU threads for parallel contrast testing (default: 1)",
    )
    
    parser.add_argument(
        "--force_rerun",
        action="store_true",
        help="Force rerun of all steps even if output files exist (disables resume)",
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    logger.info("=== Differential Splicing Analysis Pipeline ===")
    logger.info(f"Input matrix: {args.matrix}")
    logger.info(f"Sample metadata: {args.samples}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Cluster types: {args.cluster_types}")
    
    if args.min_delta_psi:
        logger.info(f"PSI filtering: |delta_PSI| >= {args.min_delta_psi} (with FDR recalculation on filtered set)")
    else:
        logger.info("PSI filtering: disabled (all introns will be included)")
    
    if args.force_rerun:
        logger.info("Force rerun enabled - will regenerate all outputs")
    else:
        logger.info("Resume mode enabled - will skip completed steps")
    
    # Prepare parameter dicts
    filter_params = {
        "min_intron_count": args.min_intron_count,
        "min_intron_samples": args.min_intron_samples,
        "min_cluster_count": args.min_cluster_count,
        "min_cluster_samples": args.min_cluster_samples,
        "keep_noncanonical": args.keep_noncanonical,
    }
    
    edgeR_params = {
        "group_col": args.group_col,
        "batch_col": args.batch_col,
        "contrast": args.contrast,
        "fdr_threshold": args.fdr_threshold,
        "min_logFC": args.min_logFC,
    }
    
    agg_params = {
        "method": args.agg_method,
        "cluster_fdr": args.cluster_fdr,
    }
    
    # Run analysis for each cluster type
    for cluster_type in args.cluster_types:
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing {cluster_type.upper()} clusters")
        logger.info(f"{'='*60}\n")
        
        # Create subdirectory for this cluster type
        cluster_dir = os.path.join(args.output_dir, cluster_type)
        os.makedirs(cluster_dir, exist_ok=True)
        
        try:
            # 1. Cluster and filter
            filtered_file = cluster_and_filter(
                args.matrix, cluster_type, cluster_dir, filter_params, 
                force_rerun=args.force_rerun
            )
            
            # 2. Compute offsets and prepare edgeR inputs
            edgeR_inputs = compute_offsets_and_prepare(
                filtered_file, cluster_type, cluster_dir,
                force_rerun=args.force_rerun
            )
            
            # 3. Run edgeR (on all data for proper dispersion estimation)
            intron_results = run_edgeR(
                edgeR_inputs, args.samples, cluster_type, cluster_dir, edgeR_params,
                force_rerun=args.force_rerun,
                cpu=args.cpu
            )
            
            # 4. Compute PSI values
            psi_file = compute_psi_for_results(
                edgeR_inputs, args.samples, cluster_type, cluster_dir,
                edgeR_params, force_rerun=args.force_rerun
            )
            
            # 5. Add PSI to results and optionally filter by delta PSI with FDR recalculation
            intron_results_with_psi = add_psi_and_filter(
                intron_results, psi_file, cluster_dir, cluster_type,
                min_delta_psi=args.min_delta_psi, force_rerun=args.force_rerun
            )
            
            # 6. Aggregate to cluster level (use PSI-enhanced results if available)
            aggregate_results(
                intron_results_with_psi, cluster_type, cluster_dir, agg_params,
                force_rerun=args.force_rerun
            )
            
            logger.info(f"\n{cluster_type.upper()} cluster analysis complete!")
            logger.info(f"Results in: {cluster_dir}")
            
        except Exception as e:
            logger.error(f"Error processing {cluster_type} clusters: {e}")
            raise
    
    # Integrate donor and acceptor results if both were run
    if len(args.cluster_types) > 1:
        logger.info(f"\n{'='*60}")
        logger.info("Integrating Results Across Clustering Strategies")
        logger.info(f"{'='*60}\n")
        
        try:
            integrate_donor_acceptor_results(
                args.output_dir,
                args.cluster_types,
                edgeR_params,
                gtf=args.gtf,
                output_prefix=getattr(args, 'output_prefix', 'integrated'),
                force_rerun=args.force_rerun
            )
        except Exception as e:
            logger.error(f"Error integrating results: {e}")
            logger.warning("Integration failed, but individual analyses completed successfully")
    
    logger.info("\n" + "="*60)
    logger.info("PIPELINE COMPLETE!")
    logger.info("="*60)
    logger.info(f"\nAll results saved to: {args.output_dir}")
    
    # Print summary of key output files
    logger.info("\nKey output files:")
    for cluster_type in args.cluster_types:
        cluster_dir = os.path.join(args.output_dir, cluster_type)
        logger.info(f"\n{cluster_type.upper()} clusters:")
        logger.info(f"  - Intron results: {cluster_dir}/{cluster_type}_edgeR_results.intron_results.tsv")
        logger.info(f"  - Cluster results: {cluster_dir}/{cluster_type}_aggregated.cluster_results.tsv")
        logger.info(f"  - Diagnostics: {cluster_dir}/{cluster_type}_edgeR_results.diagnostics.pdf")
    
    # Print integrated results if both analyses were run
    if len(args.cluster_types) > 1:
        prefix = getattr(args, 'output_prefix', 'integrated')
        integrated_file = os.path.join(args.output_dir, f"{prefix}.integrated_results.tsv")
        if file_exists_and_valid(integrated_file):
            logger.info(f"\nINTEGRATED results (combining donor + acceptor):")
            logger.info(f"  - All introns: {args.output_dir}/{prefix}.integrated_results.tsv")
            logger.info(f"  - Significant only: {args.output_dir}/{prefix}.significant_integrated.tsv")
            logger.info(f"  - Summary stats: {args.output_dir}/{prefix}.integration_summary.tsv")


if __name__ == "__main__":
    main()
