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
    
    donor_results = os.path.join(output_dir, "donor", "donor_edgeR_results.intron_results.tsv")
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
            
            # 3. Run edgeR
            intron_results = run_edgeR(
                edgeR_inputs, args.samples, cluster_type, cluster_dir, edgeR_params,
                force_rerun=args.force_rerun,
                cpu=args.cpu
            )
            
            # 4. Aggregate to cluster level
            aggregate_results(
                intron_results, cluster_type, cluster_dir, agg_params,
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
