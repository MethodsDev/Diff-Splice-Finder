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
import numpy as np
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


def annotate_clustered_file_with_genes(clustered_file, gtf_file, force_rerun=False):
    """
    Add gene annotation columns to clustered file.
    
    Args:
        clustered_file: Path to clustered intron file
        gtf_file: Path to GTF annotation file  
        force_rerun: If True, rerun even if annotated file exists
        
    Returns:
        Path to annotated clustered file
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    sys.path.insert(0, util_dir)
    from integrate_results import parse_gtf_file, parse_intron_id, build_gene_index, find_overlapping_genes
    
    annotated_file = clustered_file.replace('.tsv', '.annotated.tsv')
    
    if not force_rerun and file_exists_and_valid(annotated_file):
        logger.info(f"Using existing annotated file: {annotated_file}")
        return annotated_file
        
    logger.info(f"Annotating {clustered_file} with gene information...")
    
    # Parse GTF
    gene_map, annotated_introns, transcript_genes, transcript_exons, intron_gene_map = parse_gtf_file(gtf_file)
    gene_index = build_gene_index(gene_map, transcript_genes, transcript_exons)
    
    # Read clustered file
    df = pd.read_csv(clustered_file, sep="\t", index_col=0)
    logger.info(f"Loaded clustered file with {len(df)} introns and columns: {list(df.columns)}")
    
    # Annotate each intron
    gene_names = []
    intron_statuses = []
    
    for intron_id in df.index:
        coords = parse_intron_id(intron_id)
        
        if coords is None:
            gene_names.append('.')
            intron_statuses.append('unknown')
            continue
        
        # Check if known
        is_known = coords in annotated_introns
        intron_statuses.append('known' if is_known else 'novel')
        
        # Find overlapping genes
        overlapping = find_overlapping_genes(coords, gene_index, intron_gene_map)
        
        if overlapping:
            gene_names.append(overlapping[0])
        else:
            gene_names.append('.')
    
    # Add columns at beginning (after intron_id)
    df.insert(0, 'gene_name', gene_names)
    df.insert(1, 'intron_status', intron_statuses)
    
    # Write annotated file
    df.to_csv(annotated_file, sep="\t", na_rep='NA')
    logger.info(f"Wrote annotated file: {annotated_file}")
    
    # Log stats
    known = sum(1 for s in intron_statuses if s == 'known')
    novel = sum(1 for s in intron_statuses if s == 'novel')
    with_genes = sum(1 for g in gene_names if g != '.')
    logger.info(f"  Known: {known}, Novel: {novel}, With genes: {with_genes}")
    
    return annotated_file


def cluster_and_filter(matrix_file, cluster_type, output_dir, filter_params, force_rerun=False, annotated_file=None):
    """
    Cluster introns and apply filtering.
    
    Args:
        matrix_file: Input intron count matrix
        cluster_type: 'donor' or 'acceptor'
        output_dir: Output directory
        filter_params: Dict with filtering parameters
        force_rerun: If True, rerun even if outputs exist
        annotated_file: If provided, filter this annotated file instead of the original clustered file
        
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
    
    # Step 2: Filter (use annotated file if provided, otherwise use clustered file)
    input_for_filtering = annotated_file if annotated_file else clustered_file
    filtered_file = os.path.join(output_dir, f"{cluster_type}_filtered.tsv")
    
    cmd = [
        "python3",
        os.path.join(util_dir, "filter_introns.py"),
        "--matrix", input_for_filtering,
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


def compute_shared_offsets(annotated_clustered_file, output_dir, force_rerun=False):
    """
    Compute shared offsets from the full clustered matrix.
    These offsets will be used by both donor and acceptor analyses.
    
    Args:
        annotated_clustered_file: Annotated clustered matrix with both donor and acceptor clusters
        output_dir: Output directory
        force_rerun: If True, rerun even if output exists
        
    Returns:
        Path to shared offsets file
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    
    shared_offsets_file = os.path.join(output_dir, "shared_offsets.tsv")
    
    if not force_rerun and file_exists_and_valid(shared_offsets_file):
        logger.info("=== Computing shared offsets ===")
        logger.info(f"SKIPPING - Shared offsets file already exists")
        return shared_offsets_file
    
    cmd = [
        "python3",
        os.path.join(util_dir, "compute_offsets.py"),
        "--matrix", annotated_clustered_file,
        "--output", shared_offsets_file,
        "--compute_offsets_only",  # New flag to only compute offsets
    ]
    
    run_command(cmd, "Computing shared offsets from full clustered matrix")
    
    return shared_offsets_file


def prepare_edgeR_inputs(filtered_file, output_dir, shared_offsets_file, samples_file=None, force_rerun=False):
    """
    Prepare edgeR input files using pre-computed shared offsets.
    
    Args:
        filtered_file: Filtered intron matrix
        output_dir: Output directory
        shared_offsets_file: Path to shared offsets file
        samples_file: Path to sample metadata file (used to filter samples)
        force_rerun: If True, rerun even if outputs exist
        
    Returns:
        Dict with paths to edgeR input files
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    
    output_prefix = os.path.join(output_dir, "edgeR_input")
    
    # Check if all output files exist
    output_files = {
        "counts": f"{output_prefix}.counts.tsv",
        "offsets": f"{output_prefix}.offsets.tsv",
        "annotations": f"{output_prefix}.annotations.tsv",
    }
    
    all_exist = all(file_exists_and_valid(f) for f in output_files.values())
    
    if not force_rerun and all_exist:
        logger.info(f"=== Preparing edgeR inputs ===")
        logger.info(f"SKIPPING - All output files already exist")
        return output_files
    
    cmd = [
        "python3",
        os.path.join(util_dir, "compute_offsets.py"),
        "--matrix", filtered_file,
        "--output_prefix", output_prefix,
        "--shared_offsets", shared_offsets_file,
    ]
    
    # Add sample filtering if metadata file provided
    if samples_file:
        cmd.extend(["--samples", samples_file])
    
    run_command(cmd, f"Preparing edgeR inputs")
    
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


def run_single_contrast(contrast, edgeR_inputs, samples_file, output_dir, edgeR_params, force_rerun):
    """
    Run edgeR analysis for a single contrast.
    
    Args:
        contrast: Contrast string (e.g., "GroupA-GroupB")
        edgeR_inputs: Dict with paths to counts, offsets, annotations
        samples_file: Sample metadata file
        output_dir: Output directory
        edgeR_params: Dict with edgeR parameters
        force_rerun: If True, rerun even if outputs exist
        
    Returns:
        Tuple of (contrast, intron_results_file)
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    
    # Create contrast-specific output prefix
    contrast_safe = contrast.replace("-", "_vs_")
    output_prefix = os.path.join(output_dir, f"{contrast_safe}_edgeR_results")
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
    run_command(cmd, f"edgeR analysis - {contrast}")
    
    return (contrast, intron_results_file)


def run_edgeR(edgeR_inputs, samples_file, output_dir, edgeR_params, force_rerun=False, cpu=1):
    """
    Run edgeR differential analysis.
    
    Args:
        edgeR_inputs: Dict with paths to counts, offsets, annotations
        samples_file: Sample metadata file
        output_dir: Output directory
        edgeR_params: Dict with edgeR parameters
        force_rerun: If True, rerun even if outputs exist
        cpu: Number of CPUs for parallel contrast execution
        
    Returns:
        Path to intron results file
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    
    output_prefix = os.path.join(output_dir, "edgeR_results")
    intron_results_file = f"{output_prefix}.intron_results.tsv"
    
    # Check if combined results already exist
    if not force_rerun and file_exists_and_valid(intron_results_file):
        logger.info("=== Running edgeR analysis ===")
        logger.info(f"SKIPPING - Results already exist: {intron_results_file}")
        return intron_results_file
    
    # Determine contrasts to run
    if edgeR_params.get("contrast"):
        # Single contrast specified
        logger.info("=== Running edgeR analysis ===")
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
        
        run_command(cmd, "Running edgeR analysis")
    else:
        # Multiple contrasts - either control-based or all pairwise comparisons
        groups = get_groups_from_samples(samples_file, edgeR_params["group_col"])
        
        # Check if control groups are specified
        if edgeR_params.get("control_groups"):
            # Parse control groups
            control_groups = [g.strip() for g in edgeR_params["control_groups"].split(",")]
            
            # Validate control groups exist
            missing_controls = [g for g in control_groups if g not in groups]
            if missing_controls:
                raise ValueError(f"Control groups not found in data: {', '.join(missing_controls)}. Available groups: {', '.join(groups)}")
            
            # Generate contrasts: all non-control groups vs control
            non_control_groups = [g for g in groups if g not in control_groups]
            
            if not non_control_groups:
                raise ValueError("No non-control groups found. All groups were specified as controls.")
            
            # For each non-control group, compare against controls
            # Format: TreatmentGroup-ControlGroup1,ControlGroup2
            if len(control_groups) == 1:
                contrasts = [f"{g}-{control_groups[0]}" for g in non_control_groups]
            else:
                control_str = ",".join(control_groups)
                contrasts = [f"{g}-{control_str}" for g in non_control_groups]
            
            logger.info("=== Running edgeR analysis ===")
            logger.info(f"Control-based comparisons: {len(contrasts)} contrasts")
            logger.info(f"Control groups: {', '.join(control_groups)}")
            logger.info(f"Treatment groups: {', '.join(non_control_groups)}")
            logger.info(f"Using {cpu} CPU(s)")
        else:
            # Original behavior: all pairwise comparisons
            contrasts = [f"{g1}-{g2}" for g1, g2 in combinations(groups, 2)]
            
            logger.info("=== Running edgeR analysis ===")
            logger.info(f"All pairwise comparisons: {len(contrasts)} contrasts among {len(groups)} groups")
            logger.info(f"Groups: {', '.join(groups)}")
            logger.info(f"Using {cpu} CPU(s)")
        
        # Run contrasts in parallel
        run_contrast_partial = partial(
            run_single_contrast,
            edgeR_inputs=edgeR_inputs,
            samples_file=samples_file,
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
            combined.to_csv(intron_results_file, sep="\t", index=False, na_rep='NA')
            logger.info(f"Combined results written to: {intron_results_file}")
            logger.info(f"Total rows: {len(combined)}")
        else:
            raise RuntimeError("No contrast results were generated")
    
    return intron_results_file


def compute_psi_for_results(edgeR_inputs, samples_file, output_dir, edgeR_params, shared_offsets_file=None, force_rerun=False):
    """
    Compute PSI values after edgeR analysis.
    
    Args:
        edgeR_inputs: Dict with paths to counts, offsets, annotations
        samples_file: Sample metadata file
        output_dir: Output directory
        edgeR_params: Dict with edgeR parameters (for contrast info)
        shared_offsets_file: Path to shared offsets file (raw cluster totals)
        force_rerun: If True, rerun even if outputs exist
        
    Returns:
        Path to PSI file
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    
    # Import compute_psi utilities
    sys.path.insert(0, util_dir)
    from compute_psi import compute_psi_values
    
    output_prefix = os.path.join(output_dir, "psi")
    psi_file = f"{output_prefix}.psi_values.tsv"
    
    # Check if PSI file already exists
    if not force_rerun and file_exists_and_valid(psi_file):
        logger.info("=== Computing PSI ===")
        logger.info(f"SKIPPING - PSI file already exists: {psi_file}")
        return psi_file
    
    logger.info("=== Computing PSI ===")
    
    try:
        # Load required data
        counts_df = pd.read_csv(edgeR_inputs["counts"], sep="\t", index_col=0)
        annotations_df = pd.read_csv(edgeR_inputs["annotations"], sep="\t", index_col=0)
        samples_df = pd.read_csv(samples_file, sep="\t", comment='#')
        
        # Load shared offsets if provided (raw cluster totals, not log-transformed)
        shared_cluster_totals = None
        if shared_offsets_file:
            logger.info(f"Using shared offsets for consistent PSI denominators: {shared_offsets_file}")
            shared_cluster_totals = pd.read_csv(shared_offsets_file, sep="\t", index_col=0)
        
        # Determine cluster column name (check which one is available)
        # When using shared offsets, the cluster_col is not actually used, but we still need to provide it
        if 'donor_cluster' in annotations_df.columns:
            cluster_col = 'donor_cluster'
        elif 'acceptor_cluster' in annotations_df.columns:
            cluster_col = 'acceptor_cluster'
        else:
            raise ValueError("No cluster column found in annotations (expected 'donor_cluster' or 'acceptor_cluster')")
        
        # Get contrast from parameters
        contrast = edgeR_params.get("contrast")
        
        # Compute PSI values
        psi_df = compute_psi_values(
            counts_df,
            annotations_df,
            samples_df,
            cluster_col=cluster_col,
            group_col=edgeR_params["group_col"],
            shared_cluster_totals=shared_cluster_totals,
            contrast=contrast
        )
        
        # Write PSI values
        logger.info(f"Writing PSI values to {psi_file}")
        psi_df.to_csv(psi_file, sep="\t", na_rep='NA')
        
        logger.info("PSI computation complete!")
        return psi_file
        
    except Exception as e:
        logger.error(f"Error computing PSI: {e}")
        logger.warning("PSI computation failed")
        return None


def add_psi_and_filter(intron_results_file, psi_file, output_dir, min_delta_psi=None, force_rerun=False):
    """
    Add PSI values to edgeR results and optionally filter by delta PSI with FDR recalculation.
    
    Always creates an unfiltered file with PSI values. If min_delta_psi is specified, 
    also creates a filtered version.
    
    Args:
        intron_results_file: Path to intron results from edgeR
        psi_file: Path to PSI values file
        output_dir: Output directory
        min_delta_psi: Minimum absolute delta PSI to include (with FDR recalculation)
        force_rerun: If True, rerun even if outputs exist
        
    Returns:
        Path to final results file (filtered if threshold specified, otherwise unfiltered)
    """
    if not psi_file or not file_exists_and_valid(psi_file):
        logger.warning("No PSI file available, skipping PSI annotation")
        return intron_results_file
    
    output_prefix = os.path.join(output_dir, "edgeR_results")
    
    # Always create unfiltered file
    unfiltered_file = f"{output_prefix}.intron_results_with_psi.tsv"
    
    # Determine final output file
    if min_delta_psi:
        final_file = f"{output_prefix}.intron_results_with_psi.psi_filtered.tsv"
    else:
        final_file = unfiltered_file
    
    # Check if final output already exists (skip all if so)
    if not force_rerun and file_exists_and_valid(final_file):
        logger.info("=== Adding PSI to results ===")
        logger.info(f"SKIPPING - Results already exist: {final_file}")
        return final_file
    
    logger.info("=== Adding PSI to results ===")
    
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
        
        # Write unfiltered PSI-enhanced results first
        logger.info(f"Writing unfiltered results with PSI to {unfiltered_file}")
        results_with_psi.to_csv(unfiltered_file, sep="\t", index=False, na_rep='NA')
        logger.info(f"Added {len(psi_summary_cols)} PSI columns to results")
        
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
                logger.warning(f"Filtered file will be empty but unfiltered file is available: {unfiltered_file}")
                # Create empty filtered file
                pd.DataFrame(columns=results_with_psi.columns).to_csv(final_file, sep="\t", index=False)
                return final_file
            
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
            
            # Write filtered results
            logger.info(f"Writing PSI-filtered results to {final_file}")
            filtered_results.to_csv(final_file, sep="\t", index=False, na_rep='NA')
        
        return final_file
        
    except Exception as e:
        logger.error(f"Error adding PSI to results: {e}")
        logger.warning("Failed to add PSI, will use original results")
        return intron_results_file




# ============================================================================
# DEPRECATED FUNCTIONS
# ============================================================================
# The following functions (aggregate_results and integrate_donor_acceptor_results)
# are no longer used in the intron-level analysis approach but are kept for reference.
# 
# Previous approach: Ran separate donor and acceptor analyses, aggregated cluster-level
# results, then integrated the two analyses.
# 
# Current approach: Single intron-level analysis with shared offsets computed from
# max(donor_cluster_total, acceptor_cluster_total), eliminating need for aggregation
# and integration steps.
# ============================================================================

def aggregate_results(intron_results_file, cluster_type, output_dir, agg_params, force_rerun=False):
    """
    DEPRECATED: Aggregate intron-level results to cluster level.
    
    This function is no longer used in the intron-level analysis approach.
    
    Args:
        intron_results_file: Path to intron results from edgeR
        cluster_type: 'donor' or 'acceptor'
        output_dir: Output directory
        agg_params: Dict with aggregation parameters
        force_rerun: If True, rerun even if outputs exist
    """
    raise NotImplementedError("Cluster aggregation has been removed. The pipeline now uses intron-level analysis only.")


def integrate_donor_acceptor_results(output_dir, cluster_types, edgeR_params, gtf=None, output_prefix="integrated", force_rerun=False):
    """
    DEPRECATED: Integrate results from donor and acceptor analyses.
    
    This function is no longer used in the intron-level analysis approach.
    
    Args:
        output_dir: Output directory
        cluster_types: List of cluster types that were run
        edgeR_params: Dict with edgeR parameters
        gtf: Optional GTF file for gene annotation
        output_prefix: Prefix for integrated output files
        force_rerun: If True, rerun even if outputs exist
    """
    raise NotImplementedError("Integration of donor/acceptor results has been removed. The pipeline now uses a single intron-level analysis.")


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
        "--control_groups",
        type=str,
        default=None,
        help="Comma-separated list of control group names to compare all other groups against (e.g., 'control,wildtype'). If specified, all non-control groups will be compared to the pooled control groups instead of all pairwise comparisons.",
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
    
    # Validate mutually exclusive options
    if args.contrast and args.control_groups:
        parser.error("Cannot use --contrast and --control_groups together. Use one or the other.")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    logger.info("=== Differential Splicing Analysis Pipeline ===")
    logger.info(f"Input matrix: {args.matrix}")
    logger.info(f"Sample metadata: {args.samples}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info("Analysis mode: Intron-level with shared offsets")
    
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
        "control_groups": args.control_groups,
        "fdr_threshold": args.fdr_threshold,
        "min_logFC": args.min_logFC,
    }
    
    # Step 1: Cluster introns with BOTH donor and acceptor (done once)
    logger.info(f"\n{'='*60}")
    logger.info("Clustering introns by BOTH donor and acceptor")
    logger.info(f"{'='*60}\n")
    
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    clustered_file = os.path.join(args.output_dir, "introns_clustered.tsv")
    
    cmd = [
        "python3",
        os.path.join(util_dir, "cluster_introns.py"),
        "--matrix", args.matrix,
        "--output_donor", clustered_file,  # Will write after both donor and acceptor clustering
        "--cluster_type", "both",
    ]
    run_command(
        cmd,
        "Clustering introns by donor and acceptor",
        skip_if_exists=None if args.force_rerun else clustered_file
    )
    
    # Step 2: Annotate clustered file with genes (if GTF provided)
    annotated_clustered = None
    if args.gtf:
        annotated_clustered = annotate_clustered_file_with_genes(
            clustered_file, args.gtf, force_rerun=args.force_rerun
        )
        logger.info(f"Gene-annotated clustered file: {annotated_clustered}")
    else:
        annotated_clustered = clustered_file
    
    # Step 3: Compute shared offsets from full clustered matrix
    logger.info(f"\n{'='*60}")
    logger.info("Computing shared offsets")
    logger.info(f"{'='*60}\n")
    
    shared_offsets_file = compute_shared_offsets(
        annotated_clustered, args.output_dir, force_rerun=args.force_rerun
    )
    logger.info(f"Shared offsets file: {shared_offsets_file}")
    
    # Step 4: Filter introns (require thresholds for both donor and acceptor clusters)
    logger.info(f"\n{'='*60}")
    logger.info("Filtering introns")
    logger.info(f"{'='*60}\n")
    
    filtered_file = os.path.join(args.output_dir, "introns_filtered.tsv")
    
    cmd = [
        "python3",
        os.path.join(util_dir, "filter_introns.py"),
        "--matrix", annotated_clustered,
        "--output", filtered_file,
        "--min_intron_count", str(filter_params["min_intron_count"]),
        "--min_intron_samples", str(filter_params["min_intron_samples"]),
        "--min_cluster_count", str(filter_params["min_cluster_count"]),
        "--min_cluster_samples", str(filter_params["min_cluster_samples"]),
    ]
    if filter_params.get("keep_noncanonical", False):
        cmd.append("--keep_noncanonical")
    
    run_command(
        cmd,
        "Filtering introns",
        skip_if_exists=None if args.force_rerun else filtered_file
    )
    
    # Step 5: Prepare edgeR inputs using shared offsets
    logger.info(f"\n{'='*60}")
    logger.info("Preparing edgeR inputs")
    logger.info(f"{'='*60}\n")
    
    edgeR_inputs = prepare_edgeR_inputs(
        filtered_file, args.output_dir,
        shared_offsets_file,
        samples_file=args.samples,
        force_rerun=args.force_rerun
    )
    
    # Step 6: Run edgeR analysis
    logger.info(f"\n{'='*60}")
    logger.info("Running edgeR analysis")
    logger.info(f"{'='*60}\n")
    
    intron_results = run_edgeR(
        edgeR_inputs, args.samples, args.output_dir, edgeR_params,
        force_rerun=args.force_rerun,
        cpu=args.cpu
    )
    
    # Step 7: Compute PSI values using shared offsets
    logger.info(f"\n{'='*60}")
    logger.info("Computing PSI values")
    logger.info(f"{'='*60}\n")
    
    psi_file = compute_psi_for_results(
        edgeR_inputs, args.samples, args.output_dir,
        edgeR_params, 
        shared_offsets_file=shared_offsets_file,
        force_rerun=args.force_rerun
    )
    
    # Step 8: Add PSI to results and optionally filter by delta PSI
    logger.info(f"\n{'='*60}")
    logger.info("Adding PSI to results")
    logger.info(f"{'='*60}\n")
    
    intron_results_with_psi = add_psi_and_filter(
        intron_results, psi_file, args.output_dir,
        min_delta_psi=args.min_delta_psi, force_rerun=args.force_rerun
    )
    
    logger.info("\n" + "="*60)
    logger.info("PIPELINE COMPLETE!")
    logger.info("="*60)
    logger.info(f"\nAll results saved to: {args.output_dir}")
    
    # Print summary of key output files
    logger.info("\nKey output files:")
    logger.info(f"  - Intron results: {args.output_dir}/edgeR_results.intron_results.tsv")
    
    # Always show the unfiltered PSI file
    unfiltered_psi_file = f"{args.output_dir}/edgeR_results.intron_results_with_psi.tsv"
    if file_exists_and_valid(unfiltered_psi_file):
        logger.info(f"  - Results with PSI (unfiltered): {unfiltered_psi_file}")
    
    # Show filtered PSI file if delta PSI threshold was used
    if args.min_delta_psi and intron_results_with_psi and file_exists_and_valid(intron_results_with_psi):
        logger.info(f"  - Results with PSI (filtered |delta_PSI| >= {args.min_delta_psi}): {intron_results_with_psi}")
    
    logger.info(f"  - PSI values: {args.output_dir}/psi.psi_values.tsv")
    logger.info(f"  - Diagnostics: {args.output_dir}/edgeR_results.diagnostics.pdf")


if __name__ == "__main__":
    main()
