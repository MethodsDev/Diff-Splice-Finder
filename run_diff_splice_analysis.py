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
import glob
import shlex
import re
import pandas as pd
import numpy as np
from multiprocessing import Pool
from itertools import combinations
from functools import partial

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
    stream=sys.stdout,
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


def file_is_current(output_path, input_paths):
    """
    Return True if output_path exists and is newer than all existing input_paths.
    """
    if not file_exists_and_valid(output_path):
        return False

    output_mtime = os.path.getmtime(output_path)
    for input_path in input_paths:
        if input_path and os.path.exists(input_path) and os.path.getmtime(input_path) > output_mtime:
            return False
    return True


def get_bam_paths_from_list(bam_list_file):
    """
    Return BAM paths from a TSV with sample_id and bam columns.
    """
    if not bam_list_file or not os.path.exists(bam_list_file):
        return []

    try:
        bam_df = pd.read_csv(bam_list_file, sep="\t", comment="#")
    except Exception as exc:
        logger.warning(f"Could not read BAM list for freshness check: {exc}")
        return []

    if "bam" not in bam_df.columns:
        return []

    return [
        str(path).strip()
        for path in bam_df["bam"].tolist()
        if pd.notna(path) and str(path).strip()
    ]


def sanitize_filename_token(value):
    """
    Convert a sample identifier into a conservative filename token.
    """
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


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
        print(line, end='', flush=True)  # Print to terminal in real-time
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
        sys.executable,
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
        sys.executable,
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


def compute_shared_offsets(
    annotated_clustered_file,
    output_dir,
    force_rerun=False,
    site_depth_offsets_file=None,
    offset_mode="cluster_max",
):
    """
    Compute shared offsets from the full clustered matrix.
    These offsets will be used by both donor and acceptor analyses.
    
    Args:
        annotated_clustered_file: Annotated clustered matrix with both donor and acceptor clusters
        output_dir: Output directory
        force_rerun: If True, rerun even if output exists
        site_depth_offsets_file: Optional site-depth offset matrix
        offset_mode: Offset denominator strategy
        
    Returns:
        Path to shared offsets file
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    
    shared_offsets_file = os.path.join(output_dir, "shared_offsets.tsv")
    
    shared_offset_metadata_file = (
        shared_offsets_file[:-4] + ".metadata.tsv"
        if shared_offsets_file.endswith(".tsv")
        else shared_offsets_file + ".metadata.tsv"
    )

    offsets_current = file_is_current(
        shared_offsets_file,
        [annotated_clustered_file, site_depth_offsets_file],
    )
    metadata_mode_matches = False
    if file_exists_and_valid(shared_offset_metadata_file):
        try:
            metadata_df = pd.read_csv(shared_offset_metadata_file, sep="\t", index_col=0, nrows=20)
            metadata_mode_matches = (
                "offset_mode" in metadata_df.columns and
                metadata_df["offset_mode"].dropna().astype(str).eq(offset_mode).all()
            )
        except Exception:
            metadata_mode_matches = False

    if not force_rerun and offsets_current and metadata_mode_matches:
        logger.info("=== Computing shared offsets ===")
        logger.info(f"SKIPPING - Shared offsets file already exists")
        return shared_offsets_file
    
    cmd = [
        sys.executable,
        os.path.join(util_dir, "compute_offsets.py"),
        "--matrix", annotated_clustered_file,
        "--output", shared_offsets_file,
        "--compute_offsets_only",  # New flag to only compute offsets
    ]

    if site_depth_offsets_file:
        cmd.extend(["--site_depth_offsets", site_depth_offsets_file])

    if offset_mode:
        cmd.extend(["--offset_mode", offset_mode])
    
    run_command(cmd, "Computing shared offsets from full clustered matrix")
    
    return shared_offsets_file


def compute_site_depth_offsets(
    matrix_file,
    output_dir,
    bam_list_file,
    window_radius=10,
    min_mapq=60,
    force_rerun=False,
):
    """
    Compute site-depth offsets from BAM files for the introns in matrix_file.

    Args:
        matrix_file: Intron matrix whose rows define the introns to score
        output_dir: Work directory for output
        bam_list_file: TSV with sample_id and bam columns
        window_radius: Bases on each side of splice-site coordinate
        min_mapq: Minimum mapping quality for reads
        force_rerun: If True, recompute even if output exists

    Returns:
        Path to site-depth offset matrix
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    site_depth_offsets_file = os.path.join(output_dir, "site_depth_offsets.tsv")

    if (
        not force_rerun and
        file_is_current(
            site_depth_offsets_file,
            [matrix_file, bam_list_file] + get_bam_paths_from_list(bam_list_file),
        )
    ):
        logger.info("=== Computing site-depth offsets ===")
        logger.info(f"SKIPPING - Site-depth offsets already exist: {site_depth_offsets_file}")
        return site_depth_offsets_file

    cmd = [
        sys.executable,
        os.path.join(util_dir, "compute_splice_site_depth_offsets.py"),
        "--matrix", matrix_file,
        "--bam_list", bam_list_file,
        "--output", site_depth_offsets_file,
        "--window_radius", str(window_radius),
        "--min_mapq", str(min_mapq),
    ]

    run_command(cmd, "Computing splice-site depth offsets from BAMs")
    return site_depth_offsets_file


def prepare_edgeR_inputs(
    filtered_file,
    output_dir,
    shared_offsets_file,
    samples_file=None,
    force_rerun=False,
    offset_metadata_file=None,
):
    """
    Prepare edgeR input files using pre-computed shared offsets.
    
    Args:
        filtered_file: Filtered intron matrix
        output_dir: Output directory
        shared_offsets_file: Path to shared offsets file
        samples_file: Path to sample metadata file (used to filter samples if needed, but samples
                     should already be filtered at clustering stage)
        force_rerun: If True, rerun even if outputs exist
        offset_metadata_file: Optional metadata from shared offset computation
        
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
    
    input_paths = [filtered_file, shared_offsets_file]
    if offset_metadata_file:
        input_paths.append(offset_metadata_file)
    if samples_file:
        input_paths.append(samples_file)
    all_exist = all(file_is_current(f, input_paths) for f in output_files.values())
    
    if not force_rerun and all_exist:
        logger.info(f"=== Preparing edgeR inputs ===")
        logger.info(f"SKIPPING - All output files already exist")
        return output_files
    
    cmd = [
        sys.executable,
        os.path.join(util_dir, "compute_offsets.py"),
        "--matrix", filtered_file,
        "--output_prefix", output_prefix,
        "--shared_offsets", shared_offsets_file,
    ]
    
    # Add sample filtering if metadata file provided (safety check, samples should already be filtered)
    if samples_file:
        cmd.extend(["--samples", samples_file])

    if offset_metadata_file and file_exists_and_valid(offset_metadata_file):
        cmd.extend(["--offset_metadata", offset_metadata_file])
    
    run_command(cmd, f"Preparing edgeR inputs")
    
    return output_files


def prepare_inputs_from_bam_manifest(
    samples_manifest,
    genome_fa,
    workdir,
    force_rerun=False,
    site_depth_window_radius=10,
    min_mapping_quality=60,
):
    """
    Count introns from BAMs and build count/site-depth offset matrices.

    The manifest must contain:
        sample_type    replicate_id    bam_file

    Returns:
        Dict with matrix, site_depth_offsets, samples metadata, and BAM manifest paths.
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")

    logger.info(f"Loading BAM sample manifest from {samples_manifest}")
    manifest_df = pd.read_csv(samples_manifest, sep="\t", comment="#")
    required_cols = {"sample_type", "replicate_id", "bam_file"}
    missing_cols = required_cols - set(manifest_df.columns)
    if missing_cols:
        raise ValueError(
            "BAM manifest must contain columns sample_type, replicate_id, bam_file; "
            f"missing: {', '.join(sorted(missing_cols))}"
        )

    manifest_df = manifest_df.copy()
    manifest_df["sample_type"] = manifest_df["sample_type"].astype(str)
    manifest_df["replicate_id"] = manifest_df["replicate_id"].astype(str)
    manifest_df["bam_file"] = manifest_df["bam_file"].astype(str)

    if manifest_df["replicate_id"].duplicated().any():
        duplicated = manifest_df.loc[manifest_df["replicate_id"].duplicated(), "replicate_id"].tolist()
        raise ValueError(f"replicate_id values must be unique; duplicated: {', '.join(duplicated[:5])}")

    unsafe_ids = [
        sample_id for sample_id in manifest_df["replicate_id"]
        if sanitize_filename_token(sample_id) != sample_id or "." in sample_id
    ]
    if unsafe_ids:
        raise ValueError(
            "replicate_id values must be usable as matrix sample names and filename prefixes "
            "(letters, numbers, underscore, or dash; no dots/spaces). Invalid examples: "
            + ", ".join(unsafe_ids[:5])
        )

    missing_bams = [path for path in manifest_df["bam_file"] if not os.path.exists(path)]
    if missing_bams:
        raise FileNotFoundError(
            "BAM files not found: " + ", ".join(missing_bams[:5]) +
            ("..." if len(missing_bams) > 5 else "")
        )

    if not os.path.exists(genome_fa):
        raise FileNotFoundError(f"Genome FASTA not found: {genome_fa}")

    input_dir = os.path.join(workdir, "bam_inputs")
    discovery_dir = os.path.join(input_dir, "discovery_introns")
    targeted_dir = os.path.join(input_dir, "targeted_introns")
    os.makedirs(discovery_dir, exist_ok=True)
    os.makedirs(targeted_dir, exist_ok=True)

    discovery_files = []
    targeted_files = []

    logger.info("Counting observed introns from BAMs (discovery pass)")
    for _, row in manifest_df.iterrows():
        sample_id = row["replicate_id"]
        sample_token = sanitize_filename_token(sample_id)
        bam_file = row["bam_file"]
        output_file = os.path.join(discovery_dir, f"{sample_token}.introns")
        discovery_files.append(output_file)

        cmd = (
            f"{shlex.quote(sys.executable)} "
            f"{shlex.quote(os.path.join(util_dir, 'count_introns_from_bam.py'))} "
            f"--genome_fa {shlex.quote(genome_fa)} "
            f"--bam {shlex.quote(bam_file)} "
            f"--site_depth_window_radius {int(site_depth_window_radius)} "
            f"--min_mapping_quality {int(min_mapping_quality)} "
            f"> {shlex.quote(output_file)}"
        )
        run_command(
            cmd,
            f"Counting introns for {sample_id} (discovery)",
            skip_if_exists=None if force_rerun else output_file,
        )

    discovery_matrix = os.path.join(input_dir, "discovery_intron_counts.matrix")
    discovery_offset_matrix = os.path.join(input_dir, "discovery_intron_counts.offsets.matrix")
    cmd = [
        sys.executable,
        os.path.join(util_dir, "build_intron_count_matrix.py"),
        "--intron_files",
        *discovery_files,
        "--output_matrix",
        discovery_matrix,
        "--output_offset_matrix",
        discovery_offset_matrix,
    ]
    run_command(
        cmd,
        "Building discovery intron matrix",
        skip_if_exists=None if force_rerun else discovery_matrix,
    )

    logger.info("Counting target introns from BAMs (complete site-depth offset pass)")
    for _, row in manifest_df.iterrows():
        sample_id = row["replicate_id"]
        sample_token = sanitize_filename_token(sample_id)
        bam_file = row["bam_file"]
        output_file = os.path.join(targeted_dir, f"{sample_token}.introns")
        targeted_files.append(output_file)

        cmd = (
            f"{shlex.quote(sys.executable)} "
            f"{shlex.quote(os.path.join(util_dir, 'count_introns_from_bam.py'))} "
            f"--genome_fa {shlex.quote(genome_fa)} "
            f"--bam {shlex.quote(bam_file)} "
            f"--target_introns {shlex.quote(discovery_matrix)} "
            f"--site_depth_window_radius {int(site_depth_window_radius)} "
            f"--min_mapping_quality {int(min_mapping_quality)} "
            f"> {shlex.quote(output_file)}"
        )
        run_command(
            cmd,
            f"Counting introns for {sample_id} (targeted)",
            skip_if_exists=None if force_rerun else output_file,
        )

    final_matrix = os.path.join(input_dir, "intron_counts.matrix")
    final_offset_matrix = os.path.join(input_dir, "intron_counts.offsets.matrix")
    cmd = [
        sys.executable,
        os.path.join(util_dir, "build_intron_count_matrix.py"),
        "--intron_files",
        *targeted_files,
        "--output_matrix",
        final_matrix,
        "--output_offset_matrix",
        final_offset_matrix,
    ]
    run_command(
        cmd,
        "Building final intron count and site-depth offset matrices",
        skip_if_exists=None if force_rerun else (
            final_offset_matrix if file_exists_and_valid(final_matrix) else None
        ),
    )

    downstream_samples = os.path.join(input_dir, "sample_metadata.tsv")
    if force_rerun or not file_exists_and_valid(downstream_samples):
        logger.info(f"Writing downstream sample metadata to {downstream_samples}")
        output_df = manifest_df.copy()
        output_df.insert(0, "sample_id", output_df["replicate_id"])
        output_df.insert(1, "group", output_df["sample_type"])
        output_df.to_csv(downstream_samples, sep="\t", index=False)

    bam_list = os.path.join(input_dir, "bam_list.tsv")
    if force_rerun or not file_exists_and_valid(bam_list):
        logger.info(f"Writing BAM list to {bam_list}")
        manifest_df[["replicate_id", "bam_file"]].rename(
            columns={"replicate_id": "sample_id", "bam_file": "bam"}
        ).to_csv(bam_list, sep="\t", index=False)

    return {
        "matrix": final_matrix,
        "site_depth_offsets": final_offset_matrix,
        "samples": downstream_samples,
        "bam_list": bam_list,
    }


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


def infer_single_psi_contrast(samples_file, group_col, edgeR_params):
    """
    Infer a single contrast string for PSI directionality when one is not passed
    explicitly to edgeR.

    This is only safe when the analysis design implies exactly one treatment vs
    one control/composite control comparison. For multi-contrast analyses we
    return None rather than guessing.
    """
    explicit_contrast = edgeR_params.get("contrast")
    if explicit_contrast:
        return explicit_contrast

    groups = get_groups_from_samples(samples_file, group_col)

    control_groups_raw = edgeR_params.get("control_groups")
    if control_groups_raw:
        control_groups = [g.strip() for g in control_groups_raw.split(",") if g.strip()]
        non_control_groups = [g for g in groups if g not in control_groups]
        if len(non_control_groups) == 1 and len(control_groups) >= 1:
            if len(control_groups) == 1:
                return f"{non_control_groups[0]},{control_groups[0]}"
            return f"{non_control_groups[0]},{';'.join(control_groups)}"

    if len(groups) == 2:
        return f"{groups[0]},{groups[1]}"

    return None


def apply_contrast_ordered_delta_psi(results_df, psi_df):
    """
    Recompute delta_PSI using each row's contrast direction.

    This keeps delta_PSI aligned with the logFC direction even when group labels
    do not sort the way the contrast was specified, and it also supports
    combined multi-contrast result tables.
    """
    if "contrast" not in results_df.columns:
        return results_df

    def compute_row_delta(row):
        contrast = row.get("contrast")
        if pd.isna(contrast) or "_vs_" not in contrast:
            return np.nan

        group1, group2_raw = contrast.split("_vs_", 1)
        group2_list = [g for g in group2_raw.split("_") if g]

        needed_cols = [f"{group1}_mean_PSI"] + [f"{g}_mean_PSI" for g in group2_list]
        if any(col not in row.index for col in needed_cols):
            return np.nan

        group1_value = row[f"{group1}_mean_PSI"]
        group2_values = [row[f"{g}_mean_PSI"] for g in group2_list]

        if any(pd.isna([group1_value, *group2_values])):
            return np.nan

        return group1_value - float(np.mean(group2_values))

    results_df = results_df.copy()
    results_df["delta_PSI"] = results_df.apply(compute_row_delta, axis=1)
    return results_df


def run_single_contrast(contrast, edgeR_inputs, samples_file, output_dir, edgeR_params, force_rerun):
    """
    Run edgeR analysis for a single contrast.
    
    Args:
        contrast: Contrast string (e.g., "GroupA,GroupB" where log2FC = GroupA/GroupB)
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
    contrast_safe = contrast.replace(",", "_vs_")
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
    sig_results_file = f"{output_prefix}.significant_introns.tsv"
    
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
                contrasts = [f"{g},{control_groups[0]}" for g in non_control_groups]
            else:
                control_str = ";".join(control_groups)  # Use semicolon for multiple controls
                contrasts = [f"{g},{control_str}" for g in non_control_groups]
            
            logger.info("=== Running edgeR analysis ===")
            logger.info(f"Control-based comparisons: {len(contrasts)} contrasts")
            logger.info(f"Control groups: {', '.join(control_groups)}")
            logger.info(f"Treatment groups: {', '.join(non_control_groups)}")
            logger.info(f"Using {cpu} CPU(s)")
        else:
            # Original behavior: all pairwise comparisons
            contrasts = [f"{g1},{g2}" for g1, g2 in combinations(groups, 2)]
            
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

            if 'significant' in combined.columns:
                combined_sig = combined[combined['significant']].copy()
                if len(combined_sig) > 0:
                    combined_sig.to_csv(sig_results_file, sep="\t", index=False, na_rep='NA')
                    logger.info(f"Combined significant introns written to: {sig_results_file}")
                elif os.path.exists(sig_results_file):
                    os.remove(sig_results_file)
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
        
        # Get contrast from parameters or infer a single unambiguous contrast.
        contrast = infer_single_psi_contrast(samples_file, edgeR_params["group_col"], edgeR_params)
        
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


def add_psi_and_filter(
    intron_results_file,
    psi_file,
    output_dir,
    min_delta_psi=None,
    force_rerun=False,
    fdr_threshold=0.05,
    min_logFC=0.0,
):
    """
    Add PSI values to edgeR results and optionally filter by delta PSI with FDR recalculation.
    
    Writes intermediate PSI/edgeR artifacts into output_dir and promotes only the
    primary user-facing result files to the parent directory.
    
    Args:
        intron_results_file: Path to intron results from edgeR
        psi_file: Path to PSI values file
        output_dir: Working directory for intermediate PSI-enhanced outputs
        min_delta_psi: Minimum absolute delta PSI to include (with FDR recalculation)
        force_rerun: If True, rerun even if outputs exist
        fdr_threshold: FDR threshold used to define significance
        min_logFC: Minimum absolute log2FC used to define significance
        
    Returns:
        Tuple of (all_results_file, significant_results_file)
    """
    if not psi_file or not file_exists_and_valid(psi_file):
        logger.warning("No PSI file available, skipping PSI annotation")
        return (intron_results_file, None)
    
    final_output_dir = os.path.dirname(output_dir.rstrip(os.sep))
    output_prefix = os.path.join(output_dir, "edgeR_results")
    final_output_prefix = os.path.join(final_output_dir, "edgeR_results")
    
    # Intermediate files in workdir
    unfiltered_file = f"{output_prefix}.intron_results_with_psi.tsv"
    filtered_file = f"{output_prefix}.intron_results_with_psi.psi_filtered.tsv"

    # User-facing files in main output dir
    final_all_file = f"{final_output_prefix}.all.tsv"
    final_sig_file = f"{final_output_prefix}.significant_introns.tsv"
    
    logger.info("=== Adding PSI to results ===")
    
    try:
        # Load results and PSI
        results_df = pd.read_csv(intron_results_file, sep="\t")
        psi_df = pd.read_csv(psi_file, sep="\t", index_col=0)
        
        # Keep only summary PSI columns (not per-sample)
        psi_summary_cols = [col for col in psi_df.columns
                           if 'mean_PSI' in col or 'median_PSI' in col or
                              'std_PSI' in col]
        
        # Merge PSI data with results
        results_with_index = results_df.set_index('intron_id')
        
        for col in psi_summary_cols:
            if col in psi_df.columns:
                results_with_index[col] = psi_df.loc[results_with_index.index, col]
        
        results_with_psi = results_with_index.reset_index()
        results_with_psi = apply_contrast_ordered_delta_psi(results_with_psi, psi_df)
        
        # Write intermediate and final unfiltered PSI-enhanced results
        logger.info(f"Writing PSI-enhanced results to {unfiltered_file}")
        results_with_psi.to_csv(unfiltered_file, sep="\t", index=False, na_rep='NA')
        logger.info(f"Writing primary results to {final_all_file}")
        results_with_psi.to_csv(final_all_file, sep="\t", index=False, na_rep='NA')
        logger.info(f"Added {len(psi_summary_cols)} PSI columns to results")

        sig_source = results_with_psi.copy()
        
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
                pd.DataFrame(columns=results_with_psi.columns).to_csv(filtered_file, sep="\t", index=False)
                pd.DataFrame(columns=results_with_psi.columns).to_csv(final_sig_file, sep="\t", index=False)
                return (final_all_file, final_sig_file)
            
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
                    
                    # Recalculate significance based on the PSI-filtered FDR
                    # while retaining the edgeR logFC effect-size threshold.
                    filtered_results['significant'] = (
                        (filtered_results['FDR'] <= fdr_threshold) &
                        (filtered_results['logFC'].abs() >= min_logFC)
                    )
                    
                    sig_before = (
                        (filtered_results['FDR_original'] <= fdr_threshold) &
                        (filtered_results['logFC'].abs() >= min_logFC)
                    ).sum()
                    sig_after = filtered_results['significant'].sum()
                    
                    logger.info(
                        "Significant introns before FDR recalculation "
                        f"(FDR <= {fdr_threshold}, |logFC| >= {min_logFC}): {sig_before}"
                    )
                    logger.info(
                        "Significant introns after FDR recalculation "
                        f"(FDR <= {fdr_threshold}, |logFC| >= {min_logFC}): {sig_after}"
                    )
                    logger.info(f"Gained {sig_after - sig_before} significant introns from reduced multiple testing burden")
            
            logger.info(f"Writing intermediate PSI-filtered results to {filtered_file}")
            filtered_results.to_csv(filtered_file, sep="\t", index=False, na_rep='NA')
            sig_source = filtered_results

        final_sig_results = sig_source[sig_source['significant']].copy() if 'significant' in sig_source.columns else pd.DataFrame(columns=sig_source.columns)

        total_final_tests = len(sig_source)
        pct_significant = 100 * len(final_sig_results) / total_final_tests if total_final_tests else 0
        criteria = [f"FDR <= {fdr_threshold}", f"|logFC| >= {min_logFC}"]
        if min_delta_psi and 'delta_PSI' in sig_source.columns:
            criteria.append(f"|delta_PSI| >= {min_delta_psi}")

        logger.info("=== Results Summary ===")
        logger.info(f"Final introns evaluated: {total_final_tests}")
        logger.info(
            "Significant introns "
            f"({', '.join(criteria)}): {len(final_sig_results)} ({pct_significant:.1f}%)"
        )

        if 'logFC' in final_sig_results.columns:
            sig_up = (final_sig_results['logFC'] > 0).sum()
            sig_down = (final_sig_results['logFC'] < 0).sum()
            direction_label = "usage"
            if 'contrast' in sig_source.columns:
                contrasts = sig_source['contrast'].dropna().unique()
                if len(contrasts) == 1 and "_vs_" in contrasts[0]:
                    group1, group2 = contrasts[0].split("_vs_", 1)
                    direction_label = f"usage ({group1} vs {group2})"
            logger.info(f"  Increased {direction_label}: {sig_up}")
            logger.info(f"  Decreased {direction_label}: {sig_down}")

        logger.info(f"Writing significant results to {final_sig_file}")
        final_sig_results.to_csv(final_sig_file, sep="\t", index=False, na_rep='NA')
        
        return (final_all_file, final_sig_file)
        
    except Exception as e:
        logger.error(f"Error adding PSI to results: {e}")
        logger.warning("Failed to add PSI, will use original results")
        return (intron_results_file, None)


def cleanup_legacy_top_level_outputs(output_dir):
    """
    Remove legacy top-level artifacts so only the current primary result files remain.
    """
    legacy_paths = [
        os.path.join(output_dir, "edgeR_results.intron_results.tsv"),
        os.path.join(output_dir, "edgeR_results.intron_results_with_psi.tsv"),
        os.path.join(output_dir, "edgeR_results.intron_results_with_psi.psi_filtered.tsv"),
        os.path.join(output_dir, "edgeR_results.diagnostics.pdf"),
        os.path.join(output_dir, "edgeR_results.RData"),
        os.path.join(output_dir, "psi.psi_values.tsv"),
    ]

    legacy_paths.extend(glob.glob(os.path.join(output_dir, "*_edgeR_results.intron_results.tsv")))
    legacy_paths.extend(glob.glob(os.path.join(output_dir, "*_edgeR_results.significant_introns.tsv")))
    legacy_paths.extend(glob.glob(os.path.join(output_dir, "*_edgeR_results.diagnostics.pdf")))
    legacy_paths.extend(glob.glob(os.path.join(output_dir, "*_edgeR_results.RData")))

    for path in sorted(set(legacy_paths)):
        if os.path.exists(path):
            os.remove(path)
            logger.info(f"Removed legacy top-level output: {path}")




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
        default=None,
        help=(
            "Input intron count matrix (supports .tsv or .tsv.gz). "
            "If omitted, --samples must be a BAM manifest and --genome_fa is required."
        ),
    )
    
    parser.add_argument(
        "--samples",
        type=str,
        required=True,
        help=(
            "Sample metadata TSV. Matrix mode expects sample_id and group columns. "
            "BAM mode expects sample_type, replicate_id, and bam_file columns."
        ),
    )
    
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory for results",
    )

    parser.add_argument(
        "--genome_fa",
        type=str,
        default=None,
        help="Reference genome FASTA required when --matrix is omitted and introns are counted from BAMs",
    )
    
    parser.add_argument(
        "--gtf",
        type=str,
        default=None,
        help="GTF annotation file for gene annotation and known/novel intron status (optional)",
    )

    parser.add_argument(
        "--site_depth_offsets",
        type=str,
        default=None,
        help=(
            "Optional precomputed intron x sample raw site-depth offset matrix. "
            "In matrix mode, provide this file directly, or provide --site_depth_bam_list "
            "to compute it during the run."
        ),
    )

    parser.add_argument(
        "--site_depth_bam_list",
        type=str,
        default=None,
        help=(
            "Matrix-mode TSV with sample_id and bam columns. When provided, the pipeline "
            "computes workdir/site_depth_offsets.tsv and uses it for site-depth offset modes. "
            "When --matrix is omitted, put BAM paths in --samples instead."
        ),
    )

    parser.add_argument(
        "--site_depth_window_radius",
        type=int,
        default=10,
        help="Bases on each side of each splice-site coordinate to include for site-depth offsets",
    )

    parser.add_argument(
        "--site_depth_min_mapq",
        type=int,
        default=60,
        help="Minimum mapping quality for reads counted in site-depth offsets",
    )

    parser.add_argument(
        "--offset_mode",
        choices=[
            "auto",
            "cluster_max",
            "site_depth",
            "cluster_with_site_singleton_fallback",
        ],
        default="auto",
        help=(
            "Offset denominator strategy. auto uses site_depth in BAM-manifest mode, "
            "or when --site_depth_offsets/--site_depth_bam_list is provided; "
            "otherwise auto uses cluster_max."
        ),
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
        help="Contrast to test (e.g., 'perturb,control' where log2FC = perturb/control)",
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
        default=0.05,
        help="Minimum absolute delta PSI to include in final results. Defaults to 0.05. Set to 0 to disable PSI filtering. FDR will be recalculated on the filtered set (reduces multiple testing burden). Example: 0.1 for 10%% change",
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
    if args.site_depth_offsets and args.site_depth_bam_list:
        parser.error("Use either --site_depth_offsets or --site_depth_bam_list, not both.")
    bam_input_mode = args.matrix is None
    if not args.matrix and not args.genome_fa:
        parser.error("--genome_fa is required when --matrix is omitted and BAMs are counted from --samples")
    if bam_input_mode and (args.site_depth_offsets or args.site_depth_bam_list):
        parser.error(
            "When --matrix is omitted, site-depth offsets are generated from the BAMs listed in --samples; "
            "do not also provide --site_depth_offsets or --site_depth_bam_list."
        )

    offset_mode = args.offset_mode
    if offset_mode == "auto":
        offset_mode = "site_depth" if (args.site_depth_offsets or args.site_depth_bam_list or bam_input_mode) else "cluster_max"
    if (
        offset_mode in {"site_depth", "cluster_with_site_singleton_fallback"} and
        not (args.site_depth_offsets or args.site_depth_bam_list or bam_input_mode)
    ):
        parser.error(f"--offset_mode {offset_mode} requires --site_depth_offsets or --site_depth_bam_list")
    
    # Create output directory and workdir for intermediates
    os.makedirs(args.output_dir, exist_ok=True)
    workdir = os.path.join(args.output_dir, "workdir")
    os.makedirs(workdir, exist_ok=True)

    matrix_file = args.matrix
    samples_file = args.samples
    site_depth_offsets_file = args.site_depth_offsets

    if bam_input_mode:
        logger.info("=== Preparing count and offset matrices from BAM manifest ===")
        prepared_inputs = prepare_inputs_from_bam_manifest(
            samples_manifest=args.samples,
            genome_fa=args.genome_fa,
            workdir=workdir,
            force_rerun=args.force_rerun,
            site_depth_window_radius=args.site_depth_window_radius,
            min_mapping_quality=args.site_depth_min_mapq,
        )
        matrix_file = prepared_inputs["matrix"]
        samples_file = prepared_inputs["samples"]
        site_depth_offsets_file = prepared_inputs["site_depth_offsets"]
    
    logger.info("=== Differential Splicing Analysis Pipeline ===")
    logger.info(f"Input mode: {'BAM manifest' if bam_input_mode else 'matrix'}")
    logger.info(f"Input matrix: {matrix_file}")
    logger.info(f"Sample metadata: {samples_file}")
    if site_depth_offsets_file:
        logger.info(f"Site-depth offsets: {site_depth_offsets_file}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Work directory (intermediates): {workdir}")
    logger.info("Analysis mode: Intron-level with shared offsets")
    logger.info(f"Offset mode: {offset_mode}")
    
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
    clustered_file = os.path.join(workdir, "introns_clustered.tsv")
    
    cmd = [
        sys.executable,
        os.path.join(util_dir, "cluster_introns.py"),
        "--matrix", matrix_file,
        "--output_donor", clustered_file,  # Will write after both donor and acceptor clustering
        "--cluster_type", "both",
        "--samples", samples_file,  # Filter samples at clustering step
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

    # Optional Step 2b: Compute site-depth offsets from BAMs if requested in matrix mode.
    if args.site_depth_bam_list:
        logger.info(f"\n{'='*60}")
        logger.info("Computing site-depth offsets")
        logger.info(f"{'='*60}\n")

        site_depth_offsets_file = compute_site_depth_offsets(
            annotated_clustered,
            workdir,
            args.site_depth_bam_list,
            window_radius=args.site_depth_window_radius,
            min_mapq=args.site_depth_min_mapq,
            force_rerun=args.force_rerun,
        )
        logger.info(f"Site-depth offsets file: {site_depth_offsets_file}")
    
    # Step 3: Compute shared offsets from full clustered matrix
    logger.info(f"\n{'='*60}")
    logger.info("Computing shared offsets")
    logger.info(f"{'='*60}\n")
    
    shared_offsets_file = compute_shared_offsets(
        annotated_clustered,
        workdir,
        force_rerun=args.force_rerun,
        site_depth_offsets_file=site_depth_offsets_file,
        offset_mode=offset_mode,
    )
    logger.info(f"Shared offsets file: {shared_offsets_file}")
    shared_offset_metadata_file = (
        shared_offsets_file[:-4] + ".metadata.tsv"
        if shared_offsets_file.endswith(".tsv")
        else shared_offsets_file + ".metadata.tsv"
    )
    
    # Step 4: Filter introns (require thresholds for both donor and acceptor clusters)
    logger.info(f"\n{'='*60}")
    logger.info("Filtering introns")
    logger.info(f"{'='*60}\n")
    
    filtered_file = os.path.join(workdir, "introns_filtered.tsv")
    
    cmd = [
        sys.executable,
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
        filtered_file, workdir,
        shared_offsets_file,
        samples_file=samples_file,
        force_rerun=args.force_rerun,
        offset_metadata_file=shared_offset_metadata_file,
    )
    
    # Step 6: Run edgeR analysis
    logger.info(f"\n{'='*60}")
    logger.info("Running edgeR analysis")
    logger.info(f"{'='*60}\n")
    
    intron_results = run_edgeR(
        edgeR_inputs, samples_file, workdir, edgeR_params,  # Write raw edgeR outputs to workdir
        force_rerun=args.force_rerun,
        cpu=args.cpu
    )
    
    # Step 7: Compute PSI values using shared offsets
    logger.info(f"\n{'='*60}")
    logger.info("Computing PSI values")
    logger.info(f"{'='*60}\n")
    
    psi_file = compute_psi_for_results(
        edgeR_inputs, samples_file, workdir,  # Write PSI to workdir
        edgeR_params, 
        shared_offsets_file=shared_offsets_file,
        force_rerun=args.force_rerun
    )
    
    # Step 8: Add PSI to results and optionally filter by delta PSI
    logger.info(f"\n{'='*60}")
    logger.info("Adding PSI to results")
    logger.info(f"{'='*60}\n")
    
    intron_results_with_psi, significant_results = add_psi_and_filter(
        intron_results, psi_file, workdir,  # Keep intermediates in workdir and promote final outputs
        min_delta_psi=args.min_delta_psi,
        force_rerun=args.force_rerun,
        fdr_threshold=args.fdr_threshold,
        min_logFC=args.min_logFC,
    )

    cleanup_legacy_top_level_outputs(args.output_dir)
    
    logger.info("\n" + "="*60)
    logger.info("PIPELINE COMPLETE!")
    logger.info("="*60)
    logger.info(f"\nFinal results saved to: {args.output_dir}")
    logger.info(f"Intermediate files saved to: {workdir}")
    
    # Print summary of key output files
    logger.info("\nFinal output files:")
    if file_exists_and_valid(f"{args.output_dir}/edgeR_results.all.tsv"):
        logger.info(f"  - All intron results with PSI: {args.output_dir}/edgeR_results.all.tsv")
    sig_file = f"{args.output_dir}/edgeR_results.significant_introns.tsv"
    if file_exists_and_valid(sig_file):
        logger.info(f"  - PSI-filtered significant introns: {sig_file}")
    
    logger.info(f"\nIntermediate files (in workdir/):")
    logger.info(f"  - Clustered introns: {workdir}/introns_clustered.tsv")
    logger.info(f"  - Filtered introns: {workdir}/introns_filtered.tsv")
    logger.info(f"  - Shared offsets: {workdir}/shared_offsets.tsv")
    logger.info(f"  - Raw edgeR results: {workdir}/edgeR_results.intron_results.tsv")
    logger.info(f"  - PSI values: {workdir}/psi.psi_values.tsv")
    if file_exists_and_valid(f"{workdir}/edgeR_results.intron_results_with_psi.psi_filtered.tsv"):
        logger.info(f"  - PSI-filtered full results: {workdir}/edgeR_results.intron_results_with_psi.psi_filtered.tsv")
    logger.info(f"  - Diagnostics plots: {workdir}/edgeR_results.diagnostics.pdf")


if __name__ == "__main__":
    main()
