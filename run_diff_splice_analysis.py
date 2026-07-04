#!/usr/bin/env python3

"""
Main pipeline orchestrator for differential splicing analysis.

Coordinates the full workflow:
1. Load intron count matrix
2. Load or compute splice-site depth offsets from the supported input mode
3. Filter introns by count, selected denominator depth, and pre-test delta PSI thresholds
4. Run the selected statistical model
5. Add PSI summaries to model results
"""

import sys
import os
import argparse
import logging
import subprocess
import glob
import shlex
import re
import json
import pandas as pd
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

MATRIX_METADATA_COLS = {
    "intron_info",
    "donor_cluster",
    "acceptor_cluster",
    "gene_name",
    "intron_status",
    "overlapping_genes",
    "donor_cluster_size",
    "acceptor_cluster_size",
    "both_splice_sites_singleton",
    "offset_mode",
    "offset_source",
    "site_depth_fallback_used",
    "delta_PSI",
}

OFFSET_MODE_SPECS = {
    "splice_plus_retained": {
        "denominator_key": "max_splice_plus_retained_depth",
        "psi_denominator_label": "splice_plus_retained",
        "stat_engine": "edgeR",
    },
}


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


def parse_contrast(contrast):
    """
    Parse the single supported contrast format: GroupA,GroupB.
    """
    if not contrast:
        raise ValueError("--contrast is required")

    parts = [part.strip() for part in contrast.split(",", 1)]
    if len(parts) != 2 or not all(parts):
        raise ValueError("--contrast must be in format 'GroupA,GroupB'")

    group1, group2 = parts
    if ";" in group1 or ";" in group2:
        raise ValueError("Only individual contrasts are supported; do not use semicolon-separated groups")

    return group1, group2


def get_sample_columns(df):
    """
    Identify sample columns in a count matrix.
    """
    return [col for col in df.columns if col not in MATRIX_METADATA_COLS]


def load_and_align_counts_offsets(matrix_file, site_depth_offsets_file, samples_file, group_col):
    """
    Load count/offset matrices and subset both to samples listed in metadata.
    """
    logger.info(f"Loading intron count matrix from {matrix_file}")
    counts_df = pd.read_csv(matrix_file, sep="\t", index_col=0)
    sample_cols = get_sample_columns(counts_df)
    if not sample_cols:
        raise ValueError(f"No sample columns found in count matrix: {matrix_file}")

    logger.info(f"Loading denominator matrix from {site_depth_offsets_file}")
    offsets_df = pd.read_csv(site_depth_offsets_file, sep="\t", index_col=0)

    logger.info(f"Loading sample metadata from {samples_file}")
    samples_df = pd.read_csv(samples_file, sep="\t", comment="#")
    required_cols = {"sample_id", group_col}
    missing_cols = required_cols - set(samples_df.columns)
    if missing_cols:
        raise ValueError(
            f"Sample metadata must contain columns: {', '.join(sorted(required_cols))}; "
            f"missing: {', '.join(sorted(missing_cols))}"
        )

    samples_df = samples_df.copy()
    samples_df["sample_id"] = samples_df["sample_id"].astype(str)
    samples_df[group_col] = samples_df[group_col].astype(str)

    metadata_samples = [
        sample_id for sample_id in samples_df["sample_id"].tolist()
        if sample_id and sample_id.lower() != "nan"
    ]
    missing_count_samples = [sample for sample in metadata_samples if sample not in sample_cols]
    missing_offset_samples = [sample for sample in metadata_samples if sample not in offsets_df.columns]
    if missing_count_samples:
        raise ValueError(
            "Samples from metadata missing in count matrix: "
            + ", ".join(missing_count_samples[:10])
        )
    if missing_offset_samples:
        raise ValueError(
            "Samples from metadata missing in denominator matrix: "
            + ", ".join(missing_offset_samples[:10])
        )

    missing_introns = counts_df.index.difference(offsets_df.index)
    if len(missing_introns) > 0:
        raise ValueError(
            "Site-depth offsets are missing introns from count matrix: "
            + ", ".join(missing_introns[:5])
            + ("..." if len(missing_introns) > 5 else "")
        )

    counts_df = counts_df.loc[:, metadata_samples].apply(pd.to_numeric, errors="coerce").fillna(0)
    offsets_df = offsets_df.loc[counts_df.index, metadata_samples].apply(pd.to_numeric, errors="coerce").fillna(0)

    logger.info(f"Using {len(counts_df)} introns and {len(metadata_samples)} samples")
    return counts_df, offsets_df, samples_df, metadata_samples


def build_basic_intron_annotations(intron_ids):
    """
    Build annotation columns that do not require clustering.
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    sys.path.insert(0, util_dir)
    from cluster_introns import parse_intron_id

    records = []
    for intron_id in intron_ids:
        info = parse_intron_id(intron_id)
        records.append(
            {
                "intron_id": intron_id,
                "chr": info["chr"],
                "start": info["start"],
                "end": info["end"],
                "strand": info["strand"],
                "donor": info["donor"],
                "acceptor": info["acceptor"],
                "splice_pair": info["splice_pair"],
                "splice_flag": info["flag"],
            }
        )

    return pd.DataFrame.from_records(records, index=intron_ids)


def add_gene_annotations(annotation_df, gtf_file):
    """
    Add gene_name/intron_status/overlapping_genes to direct intron annotations.
    """
    if not gtf_file:
        annotation_df["gene_name"] = "."
        annotation_df["intron_status"] = "unknown"
        annotation_df["overlapping_genes"] = "."
        return annotation_df

    util_dir = os.path.join(os.path.dirname(__file__), "util")
    sys.path.insert(0, util_dir)
    from integrate_results import parse_gtf_file, parse_intron_id, build_gene_index, find_overlapping_genes

    logger.info(f"Annotating introns with genes from {gtf_file}")
    gene_map, annotated_introns, transcript_genes, transcript_exons, intron_gene_map = parse_gtf_file(gtf_file)
    gene_index = build_gene_index(gene_map, transcript_genes, transcript_exons)

    gene_names = []
    intron_statuses = []
    overlapping_values = []

    for intron_id in annotation_df.index:
        coords = parse_intron_id(intron_id)
        if coords is None:
            gene_names.append(".")
            intron_statuses.append("unknown")
            overlapping_values.append(".")
            continue

        intron_statuses.append("known" if coords in annotated_introns else "novel")
        overlapping = find_overlapping_genes(coords, gene_index, intron_gene_map)
        if overlapping:
            gene_names.append(overlapping[0])
            overlapping_values.append(",".join(overlapping))
        else:
            gene_names.append(".")
            overlapping_values.append(".")

    annotation_df["gene_name"] = gene_names
    annotation_df["intron_status"] = intron_statuses
    annotation_df["overlapping_genes"] = overlapping_values
    logger.info(
        f"Gene annotation complete: known={sum(x == 'known' for x in intron_statuses)}, "
        f"novel={sum(x == 'novel' for x in intron_statuses)}, "
        f"with_genes={sum(x != '.' for x in gene_names)}"
    )
    return annotation_df


def compute_site_depth_psi_values(counts_df, offsets_df, samples_df, group_col, contrast):
    """
    Compute per-sample PSI and group-level delta PSI from selected denominators.
    """
    group1, group2 = parse_contrast(contrast)
    groups = set(samples_df[group_col].astype(str))
    missing_groups = [group for group in (group1, group2) if group not in groups]
    if missing_groups:
        raise ValueError(
            "Contrast groups not found in sample metadata: "
            + ", ".join(missing_groups)
            + f". Available groups: {', '.join(sorted(groups))}"
        )

    sample_cols = list(counts_df.columns)
    psi_df = pd.DataFrame(index=counts_df.index)
    for sample in sample_cols:
        psi_df[f"{sample}_PSI"] = np.where(
            offsets_df[sample] > 0,
            counts_df[sample] / offsets_df[sample],
            0.0,
        )

    for group in samples_df[group_col].dropna().astype(str).unique():
        group_samples = samples_df.loc[samples_df[group_col].astype(str) == group, "sample_id"].astype(str).tolist()
        group_psi_cols = [f"{sample}_PSI" for sample in group_samples if f"{sample}_PSI" in psi_df.columns]
        if group_psi_cols:
            psi_df[f"{group}_mean_PSI"] = psi_df[group_psi_cols].mean(axis=1)
            psi_df[f"{group}_median_PSI"] = psi_df[group_psi_cols].median(axis=1)
            psi_df[f"{group}_std_PSI"] = psi_df[group_psi_cols].std(axis=1)

    psi_df["delta_PSI"] = psi_df[f"{group1}_mean_PSI"] - psi_df[f"{group2}_mean_PSI"]
    logger.info(f"Computed delta_PSI for contrast: {group1} - {group2}")
    return psi_df


def prepare_site_depth_edgeR_inputs(
    matrix_file,
    site_depth_offsets_file,
    samples_file,
    output_dir,
    filter_params,
    edgeR_params,
    gtf_file=None,
    force_rerun=False,
    offset_mode_label="site_depth",
    psi_denominator_label="site_depth",
    count_unit_label="read",
):
    """
    Filter introns and prepare model input files using selected denominator offsets.
    """
    output_prefix = os.path.join(output_dir, "edgeR_input")
    output_files = {
        "counts": f"{output_prefix}.counts.tsv",
        "offsets": f"{output_prefix}.offsets.tsv",
        "annotations": f"{output_prefix}.annotations.tsv",
        "raw_offsets": os.path.join(output_dir, "site_depth_offsets.filtered.tsv"),
        "filtered_matrix": os.path.join(output_dir, "introns_filtered.tsv"),
        "psi": os.path.join(output_dir, "psi.psi_values.tsv"),
        "params": f"{output_prefix}.filter_params.json",
    }

    input_paths = [matrix_file, site_depth_offsets_file, samples_file, gtf_file]
    params_record = {
        "filter_params": filter_params,
        "edgeR_params": {
            "group_col": edgeR_params.get("group_col"),
            "contrast": edgeR_params.get("contrast"),
        },
        "gtf_file": gtf_file,
        "offset_mode_label": offset_mode_label,
        "psi_denominator_label": psi_denominator_label,
        "count_unit_label": count_unit_label,
    }
    params_match = False
    if file_exists_and_valid(output_files["params"]):
        try:
            with open(output_files["params"], "rt") as params_fh:
                params_match = json.load(params_fh) == params_record
        except Exception:
            params_match = False

    all_exist = all(file_is_current(path, input_paths) for key, path in output_files.items() if key != "params")
    if not force_rerun and all_exist and params_match:
        logger.info("=== Preparing model inputs ===")
        logger.info("SKIPPING - All model input files already exist")
        return output_files

    counts_df, offsets_df, samples_df, sample_cols = load_and_align_counts_offsets(
        matrix_file,
        site_depth_offsets_file,
        samples_file,
        edgeR_params["group_col"],
    )

    annotations_df = build_basic_intron_annotations(counts_df.index)
    annotations_df = add_gene_annotations(annotations_df, gtf_file)
    annotations_df["offset_mode"] = offset_mode_label
    annotations_df["psi_denominator_mode"] = psi_denominator_label
    annotations_df["count_unit"] = count_unit_label

    n_start = len(counts_df)
    keep_mask = pd.Series(True, index=counts_df.index)

    if not filter_params.get("keep_noncanonical", False):
        canonical_mask = annotations_df["splice_flag"].eq("OK")
        logger.info(
            f"Canonical splice filter: kept {int(canonical_mask.sum())}/{n_start} introns"
        )
        keep_mask &= canonical_mask

    total_counts = counts_df.sum(axis=1)
    nonzero_samples = counts_df.gt(0).sum(axis=1)
    count_mask = total_counts.ge(filter_params["min_intron_count"])
    sample_mask = nonzero_samples.ge(filter_params["min_intron_samples"])
    logger.info(
        "Intron count filter: "
        f"{int((count_mask & sample_mask).sum())}/{n_start} pass "
        f"total_count >= {filter_params['min_intron_count']} and "
        f"nonzero_samples >= {filter_params['min_intron_samples']}"
    )
    keep_mask &= count_mask & sample_mask

    offset_samples = offsets_df.ge(filter_params["min_offset_depth"]).sum(axis=1)
    offset_mask = offset_samples.ge(filter_params["min_offset_samples"])
    logger.info(
        "Denominator depth filter: "
        f"{int(offset_mask.sum())}/{n_start} pass "
        f"offset >= {filter_params['min_offset_depth']} in "
        f">= {filter_params['min_offset_samples']} samples"
    )
    keep_mask &= offset_mask

    psi_df = compute_site_depth_psi_values(
        counts_df,
        offsets_df,
        samples_df,
        edgeR_params["group_col"],
        edgeR_params["contrast"],
    )
    min_delta_psi = filter_params.get("min_delta_psi", 0)
    if min_delta_psi:
        delta_mask = psi_df["delta_PSI"].abs().ge(min_delta_psi)
        logger.info(
            "Pre-test delta PSI filter: "
            f"{int(delta_mask.sum())}/{n_start} pass |delta_PSI| >= {min_delta_psi}"
        )
        keep_mask &= delta_mask
    else:
        logger.info("Pre-test delta PSI filter disabled")

    kept_introns = keep_mask[keep_mask].index
    if len(kept_introns) == 0:
        raise ValueError("No introns passed the count/site-depth/delta-PSI filters")

    logger.info(
        f"Final pre-test filter: kept {len(kept_introns)}/{n_start} introns "
        f"({100.0 * len(kept_introns) / n_start:.1f}%)"
    )

    filtered_counts = counts_df.loc[kept_introns, sample_cols]
    filtered_offsets = offsets_df.loc[kept_introns, sample_cols]
    filtered_annotations = annotations_df.loc[kept_introns].copy()
    filtered_psi = psi_df.loc[kept_introns].copy()

    test_offsets = filtered_offsets
    filtered_annotations["offset_source"] = psi_denominator_label
    filtered_matrix = pd.concat([filtered_annotations, filtered_counts], axis=1)

    logger.info(f"Writing filtered intron matrix to {output_files['filtered_matrix']}")
    filtered_matrix.to_csv(output_files["filtered_matrix"], sep="\t", na_rep="NA")

    logger.info(f"Writing raw filtered denominator offsets to {output_files['raw_offsets']}")
    filtered_offsets.to_csv(output_files["raw_offsets"], sep="\t", na_rep="NA")

    logger.info(f"Writing pre-test PSI values to {output_files['psi']}")
    filtered_psi.to_csv(output_files["psi"], sep="\t", na_rep="NA")

    logger.info(f"Writing model count matrix to {output_files['counts']}")
    filtered_counts.to_csv(output_files["counts"], sep="\t", na_rep="NA")

    logger.info(f"Writing model log-offset matrix to {output_files['offsets']}")
    np.log(test_offsets + 0.5).to_csv(output_files["offsets"], sep="\t", na_rep="NA")

    logger.info(f"Writing model annotations to {output_files['annotations']}")
    filtered_annotations.to_csv(output_files["annotations"], sep="\t", na_rep="NA")

    logger.info(f"Writing filter parameter metadata to {output_files['params']}")
    with open(output_files["params"], "wt") as params_fh:
        json.dump(params_record, params_fh, indent=2, sort_keys=True)

    return output_files


def prepare_inputs_from_bam_manifest(
    samples_manifest,
    genome_fa,
    workdir,
    force_rerun=False,
    site_depth_window_radius=10,
    retained_depth_inner_offset=20,
    retained_depth_window_radius=None,
    min_mapping_quality=60,
    site_depth_strand_mode="unstranded",
    strict_overhang=6,
):
    """
    Count introns from BAMs and build count/depth denominator matrices.

    The manifest must contain:
        sample_type    replicate_id    bam_file

    Returns:
        Dict with matrix, depth denominator matrices, sample metadata, and BAM manifest paths.
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

    retained_depth_window_arg = ""
    if retained_depth_window_radius is not None:
        retained_depth_window_arg = (
            f"--retained_depth_window_radius {int(retained_depth_window_radius)} "
        )

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
            f"--retained_depth_inner_offset {int(retained_depth_inner_offset)} "
            f"{retained_depth_window_arg}"
            f"--min_mapping_quality {int(min_mapping_quality)} "
            f"--site_depth_strand_mode {shlex.quote(site_depth_strand_mode)} "
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

    logger.info("Counting target introns from BAMs (complete depth-column pass)")
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
            f"--retained_depth_inner_offset {int(retained_depth_inner_offset)} "
            f"{retained_depth_window_arg}"
            f"--min_mapping_quality {int(min_mapping_quality)} "
            f"--site_depth_strand_mode {shlex.quote(site_depth_strand_mode)} "
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
        "Building final intron count and depth matrices",
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

    depth_prefix = final_matrix[:-7] if final_matrix.endswith(".matrix") else final_matrix

    return {
        "matrix": final_matrix,
        "samples": downstream_samples,
        "bam_list": bam_list,
        "max_splice_plus_retained_depth": f"{depth_prefix}.max_splice_plus_retained_depth.matrix",
    }


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


def run_edgeR(edgeR_inputs, samples_file, output_dir, edgeR_params, force_rerun=False, cpu=1):
    """
    Run edgeR differential analysis.
    
    Args:
        edgeR_inputs: Dict with paths to counts, offsets, annotations
        samples_file: Sample metadata file
        output_dir: Output directory
        edgeR_params: Dict with edgeR parameters
        force_rerun: If True, rerun even if outputs exist
        cpu: Accepted for compatibility; not used by the single-contrast runner
        
    Returns:
        Path to intron results file
    """
    util_dir = os.path.join(os.path.dirname(__file__), "util")
    
    output_prefix = os.path.join(output_dir, "edgeR_results")
    intron_results_file = f"{output_prefix}.intron_results.tsv"
    params_file = f"{output_prefix}.params.json"
    
    contrast = edgeR_params.get("contrast")
    if not contrast:
        raise ValueError("--contrast is required; only individual contrasts are supported")

    params_record = {
        "edgeR_params": edgeR_params,
    }
    params_match = False
    if file_exists_and_valid(params_file):
        try:
            with open(params_file, "rt") as params_fh:
                params_match = json.load(params_fh) == params_record
        except Exception:
            params_match = False

    result_current = file_is_current(
        intron_results_file,
        [
            edgeR_inputs["counts"],
            edgeR_inputs["offsets"],
            edgeR_inputs["annotations"],
            samples_file,
        ],
    )
    if not force_rerun and result_current and params_match:
        logger.info("=== Running edgeR analysis ===")
        logger.info(f"SKIPPING - Results already exist: {intron_results_file}")
        return intron_results_file

    logger.info("=== Running edgeR analysis ===")
    logger.info(f"Single contrast: {contrast}")

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

    run_command(cmd, "Running edgeR analysis")

    with open(params_file, "wt") as params_fh:
        json.dump(params_record, params_fh, indent=2, sort_keys=True)
    
    return intron_results_file


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
    Add precomputed PSI values to model results.
    
    Writes intermediate PSI/model artifacts into output_dir and promotes only the
    primary user-facing result files to the parent directory.
    
    Args:
        intron_results_file: Path to intron results from the selected model
        psi_file: Path to PSI values file
        output_dir: Working directory for intermediate PSI-enhanced outputs
        min_delta_psi: Pre-test delta PSI threshold, used only for reporting
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

        final_sig_results = sig_source[sig_source['significant']].copy() if 'significant' in sig_source.columns else pd.DataFrame(columns=sig_source.columns)

        total_final_tests = len(sig_source)
        pct_significant = 100 * len(final_sig_results) / total_final_tests if total_final_tests else 0
        criteria = [f"FDR <= {fdr_threshold}", f"|logFC| >= {min_logFC}"]
        if min_delta_psi and 'delta_PSI' in sig_source.columns:
            criteria.append(f"pre-test |delta_PSI| >= {min_delta_psi}")

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
        "--offset_matrix",
        type=str,
        default=None,
        help="Precomputed intron x sample raw denominator matrix. Required in matrix mode.",
    )

    parser.add_argument(
        "--site_depth_window_radius",
        type=int,
        default=10,
        help="Bases in each exon-adjacent depth window",
    )

    parser.add_argument(
        "--retained_depth_inner_offset",
        type=int,
        default=20,
        help="Bases inside the intron to skip before retained-depth windows begin",
    )

    parser.add_argument(
        "--retained_depth_window_radius",
        type=int,
        default=None,
        help="Bases in each intron-interior retained-depth window; defaults to --site_depth_window_radius",
    )

    parser.add_argument(
        "--site_depth_min_mapq",
        type=int,
        default=60,
        help="Minimum mapping quality for reads counted in depth denominators",
    )

    parser.add_argument(
        "--site_depth_strand_mode",
        choices=["unstranded", "F", "R", "FR", "RF"],
        default="unstranded",
        help=(
            "Strand mode for BAM-manifest depth denominators. "
            "F/R are single-end modes; FR/RF are paired-end modes describing "
            "read1/read2 orientations relative to the transcript."
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
        "--min_offset_depth",
        type=int,
        default=20,
        help="Minimum selected denominator depth required in a sample",
    )

    parser.add_argument(
        "--min_offset_samples",
        type=int,
        default=3,
        help="Minimum samples meeting --min_offset_depth",
    )

    parser.add_argument(
        "--keep_noncanonical",
        action="store_true",
        help="Keep non-canonical splice sites",
    )
    
    # Statistical model parameters
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
        required=True,
        help="Single contrast to test (e.g., 'perturb,control' where log2FC = perturb/control)",
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
        help="Minimum absolute delta PSI required before statistical testing. Set to 0 to disable this prefilter.",
    )
    
    # Pipeline control
    parser.add_argument(
        "--cpu",
        type=int,
        default=1,
        help="Accepted for compatibility; the streamlined pipeline runs one contrast per invocation",
    )
    
    parser.add_argument(
        "--force_rerun",
        action="store_true",
        help="Force rerun of all steps even if output files exist (disables resume)",
    )

    parser.add_argument(
        "--offset_mode",
        choices=sorted(OFFSET_MODE_SPECS),
        default="splice_plus_retained",
        help=(
            "Depth mode for PSI denominators and statistical model inputs. "
            "splice_plus_retained uses max(splice-depth + retained intron-side depth) "
            "as both PSI denominator and edgeR exposure."
        ),
    )

    parser.add_argument(
        "--strict_overhang",
        type=int,
        default=6,
        help=argparse.SUPPRESS,
    )
    
    args = parser.parse_args()
    
    # Validate mutually exclusive options
    try:
        parse_contrast(args.contrast)
    except ValueError as exc:
        parser.error(str(exc))

    bam_input_mode = args.matrix is None
    if not args.matrix and not args.genome_fa:
        parser.error("--genome_fa is required when --matrix is omitted and BAMs are counted from --samples")
    if bam_input_mode and args.offset_matrix:
        parser.error(
            "When --matrix is omitted, depth matrices are generated from the BAMs listed in --samples; "
            "do not also provide --offset_matrix."
        )
    if not bam_input_mode and not args.offset_matrix:
        parser.error("Matrix mode requires --offset_matrix")
    if args.keep_noncanonical:
        parser.error("Noncanonical introns are not supported in the offset-mode refactor; omit --keep_noncanonical")
    
    # Create output directory and workdir for intermediates
    os.makedirs(args.output_dir, exist_ok=True)
    workdir = os.path.join(args.output_dir, "workdir")
    os.makedirs(workdir, exist_ok=True)

    matrix_file = args.matrix
    samples_file = args.samples
    offset_matrix_file = args.offset_matrix

    if bam_input_mode:
        logger.info("=== Preparing count and offset matrices from BAM manifest ===")
        prepared_inputs = prepare_inputs_from_bam_manifest(
            samples_manifest=args.samples,
            genome_fa=args.genome_fa,
            workdir=workdir,
            force_rerun=args.force_rerun,
            site_depth_window_radius=args.site_depth_window_radius,
            retained_depth_inner_offset=args.retained_depth_inner_offset,
            retained_depth_window_radius=args.retained_depth_window_radius,
            min_mapping_quality=args.site_depth_min_mapq,
            site_depth_strand_mode=args.site_depth_strand_mode,
            strict_overhang=args.strict_overhang,
        )
        matrix_file = prepared_inputs["matrix"]
        samples_file = prepared_inputs["samples"]
        offset_matrix_file = prepared_inputs[OFFSET_MODE_SPECS[args.offset_mode]["denominator_key"]]

    if not file_exists_and_valid(offset_matrix_file):
        parser.error(
            f"Offset denominator matrix not found or empty for mode {args.offset_mode}: "
            f"{offset_matrix_file}. If resuming an older BAM-mode run, rerun with --force_rerun."
        )
    
    # Resolve count / PSI-denominator / test-offset sources from the selected mode.
    mode_spec = OFFSET_MODE_SPECS[args.offset_mode]
    count_file = matrix_file
    psi_denominator_file = offset_matrix_file
    psi_denominator_label = mode_spec["psi_denominator_label"]
    logger.info("=== Differential Splicing Analysis Pipeline ===")
    logger.info(f"Input mode: {'BAM manifest' if bam_input_mode else 'matrix'}")
    logger.info(f"Input matrix: {matrix_file}")
    logger.info(f"Sample metadata: {samples_file}")
    if offset_matrix_file:
        logger.info(f"Offset denominator matrix: {offset_matrix_file}")
    logger.info(f"Offset mode: {args.offset_mode}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Work directory (intermediates): {workdir}")
    logger.info("Analysis mode: Intron-level with selected depth offsets")
    if bam_input_mode:
        logger.info(f"Site-depth strand mode: {args.site_depth_strand_mode}")
    
    if args.min_delta_psi:
        logger.info(f"Pre-test PSI filtering: |delta_PSI| >= {args.min_delta_psi}")
    else:
        logger.info("Pre-test PSI filtering: disabled")
    
    if args.force_rerun:
        logger.info("Force rerun enabled - will regenerate all outputs")
    else:
        logger.info("Resume mode enabled - will skip completed steps")
    
    # Prepare parameter dicts
    filter_params = {
        "min_intron_count": args.min_intron_count,
        "min_intron_samples": args.min_intron_samples,
        "min_offset_depth": args.min_offset_depth,
        "min_offset_samples": args.min_offset_samples,
        "min_delta_psi": args.min_delta_psi,
        "keep_noncanonical": args.keep_noncanonical,
    }
    
    edgeR_params = {
        "group_col": args.group_col,
        "batch_col": args.batch_col,
        "contrast": args.contrast,
        "fdr_threshold": args.fdr_threshold,
        "min_logFC": args.min_logFC,
    }
    
    # Step 1: Filter introns and prepare model inputs from selected denominator offsets.
    logger.info(f"\n{'='*60}")
    logger.info("Preparing model inputs")
    logger.info(f"{'='*60}\n")
    
    prepared_files = prepare_site_depth_edgeR_inputs(
        count_file,
        psi_denominator_file,
        samples_file=samples_file,
        output_dir=workdir,
        filter_params=filter_params,
        edgeR_params=edgeR_params,
        gtf_file=args.gtf,
        force_rerun=args.force_rerun,
        offset_mode_label=args.offset_mode,
        psi_denominator_label=psi_denominator_label,
        count_unit_label="read",
    )
    edgeR_inputs = {
        "counts": prepared_files["counts"],
        "offsets": prepared_files["offsets"],
        "raw_offsets": prepared_files["raw_offsets"],
        "annotations": prepared_files["annotations"],
    }
    
    # Step 2: Run statistical analysis
    logger.info(f"\n{'='*60}")
    logger.info("Running statistical analysis")
    logger.info(f"{'='*60}\n")

    intron_results = run_edgeR(
        edgeR_inputs, samples_file, workdir, edgeR_params,
        force_rerun=args.force_rerun,
        cpu=args.cpu
    )

    # Step 3: Add precomputed PSI values to results.
    logger.info(f"\n{'='*60}")
    logger.info("Adding PSI to results")
    logger.info(f"{'='*60}\n")
    
    intron_results_with_psi, significant_results = add_psi_and_filter(
        intron_results, prepared_files["psi"], workdir,  # Keep intermediates in workdir and promote final outputs
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
        logger.info(f"  - Significant introns: {sig_file}")
    
    logger.info(f"\nIntermediate files (in workdir/):")
    logger.info(f"  - Filtered introns: {workdir}/introns_filtered.tsv")
    logger.info(f"  - Filtered raw denominator offsets: {workdir}/site_depth_offsets.filtered.tsv")
    logger.info(f"  - Raw model results: {workdir}/edgeR_results.intron_results.tsv")
    logger.info(f"  - PSI values: {workdir}/psi.psi_values.tsv")
    diagnostics_pdf = f"{workdir}/edgeR_results.diagnostics.pdf"
    if file_exists_and_valid(diagnostics_pdf):
        logger.info(f"  - Diagnostics plots: {diagnostics_pdf}")


if __name__ == "__main__":
    main()
