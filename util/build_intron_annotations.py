#!/usr/bin/env python3

"""Build a deterministic matrix-aligned intron annotation artifact."""

import argparse
import gzip
import os
import tempfile

from integrate_results import (
    build_gene_index,
    find_overlapping_genes,
    parse_gtf_file,
    parse_intron_id,
)


def iter_matrix_intron_ids(matrix_path):
    """Yield column one from a tab-delimited matrix in row order."""
    opener = gzip.open if matrix_path.endswith(".gz") else open
    with opener(matrix_path, "rt") as matrix_fh:
        header = matrix_fh.readline()
        if not header:
            raise ValueError(f"Matrix is empty: {matrix_path}")
        for line_number, line in enumerate(matrix_fh, start=2):
            intron_id = line.rstrip("\r\n").split("\t", 1)[0]
            if not intron_id:
                raise ValueError(f"Matrix row {line_number} has an empty intron ID")
            yield intron_id


def iter_annotations(intron_ids, gtf_path):
    """Yield raw-GTF annotations without creating a GTF sidecar cache."""
    (
        gene_map,
        annotated_introns,
        transcript_genes,
        transcript_exons,
        intron_gene_map,
    ) = parse_gtf_file(gtf_path, use_cache=False)
    gene_index = build_gene_index(gene_map, transcript_genes, transcript_exons)

    for intron_id in intron_ids:
        coords = parse_intron_id(intron_id)
        if coords is None:
            yield intron_id, ".", "unknown", "."
            continue

        status = "known" if coords in annotated_introns else "novel"
        genes = find_overlapping_genes(coords, gene_index, intron_gene_map)
        if genes:
            yield intron_id, genes[0], status, ",".join(genes)
        else:
            yield intron_id, ".", status, "."


def build_annotations(intron_ids, gtf_path):
    """Return annotations for callers that need an in-memory collection."""
    return list(iter_annotations(intron_ids, gtf_path))


def write_annotations_atomic(output_path, rows):
    """Publish a complete plain TSV with an atomic same-directory rename."""
    output_path = os.path.abspath(output_path)
    output_dir = os.path.dirname(output_path)
    os.makedirs(output_dir, exist_ok=True)
    temp_fd, temp_path = tempfile.mkstemp(
        prefix=f".{os.path.basename(output_path)}.", suffix=".tmp", dir=output_dir
    )
    try:
        with os.fdopen(temp_fd, "w") as output_fh:
            output_fh.write("intron_id\tgene_name\tintron_status\toverlapping_genes\n")
            for row in rows:
                output_fh.write("\t".join(row) + "\n")
            output_fh.flush()
            os.fsync(output_fh.fileno())
        os.replace(temp_path, output_path)
    except Exception:
        try:
            os.unlink(temp_path)
        except FileNotFoundError:
            pass
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Build exact GTF annotations in intron-matrix row order"
    )
    parser.add_argument("--matrix", required=True, help="Intron count matrix TSV or TSV.GZ")
    parser.add_argument("--gtf", required=True, help="Reference GTF annotation")
    parser.add_argument("--output", required=True, help="Output annotation TSV")
    args = parser.parse_args()

    rows = iter_annotations(iter_matrix_intron_ids(args.matrix), args.gtf)
    write_annotations_atomic(args.output, rows)


if __name__ == "__main__":
    main()
