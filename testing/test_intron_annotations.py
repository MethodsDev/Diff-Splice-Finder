#!/usr/bin/env python3
"""Focused tests for exact GTF annotation and matrix-aligned artifacts."""

import csv
import importlib.util
import os
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UTIL = os.path.join(ROOT, "util")
sys.path.insert(0, UTIL)

import build_intron_annotations as builder  # noqa: E402
import integrate_results as annotations  # noqa: E402

DSF_SPEC = importlib.util.spec_from_file_location("dsf_main", os.path.join(ROOT, "DSF.py"))
DSF = importlib.util.module_from_spec(DSF_SPEC)
DSF_SPEC.loader.exec_module(DSF)


GTF_TEXT = """\
chr1\tTEST\tgene\t1\t1000\t.\t+\t.\tgene_id "A"; gene_name "Alpha";
chr1\tTEST\tgene\t50\t900\t.\t+\t.\tgene_id "B"; gene_name "Beta";
chr1\tTEST\texon\t100\t199\t.\t+\t.\tgene_id "B"; transcript_id "TB"; gene_name "Beta";
chr1\tTEST\texon\t301\t400\t.\t+\t.\tgene_id "B"; transcript_id "TB"; gene_name "Beta";
chr1\tTEST\texon\t2000\t2100\t.\t+\t.\tgene_id "C"; transcript_id "TC"; gene_name "TranscriptOnly";
chr1\tTEST\texon\t2200\t2300\t.\t+\t.\tgene_id "C"; transcript_id "TC"; gene_name "TranscriptOnly";
chr1\tTEST\tgene\t3000\t3100\t.\t+\t.\tgene_id "D"; gene_name "Boundary";
chr1\tTEST\tgene\t4000\t5000\t.\t+\t.\tgene_id "Z"; gene_name "Zeta";
chr1\tTEST\tgene\t4050\t4950\t.\t+\t.\tgene_id "AA"; gene_name "AlphaMulti";
"""


class AnnotationSemanticsTest(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.gtf = os.path.join(self.tempdir.name, "fixture.gtf")
        with open(self.gtf, "w") as gtf_fh:
            gtf_fh.write(GTF_TEXT)
        parsed = annotations.parse_gtf_file(self.gtf, use_cache=False)
        self.gene_map, self.known, self.tx_genes, self.tx_exons, self.membership = parsed
        self.index = annotations.build_gene_index(
            self.gene_map, self.tx_genes, self.tx_exons
        )

    def tearDown(self):
        self.tempdir.cleanup()

    def genes(self, coords):
        return annotations.find_overlapping_genes(coords, self.index, self.membership)

    def test_known_intron_uses_transcript_membership_before_spatial_overlap(self):
        coords = ("chr1", 200, 300, "+")
        self.assertIn(coords, self.known)
        self.assertEqual(self.genes(coords), ["Beta"])

    def test_inclusive_overlap_boundaries(self):
        self.assertEqual(self.genes(("chr1", 2990, 3000, "+")), ["Boundary"])
        self.assertEqual(self.genes(("chr1", 3100, 3110, "+")), ["Boundary"])

    def test_transcript_span_without_gene_feature_is_indexed(self):
        self.assertEqual(
            self.genes(("chr1", 2150, 2160, "+")), ["TranscriptOnly"]
        )

    def test_overlapping_gene_names_are_sorted_and_deduplicated(self):
        self.assertEqual(
            self.genes(("chr1", 4500, 4600, "+")), ["AlphaMulti", "Zeta"]
        )

    def test_legacy_cache_is_rebuilt_with_transcript_spans_and_memberships(self):
        cache_path = annotations.get_cached_gtf_filename(self.gtf)
        with open(cache_path, "w") as cache_fh:
            cache_fh.write(
                "## GTF Intron Cache File\n"
                "## Annotated Introns\n"
                "chr\tstart\tend\tstrand\tgene_names\n"
                "chr1\t200\t300\t+\tBeta\n"
                "## Gene Regions\n"
                "chr\tstart\tend\tstrand\tgene_id\tgene_name\n"
                "chr1\t1\t1000\t+\tA\tAlpha\n"
            )

        annotations.parse_gtf_file(self.gtf, use_cache=True)
        parsed = annotations.parse_gtf_file(self.gtf, use_cache=True)
        gene_map, known, tx_genes, tx_exons, membership = parsed
        index = annotations.build_gene_index(gene_map, tx_genes, tx_exons)

        self.assertEqual(
            annotations.find_overlapping_genes(
                ("chr1", 200, 300, "+"), index, membership
            ),
            ["Beta"],
        )
        self.assertEqual(
            annotations.find_overlapping_genes(
                ("chr1", 2150, 2160, "+"), index, membership
            ),
            ["TranscriptOnly"],
        )
        self.assertIn(("chr1", 200, 300, "+"), known)
        with open(cache_path) as cache_fh:
            self.assertIn("## Cache Schema Version: 2", cache_fh.read())

    def test_malformed_intron_id_is_unknown(self):
        rows = builder.build_annotations(["not-an-intron"], self.gtf)
        self.assertEqual(rows, [("not-an-intron", ".", "unknown", ".")])

    def test_builder_cli_preserves_matrix_order_and_does_not_create_gtf_cache(self):
        matrix = os.path.join(self.tempdir.name, "counts.tsv")
        output = os.path.join(self.tempdir.name, "annotations.tsv")
        intron_ids = [
            "chr1:4500-4600^GT--AG^OK",
            "not-an-intron",
            "chr1:200-300^GT--AG^OK",
        ]
        with open(matrix, "w") as matrix_fh:
            matrix_fh.write("intron_id\tS1\n")
            for intron_id in intron_ids:
                matrix_fh.write(f"{intron_id}\t1\n")

        subprocess.run(
            [
                sys.executable,
                os.path.join(UTIL, "build_intron_annotations.py"),
                "--matrix",
                matrix,
                "--gtf",
                self.gtf,
                "--output",
                output,
            ],
            check=True,
        )

        with open(output, newline="") as output_fh:
            rows = list(csv.DictReader(output_fh, delimiter="\t"))
        self.assertEqual([row["intron_id"] for row in rows], intron_ids)
        self.assertEqual(rows[0]["overlapping_genes"], "AlphaMulti,Zeta")
        self.assertEqual(rows[1]["intron_status"], "unknown")
        self.assertEqual(rows[2]["overlapping_genes"], "Beta")
        self.assertFalse(os.path.exists(annotations.get_cached_gtf_filename(self.gtf)))


class DSFAnnotationValidationTest(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.matrix_ids = ["intron-a", "intron-b"]

    def tearDown(self):
        self.tempdir.cleanup()

    def write_artifact(self, ids):
        path = os.path.join(self.tempdir.name, f"annotations-{len(os.listdir(self.tempdir.name))}.tsv")
        with open(path, "w") as output_fh:
            output_fh.write("intron_id\tgene_name\tintron_status\toverlapping_genes\n")
            for intron_id in ids:
                output_fh.write(f"{intron_id}\t.\tunknown\t.\n")
        return path

    def assert_rejected(self, ids, message):
        with self.assertRaisesRegex(ValueError, message):
            DSF.load_intron_annotations(self.write_artifact(ids), self.matrix_ids)

    def test_rejects_missing_annotation_id(self):
        self.assert_rejected(["intron-a"], "missing")

    def test_rejects_extra_annotation_id(self):
        self.assert_rejected(["intron-a", "intron-b", "intron-c"], "extra")

    def test_rejects_reordered_annotation_ids(self):
        self.assert_rejected(["intron-b", "intron-a"], "reordered")

    def test_rejects_duplicate_annotation_ids(self):
        self.assert_rejected(["intron-a", "intron-a"], "unique")


if __name__ == "__main__":
    unittest.main()
