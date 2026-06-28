#!/usr/bin/env python3

import os
import sys
import tempfile

import pysam


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "util")))
from site_depth import compute_site_depth_offsets, write_strand_partitioned_bams


def make_read(name, flag, start, length=20, mate_start=0, tlen=0, seq_base="A"):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = seq_base * length
    read.flag = flag
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigar = [(0, length)]
    read.next_reference_id = 0
    read.next_reference_start = mate_start
    read.template_length = tlen
    read.query_qualities = pysam.qualitystring_to_array("I" * length)
    return read


def build_test_bam(workdir):
    bam_path = os.path.join(workdir, "strand_depth.bam")
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 200}],
    }

    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        # Single-end forward read covering 1-based site 30.
        out.write(make_read("singleF", 0, 24))

        # Proper FR pair covering 1-based site 50 with overlapping mates.
        out.write(make_read("pairFR", 99, 44, mate_start=47, tlen=23))
        out.write(make_read("pairFR", 147, 47, mate_start=44, tlen=-23, seq_base="T"))

    pysam.index(bam_path)
    return bam_path


def main():
    with tempfile.TemporaryDirectory(prefix="dsf_site_depth_test_") as workdir:
        bam_path = build_test_bam(workdir)

        single_sites = {
            "chr1": {
                "plus_single": [(30, "+")],
                "minus_single": [(30, "-")],
            }
        }
        paired_sites = {
            "chr1": {
                "plus_pair": [(50, "+")],
                "minus_pair": [(50, "-")],
            }
        }
        unstranded_sites = {"chr1": {"unstranded_pair": [(50, "+")]}}

        f_offsets = compute_site_depth_offsets(
            bam_path, single_sites, window_radius=0, min_mapping_quality=0, strand_mode="F"
        )
        r_offsets = compute_site_depth_offsets(
            bam_path, single_sites, window_radius=0, min_mapping_quality=0, strand_mode="R"
        )
        fr_offsets = compute_site_depth_offsets(
            bam_path, paired_sites, window_radius=0, min_mapping_quality=0, strand_mode="FR"
        )
        rf_offsets = compute_site_depth_offsets(
            bam_path, paired_sites, window_radius=0, min_mapping_quality=0, strand_mode="RF"
        )
        unstranded_offsets = compute_site_depth_offsets(
            bam_path, unstranded_sites, window_radius=0, min_mapping_quality=0, strand_mode="unstranded"
        )

        assert f_offsets["plus_single"] == 1 and f_offsets["minus_single"] == 0, f_offsets
        assert r_offsets["plus_single"] == 0 and r_offsets["minus_single"] == 1, r_offsets
        assert fr_offsets["plus_pair"] == 1 and fr_offsets["minus_pair"] == 0, fr_offsets
        assert rf_offsets["plus_pair"] == 0 and rf_offsets["minus_pair"] == 1, rf_offsets
        assert unstranded_offsets["unstranded_pair"] == 1, unstranded_offsets

        prefix = os.path.join(workdir, "sample")
        write_strand_partitioned_bams(bam_path, prefix, "FR", min_mapping_quality=0)
        plus_bam = prefix + ".transcript_plus.bam"
        minus_bam = prefix + ".transcript_minus.bam"
        with pysam.AlignmentFile(plus_bam, "rb") as plus, pysam.AlignmentFile(minus_bam, "rb") as minus:
            assert sum(1 for _ in plus.fetch(until_eof=True)) == 2
            assert sum(1 for _ in minus.fetch(until_eof=True)) == 0
        assert os.path.exists(plus_bam + ".bai")
        assert os.path.exists(minus_bam + ".bai")

    print("site-depth strand tests passed")


if __name__ == "__main__":
    main()
