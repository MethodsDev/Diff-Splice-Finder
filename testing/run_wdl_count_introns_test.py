#!/usr/bin/env python3
"""Execute WDL/count_introns.wdl end-to-end with miniwdl against a tiny fixture.

Unlike run_bam_introns_test.py (which drives the Python script directly), this
test exercises the *workflow wrapper*: input localization, the BAM/FASTA index
co-location symlinks, the gzip pipe, the strand-BAM globs, and the output
wiring. It runs the containerized tool, so it needs `miniwdl` and a working
`docker` daemon; when either is missing the test skips (exit 0) rather than
failing, so it stays opt-in for environments without Docker.

The container image defaults to the published `latest` tag and can be
overridden with the DSF_WDL_TEST_IMAGE environment variable (e.g. to point at a
locally built tag in CI).
"""

import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
WDL_PATH = os.path.join(REPO_ROOT, "WDL", "count_introns.wdl")
FIXTURE_DIR = os.path.join(SCRIPT_DIR, "mode_refactor_inputs")

BAM = os.path.join(FIXTURE_DIR, "source_bams", "A1.bam")
BAM_INDEX = os.path.join(FIXTURE_DIR, "source_bams", "A1.bam.bai")
GENOME_FASTA = os.path.join(FIXTURE_DIR, "chrTiny.fa")
GENOME_FASTA_INDEX = os.path.join(FIXTURE_DIR, "chrTiny.fa.fai")

DEFAULT_IMAGE = "us-central1-docker.pkg.dev/methods-dev-lab/diff-splice-finder/diff-splice-finder:latest"
IMAGE = os.environ.get("DSF_WDL_TEST_IMAGE", DEFAULT_IMAGE)

SAMPLE_ID = "A1"
# Junction read counts are parameter-independent, so they are safe to pin.
EXPECTED_COUNTS = {
    "chrTiny:101-200": 8,
    "chrTiny:101-250": 4,
    "chrTiny:151-200": 2,
}
EXPECTED_HEADER_FIELDS = {
    "intron",
    "splice_pair",
    "splice_flag",
    "count",
    "site_depth_offset",
}


def skip(reason):
    print(f"WDL count_introns test SKIPPED: {reason}")
    sys.exit(0)


def check_prerequisites():
    if shutil.which("miniwdl") is None:
        skip("miniwdl not found on PATH")
    if shutil.which("docker") is None:
        skip("docker not found on PATH")
    try:
        subprocess.run(
            ["docker", "info"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=30,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError):
        skip("docker daemon not available")
    for path in (WDL_PATH, BAM, BAM_INDEX, GENOME_FASTA, GENOME_FASTA_INDEX):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing required fixture: {path}")


def run_workflow(run_dir, strand_mode):
    """Run the workflow and return the parsed outputs.json object."""
    inputs = {
        "CountIntrons.sample_id": SAMPLE_ID,
        "CountIntrons.bam_file": BAM,
        "CountIntrons.bam_index": BAM_INDEX,
        "CountIntrons.genome_fasta": GENOME_FASTA,
        "CountIntrons.genome_fasta_index": GENOME_FASTA_INDEX,
        "CountIntrons.site_depth_strand_mode": strand_mode,
        "CountIntrons.docker": IMAGE,
    }
    inputs_path = os.path.join(run_dir, "inputs.json")
    with open(inputs_path, "wt") as fh:
        json.dump(inputs, fh)

    result = subprocess.run(
        [
            "miniwdl",
            "run",
            WDL_PATH,
            "-i",
            inputs_path,
            "--dir",
            run_dir,
            "--no-color",
        ],
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"miniwdl run failed (strand_mode={strand_mode}) with exit "
            f"{result.returncode}:\n{result.stdout}"
        )

    run_subdirs = sorted(glob.glob(os.path.join(run_dir, "2*")))
    if not run_subdirs:
        raise RuntimeError(f"No miniwdl run directory produced under {run_dir}")
    outputs_path = os.path.join(run_subdirs[-1], "outputs.json")
    with open(outputs_path, "rt") as fh:
        return json.load(fh)


def read_introns(gz_path):
    import gzip

    with gzip.open(gz_path, "rt") as fh:
        lines = [line.rstrip("\n") for line in fh if line.strip()]
    header = lines[0].split("\t")
    rows = {}
    for line in lines[1:]:
        fields = line.split("\t")
        record = dict(zip(header, fields))
        rows[record["intron"]] = record
    return header, rows


def validate_unstranded(outputs):
    introns_gz = outputs["CountIntrons.introns"]
    if not os.path.exists(introns_gz):
        raise AssertionError(f"introns output missing: {introns_gz}")
    if not introns_gz.endswith(".gz"):
        raise AssertionError(f"introns output is not gzipped: {introns_gz}")

    header, rows = read_introns(introns_gz)
    missing_cols = EXPECTED_HEADER_FIELDS - set(header)
    if missing_cols:
        raise AssertionError(f"introns header missing columns: {sorted(missing_cols)}")

    for intron, expected_count in EXPECTED_COUNTS.items():
        if intron not in rows:
            raise AssertionError(f"expected intron {intron} not reported")
        record = rows[intron]
        if int(record["count"]) != expected_count:
            raise AssertionError(
                f"{intron}: expected count {expected_count}, got {record['count']}"
            )
        if record["splice_pair"] != "GT--AG":
            raise AssertionError(
                f"{intron}: expected GT--AG splice_pair, got {record['splice_pair']}"
            )
        if record["splice_flag"] != "OK":
            raise AssertionError(
                f"{intron}: expected OK splice_flag, got {record['splice_flag']}"
            )

    # Unstranded mode partitions nothing, so no strand BAMs should be emitted.
    if outputs["CountIntrons.strand_bams"]:
        raise AssertionError(
            "unstranded run unexpectedly produced strand_bams: "
            f"{outputs['CountIntrons.strand_bams']}"
        )
    print("  unstranded: introns.gz OK (counts, splice pairs, columns verified)")


def validate_stranded(outputs):
    strand_bams = outputs["CountIntrons.strand_bams"]
    strand_indexes = outputs["CountIntrons.strand_bam_indexes"]
    bam_names = sorted(os.path.basename(p) for p in strand_bams)
    index_names = sorted(os.path.basename(p) for p in strand_indexes)

    expected_bams = [f"{SAMPLE_ID}.transcript_minus.bam", f"{SAMPLE_ID}.transcript_plus.bam"]
    expected_indexes = [name + ".bai" for name in expected_bams]
    if bam_names != expected_bams:
        raise AssertionError(f"strand_bams mismatch: expected {expected_bams}, got {bam_names}")
    if index_names != expected_indexes:
        raise AssertionError(
            f"strand_bam_indexes mismatch: expected {expected_indexes}, got {index_names}"
        )
    for path in strand_bams + strand_indexes:
        if not os.path.exists(path):
            raise AssertionError(f"strand output missing on disk: {path}")
    if not os.path.exists(outputs["CountIntrons.introns"]):
        raise AssertionError("stranded run did not produce introns output")
    print("  stranded (F): introns.gz + transcript_{plus,minus}.bam(.bai) globs OK")


def main():
    check_prerequisites()
    print(f"Running WDL execution test with image: {IMAGE}")

    unstranded_dir = tempfile.mkdtemp(prefix="dsf_wdl_unstranded_")
    stranded_dir = tempfile.mkdtemp(prefix="dsf_wdl_stranded_")
    try:
        validate_unstranded(run_workflow(unstranded_dir, "unstranded"))
        validate_stranded(run_workflow(stranded_dir, "F"))
    finally:
        shutil.rmtree(unstranded_dir, ignore_errors=True)
        shutil.rmtree(stranded_dir, ignore_errors=True)

    print("WDL count_introns execution test passed")


if __name__ == "__main__":
    main()
