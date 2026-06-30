#!/usr/bin/env python3
"""Unit tests for util/strict_splice_depth.py on a synthetic BAM.

Covers the validation bullets from the plan:
  - a split fragment supporting intron j increments focal once;
  - a paired fragment whose mates both overlap j is not counted twice;
  - a fragment that splices boundary B is not also counted as unspliced at B;
  - a contiguous fragment crossing a boundary with overhang increments retention;
  - a contiguous fragment spanning both boundaries of a short intron does NOT double
    the intron-level depth (strict_local = max(left, right), not left + right);
  - a target intron unspliced in the sample still gets retention depth (target seeding);
  - focal_fragment_count <= strict_local_depth for every intron/sample.

Run: python testing/test_strict_splice_depth.py
"""
import os
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "util"))
import strict_splice_depth as ssd  # noqa: E402

CHROM = "chr1"
CHROM_LEN = 1000


def _seg(qname, start0, cigar, *, read1=True, paired=False):
    s = pysam.AlignedSegment()
    s.query_name = qname
    s.reference_id = 0
    s.reference_start = start0
    s.mapping_quality = 60
    s.cigartuples = cigar  # list of (op, length); 0=M, 3=N
    mlen = sum(length for op, length in cigar if op == 0)
    s.query_sequence = "A" * mlen
    s.query_qualities = [30] * mlen
    flag = 0
    if paired:
        flag |= 0x1 | 0x2
        flag |= 0x40 if read1 else 0x80
    s.flag = flag
    return s


def _build_bam(path):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}
    M, N = 0, 3
    reads = [
        # R1 focal: split intron j = (101,200): 20M100N20M @ 80 -> introns (100,200)0h -> (101,200)
        _seg("f_focal", 80, [(M, 20), (N, 100), (M, 20)]),
        # R2 paired: both mates split the SAME intron j -> fragment counts focal once
        _seg("f_pairfocal", 80, [(M, 20), (N, 100), (M, 20)], read1=True, paired=True),
        _seg("f_pairfocal", 79, [(M, 21), (N, 100), (M, 20)], read1=False, paired=True),
        # R3 alt junction sharing donor lend=101: 20M150N20M @ 80 -> (101,250)
        _seg("f_jalt", 80, [(M, 20), (N, 150), (M, 20)]),
        # R4 contiguous across acceptor boundary of j (rend=200), not spliced: 50M @ 180
        _seg("f_unsplicedR", 180, [(M, 50)]),
        # R5 contiguous spanning BOTH boundaries of short intron j2=(301,310): 40M @ 280
        _seg("f_short_both", 280, [(M, 40)]),
    ]
    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.sort("-o", path, unsorted)
    pysam.index(path)
    os.remove(unsorted)


def main():
    tmp = tempfile.mkdtemp(prefix="strict_depth_test_")
    bam = os.path.join(tmp, "s1.bam")
    _build_bam(bam)

    j = (101, 200)
    jalt = (101, 250)
    j2 = (301, 310)  # never spliced -> must come from targets (seeding test)
    targets = [
        (f"{CHROM}:{j[0]}-{j[1]}", CHROM, j[0], j[1]),
        (f"{CHROM}:{jalt[0]}-{jalt[1]}", CHROM, jalt[0], jalt[1]),
        (f"{CHROM}:{j2[0]}-{j2[1]}", CHROM, j2[0], j2[1]),
    ]
    samples = ["s1"]
    data = ssd.count_strict_depths(samples, {"s1": bam}, targets, min_mapq=60, overhang=6)

    def g(field, iid):
        return int(data[field].loc[iid, "s1"])

    jid = f"{CHROM}:{j[0]}-{j[1]}"
    jaltid = f"{CHROM}:{jalt[0]}-{jalt[1]}"
    j2id = f"{CHROM}:{j2[0]}-{j2[1]}"

    checks = []

    def check(name, cond, *vals):
        checks.append((name, cond))
        status = "ok" if cond else "FAIL"
        print(f"[{status}] {name}" + (f"  {vals}" if vals else ""))

    # focal: f_focal + f_pairfocal(once) = 2, NOT 3 (pair not double-counted)
    check("focal[j] == 2 (paired fragment counted once)", g("focal", jid) == 2, g("focal", jid))
    # donor cluster at lend=101 = j(2) + jalt(1) = 3
    check("left_cluster[j] == 3 (donor-shared junctions)", g("left_cluster", jid) == 3, g("left_cluster", jid))
    # acceptor retention at rend=200 from f_unsplicedR only (split reads excluded)
    check("right_unspliced[j] == 1 (split-dominates: spliced reads not retained)",
          g("right_unspliced", jid) == 1, g("right_unspliced", jid))
    check("left_unspliced[j] == 0", g("left_unspliced", jid) == 0, g("left_unspliced", jid))
    # strict_local[j] = max(dl=3, dr=cluster2+ret1=3) = 3
    check("strict_local[j] == 3", g("strict_local", jid) == 3, g("strict_local", jid))

    # jalt: dl = donor_tot[101]=3 ; dr = acc_tot[250]=1 -> max 3, source left
    check("strict_local[jalt] == 3", g("strict_local", jaltid) == 3, g("strict_local", jaltid))
    check("source[jalt] == left", data["source"].loc[jaltid, "s1"] == "left", data["source"].loc[jaltid, "s1"])

    # j2: never spliced (focal 0); one fragment crosses BOTH boundaries -> retL=retR=1
    # strict_local = max(1,1) = 1, NOT left+right = 2
    check("focal[j2] == 0 (target seeded, no splice)", g("focal", j2id) == 0, g("focal", j2id))
    check("left_unspliced[j2] == 1", g("left_unspliced", j2id) == 1, g("left_unspliced", j2id))
    check("right_unspliced[j2] == 1", g("right_unspliced", j2id) == 1, g("right_unspliced", j2id))
    check("strict_local[j2] == 1 (max not sum)", g("strict_local", j2id) == 1, g("strict_local", j2id))

    # global invariant: focal <= strict_local everywhere
    inv = all(g("focal", t[0]) <= g("strict_local", t[0]) for t in targets)
    check("focal <= strict_local for all targets", inv)

    failed = [n for n, c in checks if not c]
    if failed:
        print(f"\n{len(failed)} CHECK(S) FAILED: {failed}")
        sys.exit(1)
    print(f"\nAll {len(checks)} checks passed.")


if __name__ == "__main__":
    main()
