#!/usr/bin/env python3
"""Strict splice-decision depth counter (fragment-level) for regular DSF.

This is a *denominator/offset* helper. It does not call differential usage; it
produces, per target intron and sample, a focal fragment count and a strict
local splice-decision depth, plus the per-boundary components used to build them.

Motivation (see docs/dsf_strict_depth_gene_median_plan.md in the benchmark repo):
the regular DSF site-depth offset is a max-coverage summary over a +/-10 bp window
and includes exon-body coverage that did not make a splice decision; it also
self-cancels when the focal junction dominates the site. This counter measures
only molecules that actually traverse a splice-site boundary.

For an intron j with left boundary at its donor coordinate `lend` and right
boundary at its acceptor coordinate `rend`, over fragments (independent
molecules, paired mates collapsed by query_name):

    cluster_split_B   = fragments with any split alignment supporting a junction
                        that shares boundary B (donor lend, or acceptor rend)
    unspliced_B       = fragments that read contiguously across boundary B with
                        >= overhang on both sides and DO NOT splice at B
    decision_depth_B  = cluster_split_B + unspliced_B
    focal_j           = fragments splicing exactly intron j

    strict_local_depth_j = max(decision_depth_left, decision_depth_right)

`max` (not sum) avoids double counting fragments that cross both boundaries of a
short intron and keeps the denominator independent of intron length, matching the
current DSF site-depth offset which is a max over its two windows.

Splice evidence dominates retention: a fragment that splices boundary B is never
also counted as unspliced at B (prevents retained-read inflation at splice sites).

Coordinate conventions (matched to DSF count_introns_from_bam / cluster_introns so
junction ids join 1:1): intron_id = "chrom:lend-rend^splice_pair^flag" with lend the
1-based first intron base and rend the 1-based last intron base. pysam find_introns
returns (start, stop) 0-based half-open, so the fragment junction key is
(chrom, start + 1, stop) == (chrom, lend, rend).

The fragment-collapse, overhang, and retention logic is ported from the vetted
DSF-beta counter (MethodsDev/DSF-beta, dsf_beta/count_fragments.py +
dsf_beta/build_tests.py), collapsed here to one row per intron for the DSF contract.
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
from bisect import bisect_left, bisect_right
from collections import defaultdict

import pandas as pd
import pysam

logger = logging.getLogger(__name__)

DEFAULT_MIN_MAPQ = 60
DEFAULT_OVERHANG = 6


def _passes(read, min_mapq: int) -> bool:
    if (read.is_unmapped or read.is_secondary or read.is_supplementary
            or read.is_duplicate or read.is_qcfail):
        return False
    return read.mapping_quality >= min_mapq


def _merge(blocks):
    """Merge overlapping/adjacent 0-based half-open intervals."""
    if not blocks:
        return []
    blocks = sorted(blocks)
    out = [list(blocks[0])]
    for bs, be in blocks[1:]:
        if bs <= out[-1][1]:
            if be > out[-1][1]:
                out[-1][1] = be
        else:
            out.append([bs, be])
    return [(s, e) for s, e in out]


def _positions_in(sorted_pos, lo, hi):
    """Sorted positions p with lo <= p <= hi."""
    if not sorted_pos or lo > hi:
        return ()
    i = bisect_left(sorted_pos, lo)
    j = bisect_right(sorted_pos, hi)
    return sorted_pos[i:j]


def count_bam(bam_path, target_lpos=None, target_rpos=None,
              min_mapq=DEFAULT_MIN_MAPQ, overhang=DEFAULT_OVERHANG, chroms=None):
    """Count fragment-level junctions and per-boundary retention for one BAM.

    target_lpos/target_rpos: dict[chrom -> set(pos)] of donor (lend) / acceptor (rend)
    boundaries to additionally compute retention for, even when no fragment splices
    them in this sample (so a target intron absent here still gets read-through depth).

    Returns (jcounts, retL, retR):
        jcounts: dict[(chrom, lend, rend) -> int]   fragment junction counts
        retL:    dict[(chrom, lend)       -> int]    read-through of the donor boundary
        retR:    dict[(chrom, rend)       -> int]    read-through of the acceptor boundary
    """
    af = pysam.AlignmentFile(bam_path, "rb")
    chrom_set = set(chroms) if chroms else None

    # Single collect pass: group mates by query_name into fragments.
    frag = {}  # qn -> [chrom, set((lend, rend)), [blocks]]
    for read in af.fetch():
        if not _passes(read, min_mapq):
            continue
        chrom = af.get_reference_name(read.reference_id)
        if chrom_set is not None and chrom not in chrom_set:
            continue
        rec = frag.get(read.query_name)
        if rec is None:
            rec = [chrom, set(), []]
            frag[read.query_name] = rec
        elif rec[0] != chrom:
            continue  # mate on a different contig: ignore the stray (rare)
        for (s, e) in af.find_introns([read]).keys():
            rec[1].add((s + 1, e))
        rec[2].extend(read.get_blocks())
    af.close()

    # Junctions + boundary universe (one fragment counts once per junction).
    jcounts = defaultdict(int)
    lpos = defaultdict(set)
    rpos = defaultdict(set)
    for chrom, juncs, _ in frag.values():
        for (lend, rend) in juncs:
            jcounts[(chrom, lend, rend)] += 1
            lpos[chrom].add(lend)
            rpos[chrom].add(rend)
    # Seed with target boundaries so targets unspliced in this sample still get retention.
    if target_lpos:
        for c, s in target_lpos.items():
            lpos[c].update(s)
    if target_rpos:
        for c, s in target_rpos.items():
            rpos[c].update(s)
    lsorted = {c: sorted(s) for c, s in lpos.items()}
    rsorted = {c: sorted(s) for c, s in rpos.items()}

    # Retention: fragments reading through a boundary without splicing it.
    retL = defaultdict(int)
    retR = defaultdict(int)
    for chrom, juncs, blocks in frag.values():
        if not blocks:
            continue
        merged = _merge(blocks)
        spliced_l = {lend for (lend, _r) in juncs}
        spliced_r = {rend for (_l, rend) in juncs}
        ls = lsorted.get(chrom)
        rs = rsorted.get(chrom)
        seen_l = set()
        seen_r = set()
        for (bs, be) in merged:
            if ls:
                # donor boundary 0-based base a = lend - 1; need a in [bs+ovh, be-ovh]
                for lend in _positions_in(ls, bs + overhang + 1, be - overhang + 1):
                    if lend not in spliced_l:
                        seen_l.add(lend)
            if rs:
                # acceptor boundary 0-based base a = rend; need rend in [bs+ovh, be-ovh]
                for rend in _positions_in(rs, bs + overhang, be - overhang):
                    if rend not in spliced_r:
                        seen_r.add(rend)
        for lend in seen_l:
            retL[(chrom, lend)] += 1
        for rend in seen_r:
            retR[(chrom, rend)] += 1

    return dict(jcounts), dict(retL), dict(retR)


def _source(left, right):
    if left > right:
        return "left"
    if right > left:
        return "right"
    return "tied"


def count_strict_depths(samples, bam_by_sample, targets,
                        min_mapq=DEFAULT_MIN_MAPQ, overhang=DEFAULT_OVERHANG, chroms=None):
    """Compute strict depths for every target intron across samples.

    samples:       ordered list of sample ids (output column order)
    bam_by_sample: dict sample_id -> bam path
    targets:       list of (intron_id, chrom, lend, rend)

    Returns a dict of pandas DataFrames, all indexed by intron_id in `targets` order
    with columns = samples:
        focal, strict_local, strict_left, strict_right,
        left_cluster, right_cluster, left_unspliced, right_unspliced
    plus 'source' (intron_id x samples, values in {left,right,tied}).
    """
    target_lpos = defaultdict(set)
    target_rpos = defaultdict(set)
    for _iid, chrom, lend, rend in targets:
        target_lpos[chrom].add(lend)
        target_rpos[chrom].add(rend)

    iids = [t[0] for t in targets]
    cols = list(samples)
    fields = ["focal", "strict_local", "strict_left", "strict_right",
              "left_cluster", "right_cluster", "left_unspliced", "right_unspliced", "source"]
    data = {f: pd.DataFrame(index=iids, columns=cols, dtype=object if f == "source" else int)
            for f in fields}

    for sid in samples:
        bam = bam_by_sample[sid]
        logger.info(f"[strict-depth] counting {sid} ({bam})")
        jcounts, retL, retR = count_bam(
            bam, target_lpos=target_lpos, target_rpos=target_rpos,
            min_mapq=min_mapq, overhang=overhang, chroms=chroms,
        )
        donor_tot = defaultdict(int)   # (chrom,lend) -> sum c over junctions sharing donor
        acc_tot = defaultdict(int)     # (chrom,rend) -> sum c over junctions sharing acceptor
        for (chrom, lend, rend), c in jcounts.items():
            donor_tot[(chrom, lend)] += c
            acc_tot[(chrom, rend)] += c

        focal_col, sl_col, slft, srt = [], [], [], []
        lc, rc, lu, ru, src = [], [], [], [], []
        for _iid, chrom, lend, rend in targets:
            focal = jcounts.get((chrom, lend, rend), 0)
            lcl = donor_tot.get((chrom, lend), 0)
            rcl = acc_tot.get((chrom, rend), 0)
            lun = retL.get((chrom, lend), 0)
            run = retR.get((chrom, rend), 0)
            dl = lcl + lun
            dr = rcl + run
            focal_col.append(focal)
            slft.append(dl)
            srt.append(dr)
            sl_col.append(dl if dl >= dr else dr)
            lc.append(lcl)
            rc.append(rcl)
            lu.append(lun)
            ru.append(run)
            src.append(_source(dl, dr))
        data["focal"][sid] = focal_col
        data["strict_local"][sid] = sl_col
        data["strict_left"][sid] = slft
        data["strict_right"][sid] = srt
        data["left_cluster"][sid] = lc
        data["right_cluster"][sid] = rc
        data["left_unspliced"][sid] = lu
        data["right_unspliced"][sid] = ru
        data["source"][sid] = src

    return data


def _load_targets(target_introns_file):
    """Read intron ids from a count matrix (index col) and parse to (id, chrom, lend, rend)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from cluster_introns import parse_intron_id

    df = pd.read_csv(target_introns_file, sep="\t", index_col=0)
    targets = []
    for iid in df.index.astype(str):
        info = parse_intron_id(iid)
        targets.append((iid, info["chr"], int(info["start"]), int(info["end"])))
    return targets


def _read_bam_list(bam_list_file):
    df = pd.read_csv(bam_list_file, sep="\t", comment="#")
    cols = {c.lower(): c for c in df.columns}
    sid_c = cols.get("sample_id") or cols.get("replicate_id")
    bam_c = cols.get("bam") or cols.get("bam_file")
    if not sid_c or not bam_c:
        raise ValueError("bam list must have sample_id/replicate_id and bam/bam_file columns")
    samples = df[sid_c].astype(str).tolist()
    return samples, dict(zip(df[sid_c].astype(str), df[bam_c].astype(str)))


def main():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s : %(levelname)s : %(message)s", datefmt="%H:%M:%S")
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam_list", required=True, help="TSV with sample_id/replicate_id and bam/bam_file")
    ap.add_argument("--target_introns", required=True,
                    help="Intron count matrix; its row index defines the target intron universe")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min_mapping_quality", type=int, default=DEFAULT_MIN_MAPQ)
    ap.add_argument("--overhang", type=int, default=DEFAULT_OVERHANG)
    ap.add_argument("--site_depth_strand_mode", default="unstranded",
                    help="Only 'unstranded' is supported; other modes fall back to genomic "
                         "(unstranded-equivalent) clustering with a warning.")
    ap.add_argument("--chroms", nargs="*", default=None)
    args = ap.parse_args()

    if args.site_depth_strand_mode != "unstranded":
        logger.warning("strict_splice_depth only implements unstranded counting; "
                       f"'{args.site_depth_strand_mode}' falls back to genomic clustering.")

    os.makedirs(args.outdir, exist_ok=True)
    samples, bam_by_sample = _read_bam_list(args.bam_list)
    targets = _load_targets(args.target_introns)
    logger.info(f"[strict-depth] {len(targets)} target introns x {len(samples)} samples; "
                f"overhang={args.overhang} min_mapq={args.min_mapping_quality}")

    data = count_strict_depths(samples, bam_by_sample, targets,
                               min_mapq=args.min_mapping_quality, overhang=args.overhang,
                               chroms=args.chroms)

    def _write(name, key):
        path = os.path.join(args.outdir, name)
        data[key].to_csv(path, sep="\t", na_rep="0")
        logger.info(f"[strict-depth] wrote {path}")
        return path

    _write("focal_fragment_counts.matrix", "focal")
    _write("strict_local_depth.matrix", "strict_local")
    _write("strict_left_depth.matrix", "strict_left")
    _write("strict_right_depth.matrix", "strict_right")

    # Wide audit TSV: focal + all components per sample, easy to join/inspect.
    wide = pd.DataFrame(index=data["focal"].index)
    wide.index.name = "intron_id"
    comp = {
        "focal_fragment_count": "focal",
        "strict_depth_left": "strict_left",
        "strict_depth_right": "strict_right",
        "strict_local_depth": "strict_local",
        "left_split_cluster_count": "left_cluster",
        "right_split_cluster_count": "right_cluster",
        "left_unspliced_count": "left_unspliced",
        "right_unspliced_count": "right_unspliced",
        "strict_depth_source": "source",
    }
    for sid in samples:
        for label, key in comp.items():
            wide[f"{sid}_{label}"] = data[key][sid].values
    audit = os.path.join(args.outdir, "intron_counts_with_strict_depths.tsv")
    wide.to_csv(audit, sep="\t")
    logger.info(f"[strict-depth] wrote {audit}")


if __name__ == "__main__":
    main()
