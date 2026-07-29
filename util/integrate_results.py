#!/usr/bin/env python3

"""
Integrate results from donor and acceptor clustering analyses.

Creates unified summary files combining both analyses, allowing users to:
- See all tested introns across both clustering strategies
- Identify introns significant in donor, acceptor, or both
- Compare effect sizes between clustering methods
- Annotate introns with gene information from GTF (optional)
"""

import sys
import os
import argparse
import logging
import pandas as pd
import numpy as np
from collections import defaultdict
import re
import time
import tempfile

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def parse_intron_id(intron_id):
    """
    Parse intron ID to extract coordinates.
    
    Supports multiple formats:
    - chr:start-end:strand
    - chr:start-end^MOTIF^STATUS (extracts strand from motif)
    
    Args:
        intron_id: Intron identifier string
        
    Returns:
        Tuple of (chr, start, end, strand) or None if parsing fails
    """
    try:
        # Format 1: chr:start-end:strand
        match = re.match(r'^(.+):(\d+)-(\d+):([+-])$', intron_id)
        if match:
            chrom, start, end, strand = match.groups()
            return chrom, int(start), int(end), strand
        
        # Format 2: chr:start-end^MOTIF^STATUS
        # Motifs: GT--AG (plus strand), CT--AC (minus strand)
        match = re.match(r'^(.+):(\d+)-(\d+)\^([ACGT]{2})--([ACGT]{2})\^', intron_id)
        if match:
            chrom, start, end, donor, acceptor = match.groups()
            # Determine strand from splice motif
            # GT--AG is canonical for + strand
            # CT--AC is canonical for - strand (reverse complement of GT--AG)
            if donor == 'GT' and acceptor == 'AG':
                strand = '+'
            elif donor == 'CT' and acceptor == 'AC':
                strand = '-'
            else:
                # Non-canonical, try to infer or use +
                strand = '+' if donor == 'GT' else '-' if donor == 'CT' else '+'
            
            return chrom, int(start), int(end), strand
            
    except Exception as e:
        logger.warning(f"Failed to parse intron_id '{intron_id}': {e}")
    return None


GTF_CACHE_SCHEMA_VERSION = 2


class GtfCacheFormatError(ValueError):
    """Raised when a GTF cache cannot preserve current annotation semantics."""


class _IntervalNode:
    """Node in a height-balanced interval tree augmented with subtree max-end."""

    __slots__ = ("start", "end", "gene_name", "max_end", "left", "right")

    def __init__(self, intervals, low, high):
        middle = (low + high) // 2
        self.start, self.end, self.gene_name = intervals[middle]
        self.left = _IntervalNode(intervals, low, middle) if low < middle else None
        self.right = _IntervalNode(intervals, middle + 1, high) if middle + 1 < high else None
        self.max_end = max(
            self.end,
            self.left.max_end if self.left is not None else self.end,
            self.right.max_end if self.right is not None else self.end,
        )


class IntervalTree:
    """Immutable height-balanced interval tree for named, closed intervals."""

    __slots__ = ("root", "size")

    def __init__(self, intervals):
        unique_intervals = intervals if isinstance(intervals, set) else set(intervals)
        deduplicated = tuple(sorted(unique_intervals))
        self.size = len(deduplicated)
        self.root = (
            _IntervalNode(deduplicated, 0, len(deduplicated))
            if deduplicated
            else None
        )

    def query(self, start, end):
        results = set()
        if self.root is None:
            return results

        pending = [self.root]
        while pending:
            node = pending.pop()
            if node.left is not None and node.left.max_end >= start:
                pending.append(node.left)
            if node.start <= end and node.end >= start:
                results.add(node.gene_name)
            if (
                node.right is not None
                and node.start <= end
                and node.right.max_end >= start
            ):
                pending.append(node.right)
        return results


def get_cached_gtf_filename(gtf_file):
    """Generate the sidecar cache filename for a GTF annotation file."""
    base = os.path.splitext(gtf_file)[0]
    return f"{base}.intron_cache.tsv"


def _transcript_regions(transcript_genes, transcript_exons):
    """Yield deterministic transcript spans and exact gene memberships."""
    for transcript_id in sorted(transcript_exons):
        exons = transcript_exons[transcript_id]
        gene_info = transcript_genes.get(transcript_id)
        if not exons or not gene_info:
            continue
        chrom = gene_info['chr']
        strand = gene_info['strand']
        starts = [exon[1] for exon in exons]
        ends = [exon[2] for exon in exons]
        yield (
            chrom,
            min(starts),
            max(ends),
            strand,
            transcript_id,
            gene_info['gene_id'],
            gene_info['gene_name'],
        )


def write_gtf_cache(cache_file, gene_map, annotated_introns, transcript_genes, transcript_exons):
    """Atomically write a complete, deterministic version-2 GTF cache."""
    logger.info(f"Writing GTF cache to: {cache_file}")

    intron_genes = defaultdict(set)
    for transcript_id, exons in transcript_exons.items():
        if len(exons) < 2:
            continue
        gene_info = transcript_genes.get(transcript_id)
        if not gene_info:
            continue
        gene_name = gene_info['gene_name']
        exons_sorted = sorted(exons, key=lambda exon: exon[1])
        for first, second in zip(exons_sorted, exons_sorted[1:]):
            chrom1, _, end1, strand = first
            chrom2, start2, _, _ = second
            intron_start = end1 + 1
            intron_end = start2 - 1
            if chrom1 == chrom2 and intron_start < intron_end:
                intron_genes[(chrom1, intron_start, intron_end, strand)].add(gene_name)

    cache_dir = os.path.dirname(os.path.abspath(cache_file))
    temp_fd, temp_path = tempfile.mkstemp(prefix=".gtf-cache-", suffix=".tmp", dir=cache_dir)
    try:
        with os.fdopen(temp_fd, 'w') as cache_fh:
            cache_fh.write("## GTF Intron Cache File\n")
            cache_fh.write(f"## Cache Schema Version: {GTF_CACHE_SCHEMA_VERSION}\n")
            cache_fh.write("##\n## Gene Regions\n")
            cache_fh.write("chr\tstart\tend\tstrand\tgene_id\tgene_name\n")
            for gene_coords, gene_info in sorted(gene_map.items()):
                chrom, start, end, strand = gene_coords
                cache_fh.write(
                    f"{chrom}\t{start}\t{end}\t{strand}\t"
                    f"{gene_info['gene_id']}\t{gene_info['gene_name']}\n"
                )

            cache_fh.write("##\n## Transcript Regions\n")
            cache_fh.write("chr\tstart\tend\tstrand\ttranscript_id\tgene_id\tgene_name\n")
            for region in _transcript_regions(transcript_genes, transcript_exons):
                cache_fh.write("\t".join(map(str, region)) + "\n")

            cache_fh.write("##\n## Annotated Introns\n")
            cache_fh.write("chr\tstart\tend\tstrand\tgene_names\n")
            for intron_coords in sorted(annotated_introns):
                chrom, start, end, strand = intron_coords
                genes = sorted(intron_genes.get(intron_coords, set()))
                gene_str = ','.join(genes) if genes else '.'
                cache_fh.write(f"{chrom}\t{start}\t{end}\t{strand}\t{gene_str}\n")
            cache_fh.flush()
            os.fsync(cache_fh.fileno())
        os.replace(temp_path, cache_file)
    except Exception:
        try:
            os.unlink(temp_path)
        except FileNotFoundError:
            pass
        raise

    logger.info(
        f"Cached {len(annotated_introns)} introns, {len(gene_map)} genes, "
        f"and {len(transcript_exons)} transcripts"
    )


def read_gtf_cache(cache_file):
    """Read a complete version-2 cache, rejecting legacy/incomplete formats."""
    logger.info(f"Reading GTF cache from: {cache_file}")
    gene_map = {}
    annotated_introns = set()
    transcript_genes = {}
    transcript_exons = defaultdict(list)
    intron_gene_map = {}
    sections_seen = set()
    section = None
    schema_version = None

    with open(cache_file, 'r') as cache_fh:
        for raw_line in cache_fh:
            line = raw_line.rstrip('\n')
            if line.startswith("## Cache Schema Version:"):
                try:
                    schema_version = int(line.rsplit(':', 1)[1].strip())
                except ValueError as exc:
                    raise GtfCacheFormatError("Invalid GTF cache schema version") from exc
                continue
            if line.startswith('##'):
                if 'Gene Regions' in line:
                    section = 'genes'
                    sections_seen.add(section)
                elif 'Transcript Regions' in line:
                    section = 'transcripts'
                    sections_seen.add(section)
                elif 'Annotated Introns' in line:
                    section = 'introns'
                    sections_seen.add(section)
                continue
            if not line or line in {
                'chr\tstart\tend\tstrand\tgene_id\tgene_name',
                'chr\tstart\tend\tstrand\ttranscript_id\tgene_id\tgene_name',
                'chr\tstart\tend\tstrand\tgene_names',
            }:
                continue

            fields = line.split('\t')
            try:
                if section == 'genes' and len(fields) == 6:
                    chrom, start_text, end_text, strand, gene_id, gene_name = fields
                    start = int(start_text)
                    end = int(end_text)
                    gene_map[(chrom, start, end, strand)] = {
                        'gene_id': gene_id,
                        'gene_name': gene_name,
                        'chr': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                    }
                elif section == 'transcripts' and len(fields) == 7:
                    chrom, start_text, end_text, strand, transcript_id, gene_id, gene_name = fields
                    start = int(start_text)
                    end = int(end_text)
                    transcript_genes[transcript_id] = {
                        'gene_id': gene_id,
                        'gene_name': gene_name,
                        'chr': chrom,
                        'strand': strand,
                    }
                    transcript_exons[transcript_id].append((chrom, start, end, strand))
                elif section == 'introns' and len(fields) == 5:
                    chrom, start_text, end_text, strand, gene_str = fields
                    intron_coords = (chrom, int(start_text), int(end_text), strand)
                    annotated_introns.add(intron_coords)
                    intron_gene_map[intron_coords] = (
                        set(gene_str.split(',')) if gene_str != '.' else set()
                    )
                else:
                    raise GtfCacheFormatError(
                        f"Malformed row in GTF cache section {section!r}: {line}"
                    )
            except ValueError as exc:
                if isinstance(exc, GtfCacheFormatError):
                    raise
                raise GtfCacheFormatError(f"Malformed coordinates in GTF cache: {line}") from exc

    if schema_version != GTF_CACHE_SCHEMA_VERSION:
        raise GtfCacheFormatError(
            f"GTF cache schema {schema_version!r} is not supported; "
            f"expected {GTF_CACHE_SCHEMA_VERSION}"
        )
    required_sections = {'genes', 'transcripts', 'introns'}
    if sections_seen != required_sections:
        raise GtfCacheFormatError(
            "GTF cache is incomplete; missing sections: "
            + ", ".join(sorted(required_sections - sections_seen))
        )

    logger.info(
        f"Loaded {len(annotated_introns)} introns, {len(gene_map)} genes, "
        f"and {len(transcript_exons)} transcripts from cache"
    )
    return gene_map, annotated_introns, transcript_genes, transcript_exons, intron_gene_map


def parse_gtf_file(gtf_file, use_cache=True):
    """Parse a GTF, optionally using/writing its complete versioned sidecar cache."""
    cache_file = get_cached_gtf_filename(gtf_file)
    if use_cache and os.path.exists(cache_file):
        gtf_mtime = os.path.getmtime(gtf_file)
        cache_mtime = os.path.getmtime(cache_file)
        if cache_mtime >= gtf_mtime:
            try:
                logger.info(f"Using cached GTF data: {cache_file}")
                return read_gtf_cache(cache_file)
            except GtfCacheFormatError as exc:
                logger.info(f"Rebuilding incompatible GTF cache: {exc}")
        else:
            logger.info("Cache file outdated, re-parsing GTF")

    logger.info(f"Parsing GTF file: {gtf_file}")
    gene_map = {}
    transcript_exons = defaultdict(list)
    transcript_genes = {}

    with open(gtf_file, 'r') as gtf_fh:
        for line in gtf_fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            chrom = fields[0]
            feature_type = fields[2]
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue
            strand = fields[6]

            attr_dict = {}
            for attr in fields[8].split(';'):
                match = re.match(r'(\S+)\s+"([^"]+)"', attr.strip())
                if match:
                    attr_dict[match.group(1)] = match.group(2)
            gene_id = attr_dict.get('gene_id', '')
            gene_name = attr_dict.get('gene_name', gene_id)
            transcript_id = attr_dict.get('transcript_id', '')

            if feature_type == 'exon' and transcript_id:
                transcript_exons[transcript_id].append((chrom, start, end, strand))
                if transcript_id not in transcript_genes:
                    transcript_genes[transcript_id] = {
                        'gene_id': gene_id,
                        'gene_name': gene_name,
                        'chr': chrom,
                        'strand': strand,
                    }
            elif feature_type == 'gene':
                gene_map[(chrom, start, end, strand)] = {
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'chr': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                }

    annotated_introns = set()
    intron_gene_map = defaultdict(set)
    for transcript_id, exons in transcript_exons.items():
        if len(exons) < 2:
            continue
        gene_info = transcript_genes.get(transcript_id)
        if not gene_info:
            continue
        exons_sorted = sorted(exons, key=lambda exon: exon[1])
        for first, second in zip(exons_sorted, exons_sorted[1:]):
            chrom1, _, end1, strand = first
            chrom2, start2, _, _ = second
            intron_start = end1 + 1
            intron_end = start2 - 1
            if chrom1 == chrom2 and intron_start < intron_end:
                intron_coords = (chrom1, intron_start, intron_end, strand)
                annotated_introns.add(intron_coords)
                intron_gene_map[intron_coords].add(gene_info['gene_name'])

    logger.info(f"Parsed {len(gene_map)} genes and {len(annotated_introns)} annotated introns")
    if use_cache:
        try:
            write_gtf_cache(cache_file, gene_map, annotated_introns, transcript_genes, transcript_exons)
        except Exception as exc:
            logger.warning(f"Failed to write cache file: {exc}")
    return gene_map, annotated_introns, transcript_genes, transcript_exons, intron_gene_map


def build_gene_index(gene_map, transcript_genes, transcript_exons):
    """Build a balanced interval tree for each chromosome/strand."""
    logger.info("Building interval-tree index for genes...")
    regions = defaultdict(set)
    for gene_coords, gene_info in gene_map.items():
        chrom, start, end, strand = gene_coords
        regions[(chrom, strand)].add((start, end, gene_info['gene_name']))

    for chrom, start, end, strand, _, _, gene_name in _transcript_regions(
        transcript_genes, transcript_exons
    ):
        regions[(chrom, strand)].add((start, end, gene_name))

    index = {key: IntervalTree(intervals) for key, intervals in regions.items()}
    total_entries = sum(tree.size for tree in index.values())
    logger.info(
        f"Built interval-tree index with {len(index)} chromosome/strand combinations, "
        f"{total_entries} gene regions"
    )
    return index


def find_overlapping_genes(intron_coords, gene_index, intron_gene_map=None):
    """Return sorted genes using exact known membership before interval overlap."""
    if intron_gene_map is not None and intron_coords in intron_gene_map:
        return sorted(intron_gene_map[intron_coords])

    chrom, start, end, strand = intron_coords
    tree = gene_index.get((chrom, strand))
    if tree is None:
        return []
    return sorted(tree.query(start, end))


def annotate_introns_with_gtf(integrated_df, gtf_file):
    """
    Annotate introns with gene information from GTF file.
    
    Args:
        integrated_df: Integrated results DataFrame
        gtf_file: Path to GTF annotation file
        
    Returns:
        Annotated DataFrame with additional columns:
        - gene_name: Best matching gene
        - overlapping_genes: All overlapping genes (comma-separated)
        - intron_status: 'known' or 'novel'
    """
    # Check if already annotated (from earlier pipeline steps)
    if 'gene_name' in integrated_df.columns and 'intron_status' in integrated_df.columns:
        logger.info("Introns already annotated with gene information, skipping GTF annotation")
        return integrated_df
    
    logger.info("Annotating introns with GTF information...")
    
    # Parse GTF (uses cache if available)
    gene_map, annotated_introns, transcript_genes, transcript_exons, intron_gene_map = parse_gtf_file(gtf_file)
    
    # Build spatial index for fast gene lookups
    gene_index = build_gene_index(gene_map, transcript_genes, transcript_exons)
    
    # Annotate each intron
    gene_names = []
    overlapping_genes_list = []
    intron_statuses = []
    
    total_introns = len(integrated_df)
    logger.info(f"Annotating {total_introns} introns...")
    
    # Adaptive progress reporting based on dataset size
    if total_introns < 1000:
        report_every = 100
    elif total_introns < 10000:
        report_every = 1000
    elif total_introns < 100000:
        report_every = 5000
    else:
        report_every = 10000
    
    start_time = time.time()
    last_report_time = start_time
    
    for idx, intron_id in enumerate(integrated_df['intron_id']):
        # Progress indicator with timing
        if (idx + 1) % report_every == 0 or (idx + 1) == total_introns:
            elapsed = time.time() - start_time
            current_time = time.time()
            interval = current_time - last_report_time
            last_report_time = current_time
            
            rate = (idx + 1) / elapsed if elapsed > 0 else 0
            remaining = (total_introns - idx - 1) / rate if rate > 0 else 0
            
            logger.info(
                f"  Processed {idx + 1:,}/{total_introns:,} introns "
                f"({100*(idx+1)/total_introns:.1f}%) | "
                f"Rate: {rate:.0f} introns/sec | "
                f"Elapsed: {elapsed:.1f}s | "
                f"Remaining: ~{remaining:.0f}s"
            )
        
        coords = parse_intron_id(intron_id)
        
        if coords is None:
            gene_names.append('.')
            overlapping_genes_list.append('.')
            intron_statuses.append('unknown')
            continue
        
        # Check if intron is known (exact match)
        is_known = coords in annotated_introns
        intron_statuses.append('known' if is_known else 'novel')
        
        # Find overlapping genes (uses spatial index and cache)
        overlapping = find_overlapping_genes(coords, gene_index, intron_gene_map)
        
        if overlapping:
            gene_names.append(overlapping[0])  # Best match (first one)
            overlapping_genes_list.append(','.join(overlapping))
        else:
            gene_names.append('.')
            overlapping_genes_list.append('.')
    
    # Final timing summary
    total_time = time.time() - start_time
    avg_rate = total_introns / total_time if total_time > 0 else 0
    logger.info(f"Annotation completed in {total_time:.1f}s (average rate: {avg_rate:.0f} introns/sec)")
    
    # Add columns to dataframe (insert right after intron_id)
    integrated_df = integrated_df.copy()
    
    # Get column order - insert annotation columns after intron_id
    cols = list(integrated_df.columns)
    intron_id_idx = cols.index('intron_id')
    
    # Insert new columns after intron_id
    cols.insert(intron_id_idx + 1, 'gene_name')
    cols.insert(intron_id_idx + 2, 'intron_status')
    cols.insert(intron_id_idx + 3, 'overlapping_genes')
    
    # Add the data
    integrated_df['gene_name'] = gene_names
    integrated_df['overlapping_genes'] = overlapping_genes_list
    integrated_df['intron_status'] = intron_statuses
    
    # Reorder columns
    integrated_df = integrated_df[cols]
    
    # Log summary statistics
    try:
        known_count = (integrated_df['intron_status'] == 'known').sum()
        novel_count = (integrated_df['intron_status'] == 'novel').sum()
        unknown_count = (integrated_df['intron_status'] == 'unknown').sum()
        total_introns = len(integrated_df)
        
        # Convert to Python int if it's a numpy/pandas scalar
        if hasattr(known_count, 'item'):
            known_count = known_count.item()
        if hasattr(novel_count, 'item'):
            novel_count = novel_count.item()
        if hasattr(unknown_count, 'item'):
            unknown_count = unknown_count.item()
        
        logger.info(f"Annotation summary:")
        logger.info(f"  Known introns: {known_count} ({100*known_count/total_introns:.1f}%)")
        logger.info(f"  Novel introns: {novel_count} ({100*novel_count/total_introns:.1f}%)")
        if unknown_count > 0:
            logger.info(f"  Unknown (failed parsing): {unknown_count}")
        
        with_genes = (integrated_df['gene_name'] != '.').sum()
        if hasattr(with_genes, 'item'):
            with_genes = with_genes.item()
        logger.info(f"  Introns with gene assignment: {with_genes} ({100*with_genes/total_introns:.1f}%)")
    except Exception as e:
        logger.warning(f"Could not generate annotation summary: {e}")
    
    return integrated_df


def load_results(donor_file, acceptor_file):
    """
    Load donor and acceptor results files.
    
    Args:
        donor_file: Path to donor intron results
        acceptor_file: Path to acceptor intron results
        
    Returns:
        Tuple of (donor_df, acceptor_df)
    """
    logger.info(f"Loading donor results from {donor_file}")
    donor_df = pd.read_csv(donor_file, sep="\t")
    
    logger.info(f"Loading acceptor results from {acceptor_file}")
    acceptor_df = pd.read_csv(acceptor_file, sep="\t")
    
    logger.info(f"Loaded {len(donor_df)} donor introns, {len(acceptor_df)} acceptor introns")
    
    return donor_df, acceptor_df


def integrate_intron_results(donor_df, acceptor_df):
    """
    Integrate donor and acceptor intron-level results.
    
    For large datasets with multiple contrasts, integrates per-contrast to reduce memory usage.
    
    Args:
        donor_df: DataFrame with donor results
        acceptor_df: DataFrame with acceptor results
        
    Returns:
        Integrated DataFrame with columns indicating significance in each analysis
    """
    logger.info("Integrating intron-level results...")
    
    # Check if we have contrast columns (multiple contrasts per intron)
    has_contrasts = 'contrast' in donor_df.columns or 'contrast' in acceptor_df.columns
    
    if has_contrasts:
        # Get unique contrasts from both datasets
        donor_contrasts = set(donor_df['contrast'].unique()) if 'contrast' in donor_df.columns else set()
        acceptor_contrasts = set(acceptor_df['contrast'].unique()) if 'contrast' in acceptor_df.columns else set()
        all_contrasts = sorted(donor_contrasts | acceptor_contrasts)
        
        logger.info(f"Found {len(all_contrasts)} contrasts - integrating per-contrast to reduce memory usage")
        
        # Integrate each contrast separately and combine
        integrated_chunks = []
        for i, contrast in enumerate(all_contrasts, 1):
            if i % 5 == 0 or i == len(all_contrasts):
                logger.info(f"  Processing contrast {i}/{len(all_contrasts)}: {contrast}")
            
            # Filter to single contrast
            donor_subset = donor_df[donor_df['contrast'] == contrast].copy() if 'contrast' in donor_df.columns else pd.DataFrame()
            acceptor_subset = acceptor_df[acceptor_df['contrast'] == contrast].copy() if 'contrast' in acceptor_df.columns else pd.DataFrame()
            
            # Integrate this contrast
            if len(donor_subset) > 0 and len(acceptor_subset) > 0:
                chunk = integrate_single_contrast(donor_subset, acceptor_subset)
            elif len(donor_subset) > 0:
                chunk = integrate_single_contrast(donor_subset, pd.DataFrame())
            elif len(acceptor_subset) > 0:
                chunk = integrate_single_contrast(pd.DataFrame(), acceptor_subset)
            else:
                continue
            
            integrated_chunks.append(chunk)
        
        logger.info(f"Combining {len(integrated_chunks)} integrated contrasts...")
        integrated = pd.concat(integrated_chunks, ignore_index=True)
        logger.info(f"Combined into {len(integrated)} total rows")
    else:
        # Single contrast - integrate directly
        integrated = integrate_single_contrast(donor_df, acceptor_df)
    
    # Summary statistics for combined data
    logger.info(f"\n=== Integration Summary ===")
    logger.info(f"Total unique introns: {len(integrated)}")
    logger.info(f"Tested in both analyses: {(integrated['tested_in'] == 'both').sum()}")
    logger.info(f"Tested in donor only: {(integrated['tested_in'] == 'donor_only').sum()}")
    logger.info(f"Tested in acceptor only: {(integrated['tested_in'] == 'acceptor_only').sum()}")
    
    sig_counts = integrated['significant_in'].value_counts()
    logger.info(f"\nSignificance summary:")
    for status, count in sig_counts.items():
        logger.info(f"  {status}: {count}")
    
    # Direction consistency for those significant in both
    both_sig_count = (integrated['significant_in'] == 'both').sum()
    if both_sig_count > 0:
        both_sig_mask = integrated['significant_in'] == 'both'
        consistent = integrated.loc[both_sig_mask, 'direction_consistent'].sum()
        logger.info(f"\nDirection consistency (significant in both):")
        logger.info(f"  Consistent: {consistent}/{both_sig_count} ({100*consistent/both_sig_count:.1f}%)")
    
    return integrated


def integrate_single_contrast(donor_df, acceptor_df):
    """
    Integrate donor and acceptor results for a single contrast.
    
    Args:
        donor_df: DataFrame with donor results (single contrast)
        acceptor_df: DataFrame with acceptor results (single contrast)
        
    Returns:
        Integrated DataFrame
    """
    
    # Add clustering type indicator
    donor_df = donor_df.copy()
    acceptor_df = acceptor_df.copy()
    
    # Preserve contrast column if it exists
    has_contrast = 'contrast' in donor_df.columns or 'contrast' in acceptor_df.columns
    
    # Rename cluster-specific columns
    donor_df = donor_df.rename(columns={
        'logFC': 'donor_logFC',
        'logCPM': 'donor_logCPM',
        'F': 'donor_F',
        'PValue': 'donor_PValue',
        'FDR': 'donor_FDR',
        'FDR_original': 'donor_FDR_original',
        'cluster': 'donor_cluster',
        'significant': 'donor_significant',
        'contrast': 'donor_contrast'  # Preserve contrast if present
    })
    
    # Also rename PSI columns to be donor-specific
    for col in donor_df.columns:
        if '_PSI' in col and not col.startswith('donor_'):
            donor_df = donor_df.rename(columns={col: f'donor_{col}'})
    
    acceptor_df = acceptor_df.rename(columns={
        'logFC': 'acceptor_logFC',
        'logCPM': 'acceptor_logCPM',
        'F': 'acceptor_F',
        'PValue': 'acceptor_PValue',
        'FDR': 'acceptor_FDR',
        'FDR_original': 'acceptor_FDR_original',
        'cluster': 'acceptor_cluster',
        'significant': 'acceptor_significant',
        'contrast': 'acceptor_contrast'  # Preserve contrast if present
    })
    
    # Also rename PSI columns to be acceptor-specific
    for col in acceptor_df.columns:
        if '_PSI' in col and not col.startswith('acceptor_'):
            acceptor_df = acceptor_df.rename(columns={col: f'acceptor_{col}'})
    
    # Merge on intron_id (outer join to get all introns)
    # Also merge on metadata columns that should be identical (gene_name, intron_status, overlapping_genes)
    merge_on = ['intron_id']
    metadata_cols = ['gene_name', 'intron_status', 'overlapping_genes']
    
    # Track which metadata columns are in both dataframes for merging
    shared_metadata = []
    for col in metadata_cols:
        if col in donor_df.columns and col in acceptor_df.columns:
            merge_on.append(col)
            shared_metadata.append(col)
    
    integrated = pd.merge(
        donor_df,
        acceptor_df,
        on=merge_on,
        how='outer',
        suffixes=('', '_dup')
    )
    
    # Drop any remaining duplicate columns (with _dup suffix)
    dup_cols = [col for col in integrated.columns if col.endswith('_dup')]
    if dup_cols:
        integrated = integrated.drop(columns=dup_cols)
    
    # Also check for any unmerged metadata columns that might appear twice
    # (e.g., if one was missing from donor or acceptor but both have it after merge)
    for col in metadata_cols:
        if col not in shared_metadata:
            # This column wasn't in both dataframes, might appear separately
            # Look for the column in integrated result
            if col in integrated.columns:
                # Check if there's a duplicate with different name or if it appears multiple times
                col_count = integrated.columns.tolist().count(col)
                if col_count > 1:
                    logger.warning(f"Found {col_count} copies of {col} - keeping only first occurrence")
                    # Keep only first occurrence
                    cols_to_keep = []
                    seen_col = False
                    for c in integrated.columns:
                        if c == col:
                            if not seen_col:
                                cols_to_keep.append(c)
                                seen_col = True
                        else:
                            cols_to_keep.append(c)
                    integrated = integrated[cols_to_keep]
    
    # Create summary columns
    integrated['tested_in'] = np.where(
        integrated['donor_PValue'].notna() & integrated['acceptor_PValue'].notna(),
        'both',
        np.where(integrated['donor_PValue'].notna(), 'donor_only', 'acceptor_only')
    )
    
    # Significance status
    donor_sig = integrated['donor_significant'].fillna(False)
    acceptor_sig = integrated['acceptor_significant'].fillna(False)
    donor_sig = donor_sig.astype(bool)
    acceptor_sig = acceptor_sig.astype(bool)
    
    integrated['significant_in'] = np.where(
        donor_sig & acceptor_sig,
        'both',
        np.where(
            donor_sig,
            'donor_only',
            np.where(acceptor_sig, 'acceptor_only', 'neither')
        )
    )
    
    # For introns significant in both, check direction consistency
    both_sig = (donor_sig & acceptor_sig)
    if both_sig.any():
        same_direction = np.sign(integrated.loc[both_sig, 'donor_logFC']) == \
                        np.sign(integrated.loc[both_sig, 'acceptor_logFC'])
        integrated.loc[both_sig, 'direction_consistent'] = same_direction
    else:
        integrated['direction_consistent'] = np.nan
    
    # Best (most significant) analysis
    integrated['best_FDR'] = integrated[['donor_FDR', 'acceptor_FDR']].min(axis=1)
    integrated['best_analysis'] = np.where(
        integrated['donor_FDR'] <= integrated['acceptor_FDR'],
        'donor',
        'acceptor'
    )
    
    # For best analysis, get the corresponding logFC
    integrated['best_logFC'] = np.where(
        integrated['best_analysis'] == 'donor',
        integrated['donor_logFC'],
        integrated['acceptor_logFC']
    )
    
    # For best analysis, get the corresponding delta_PSI if available
    if 'donor_delta_PSI' in integrated.columns or 'acceptor_delta_PSI' in integrated.columns:
        integrated['best_delta_PSI'] = np.where(
            integrated['best_analysis'] == 'donor',
            integrated.get('donor_delta_PSI', np.nan),
            integrated.get('acceptor_delta_PSI', np.nan)
        )
    
    # Collect PSI columns from both analyses (excluding ones already in priority list)
    already_in_priority = {'best_delta_PSI', 'donor_delta_PSI', 'acceptor_delta_PSI'}
    all_psi_cols = [col for col in integrated.columns 
                    if ('_PSI' in col or 'delta_PSI' in col) and col not in already_in_priority]
    
    acceptor_psi_cols = [col for col in all_psi_cols if col.startswith('acceptor_')]
    donor_psi_cols = [col for col in all_psi_cols if col.startswith('donor_')]
    
    # Reorder columns: metadata → ALL acceptor (including PSI) → ALL donor (including PSI) → best summary
    priority_cols = [
        'intron_id',
        'gene_name',
        'intron_status',
        'overlapping_genes',
        # ALL acceptor columns grouped together (main analysis + PSI stats)
        'acceptor_contrast',
        'acceptor_cluster',
        'acceptor_logFC',
        'acceptor_PValue',
        'acceptor_FDR',
        'acceptor_FDR_original',
        'acceptor_significant',
        'acceptor_delta_PSI',
        'acceptor_logCPM',
    ]
    
    # Add acceptor PSI summary statistics right after acceptor main columns
    priority_cols.extend(acceptor_psi_cols)
    
    # Now all donor columns (main analysis + PSI stats)
    priority_cols.extend([
        'donor_contrast',
        'donor_cluster',
        'donor_logFC',
        'donor_PValue',
        'donor_FDR',
        'donor_FDR_original',
        'donor_significant',
        'donor_delta_PSI',
        'donor_logCPM',
    ])
    
    # Add donor PSI summary statistics right after donor main columns
    priority_cols.extend(donor_psi_cols)
    
    # Best value summary columns go at the very end
    priority_cols.extend([
        'tested_in',
        'significant_in',
        'best_analysis',
        'best_FDR',
        'best_logFC',
        'best_delta_PSI',
        'direction_consistent',
    ])
    
    # Keep priority columns first (only if they exist), then any remaining columns
    # Exclude F statistic columns from final output
    existing_priority_cols = [col for col in priority_cols if col in integrated.columns]
    other_cols = [col for col in integrated.columns 
                  if col not in priority_cols and col not in ['donor_F', 'acceptor_F']]
    integrated = integrated[existing_priority_cols + other_cols]
    
    # Sort by best FDR
    integrated = integrated.sort_values('best_FDR')
    
    return integrated
    
    # Direction consistency for those significant in both
    both_sig_count = (integrated['significant_in'] == 'both').sum()
    if both_sig_count > 0:
        both_sig_mask = integrated['significant_in'] == 'both'
        consistent_count = (integrated.loc[both_sig_mask, 'direction_consistent'] == True).sum()
        opposite_count = (integrated.loc[both_sig_mask, 'direction_consistent'] == False).sum()
        logger.info(f"\nOf {both_sig_count} introns significant in both:")
        logger.info(f"  Consistent direction: {consistent_count}")
        logger.info(f"  Opposite direction: {opposite_count}")
    
    return integrated


def create_summary_table(integrated_df):
    """
    Create a high-level summary table of integration results.
    
    Args:
        integrated_df: Integrated results DataFrame
        
    Returns:
        Summary DataFrame
    """
    summary = {
        'total_introns': len(integrated_df),
        'tested_both': (integrated_df['tested_in'] == 'both').sum(),
        'tested_donor_only': (integrated_df['tested_in'] == 'donor_only').sum(),
        'tested_acceptor_only': (integrated_df['tested_in'] == 'acceptor_only').sum(),
        'significant_both': (integrated_df['significant_in'] == 'both').sum(),
        'significant_donor_only': (integrated_df['significant_in'] == 'donor_only').sum(),
        'significant_acceptor_only': (integrated_df['significant_in'] == 'acceptor_only').sum(),
        'significant_neither': (integrated_df['significant_in'] == 'neither').sum(),
    }
    
    # Direction consistency for those significant in both
    both_sig = integrated_df[integrated_df['significant_in'] == 'both']
    if len(both_sig) > 0:
        # Count True/False explicitly to avoid issues with NaN and boolean negation
        consistent = (both_sig['direction_consistent'] == True).sum()
        opposite = (both_sig['direction_consistent'] == False).sum()
        summary['both_sig_consistent_direction'] = int(consistent)
        summary['both_sig_opposite_direction'] = int(opposite)
    else:
        summary['both_sig_consistent_direction'] = 0
        summary['both_sig_opposite_direction'] = 0
    
    return pd.DataFrame([summary])


def main():
    parser = argparse.ArgumentParser(
        description="Integrate donor and acceptor differential splicing results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    parser.add_argument(
        "--donor_results",
        type=str,
        required=True,
        help="Donor intron results file (.intron_results.tsv)",
    )
    
    parser.add_argument(
        "--acceptor_results",
        type=str,
        required=True,
        help="Acceptor intron results file (.intron_results.tsv)",
    )
    
    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="Output prefix for integrated results",
    )
    
    parser.add_argument(
        "--gtf",
        type=str,
        default=None,
        help="GTF annotation file for gene annotation and known/novel intron status",
    )
    
    parser.add_argument(
        "--fdr_threshold",
        type=float,
        default=0.05,
        help="FDR threshold for filtering significant results",
    )
    
    args = parser.parse_args()
    
    # Load results
    donor_df, acceptor_df = load_results(args.donor_results, args.acceptor_results)
    
    # Integrate
    integrated = integrate_intron_results(donor_df, acceptor_df)
    
    # Annotate with GTF if provided
    if args.gtf:
        integrated = annotate_introns_with_gtf(integrated, args.gtf)
    else:
        logger.info("No GTF file provided, skipping gene annotation")
    
    # Write all integrated results
    all_results_file = f"{args.output_prefix}.integrated_results.tsv"
    logger.info(f"\nWriting all integrated results to {all_results_file}")
    integrated.to_csv(all_results_file, sep="\t", index=False, na_rep='NA')
    
    # Filter to significant introns (in at least one analysis)
    significant = integrated[integrated['significant_in'] != 'neither']
    
    if len(significant) > 0:
        sig_results_file = f"{args.output_prefix}.significant_integrated.tsv"
        logger.info(f"Writing {len(significant)} significant introns to {sig_results_file}")
        significant.to_csv(sig_results_file, sep="\t", index=False, na_rep='NA')
        
        # Create summary table
        summary = create_summary_table(integrated)
        summary_file = f"{args.output_prefix}.integration_summary.tsv"
        logger.info(f"Writing integration summary to {summary_file}")
        summary.to_csv(summary_file, sep="\t", index=False, na_rep='NA')
    else:
        logger.warning("No significant introns found in either analysis")
    
    logger.info("\nIntegration complete!")


if __name__ == "__main__":
    main()
