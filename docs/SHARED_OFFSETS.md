# Shared Offsets: Understanding the Denominator Strategy

## Overview

Diff-Splice-Finder uses a **shared offset approach** where each intron's normalization factor is computed as the **maximum** of its donor and acceptor cluster totals. This document explains why this is necessary and how it works.

## The Problem: Singleton Cluster Artifacts

### Scenario: Novel Acceptor with Common Donor

Consider this real example from `chr8:79611215-79616821`:

**Control samples:**
- Intron counts: 2, 2, 3 reads
- **Acceptor cluster**: Only this intron uses this acceptor site
  - Acceptor cluster total: 2, 2, 3
- **Donor cluster**: Many introns share this donor site
  - Donor cluster total: ~3500, ~3500, ~3600

**Disease samples:**
- Intron counts: 1354, 1508, 1505 reads
- Acceptor cluster total: 1354, 1508, 1505
- Donor cluster total: ~3500, ~3500, ~3600

### The Problem with Separate Analyses

**If we only use acceptor cluster total:**
```
Control PSI = 2/2 = 1.0 (100% usage!)
Disease PSI = 1354/1354 = 1.0 (100% usage!)
Delta PSI = 0 (no change detected)
```

This is misleading! The intron went from 2 reads to 1354 reads - a massive biological change. But because it's the only intron using that acceptor, PSI = 1.0 in both conditions.

**If we only use donor cluster total:**
```
Control PSI = 2/3500 = 0.0006 (0.06% usage)
Disease PSI = 1354/3500 = 0.39 (39% usage)
Delta PSI = 0.39 (huge change!)
```

This correctly captures the biological reality: the intron went from being a rare alternative to a major isoform.

## The Solution: Shared Offsets

### Computing Shared Offsets

For each intron in each sample:
1. Calculate **donor cluster total**: Sum all introns sharing the same donor site
2. Calculate **acceptor cluster total**: Sum all introns sharing the same acceptor site
3. **Use the maximum**: `shared_offset = max(donor_total, acceptor_total)`

### Why Maximum?

The maximum represents the **expression potential** at either end of the intron:

- **Novel acceptor + common donor**: Uses donor total (appropriate - reflects gene expression)
- **Common acceptor + novel donor**: Uses acceptor total (appropriate - reflects gene expression)
- **Both sites moderately used**: Uses whichever is higher (conservative, appropriate)
- **Both sites novel**: Uses the maximum of two small numbers (appropriate - truly rare alternative)

### Implementation

```python
# Calculate both cluster totals
donor_totals = introns_df.groupby('donor_cluster')[sample_cols].sum()
acceptor_totals = introns_df.groupby('acceptor_cluster')[sample_cols].sum()

# Map to each intron
donor_offsets = introns_df[['donor_cluster']].join(donor_totals, on='donor_cluster')
acceptor_offsets = introns_df[['acceptor_cluster']].join(acceptor_totals, on='acceptor_cluster')

# Use maximum
shared_offsets = donor_offsets.combine(acceptor_offsets, np.maximum)
```

## Usage in Analysis

### edgeR GLM

```
log(μ_i,s) = X_s × β_i + log(shared_offset_i,s)
```

- Same intron gets identical offset regardless of context
- Normalizes for maximum expression potential at either splice site
- Prevents spurious significance from singleton clusters

### PSI Calculation

```
PSI = intron_count / shared_cluster_total
```

- Uses the **same denominators** that edgeR uses
- Ensures consistency between statistical model and biological interpretation
- PSI values reflect true proportional usage
- Eliminates singleton cluster artifact

## Biological Interpretation

### What Does the Shared Offset Represent?

The shared offset answers: **"What's the maximum splicing activity we've observed at either end of this intron?"**

This is biologically meaningful because:
- An intron's usage is constrained by **both** its donor and acceptor sites
- If either site has high activity, the intron has **potential** for high usage
- Novel sites paired with common sites should be normalized by the common site's expression

### Example Interpretations

**Case 1: Novel disease-specific exon skipping**
- Donor: chr1:1000 (used by 10 introns, total = 5000 reads)
- Acceptor: chr1:5000 (used by 2 introns, total = 500 reads)
- Your intron: chr1:1000-5000
- **Shared offset = max(5000, 500) = 5000**
- If your intron has 50 reads: PSI = 50/5000 = 1% (rare alternative)

**Case 2: Major isoform switch**
- Donor: chr2:2000 (used by 3 introns, total = 3000 reads)
- Acceptor: chr2:4000 (used by 3 introns, total = 3000 reads)  
- Your intron: chr2:2000-4000
- **Shared offset = max(3000, 3000) = 3000**
- If your intron has 1500 reads: PSI = 1500/3000 = 50% (major isoform)

**Case 3: Singleton cluster (the problem case)**
- Donor: chr3:3000 (used by 8 introns, total = 4000 reads)
- Acceptor: chr3:6000 (used by 1 intron - only yours!, total = 5 reads)
- Your intron: chr3:3000-6000
- **Shared offset = max(4000, 5) = 4000**
- If your intron has 5 reads: PSI = 5/4000 = 0.125% (truly rare)
- **Without shared offsets**: PSI = 5/5 = 100% (misleading!)

## Advantages of This Approach

### 1. Consistency
- Same intron has identical offset in all analyses
- No need to integrate separate donor/acceptor results
- Simpler interpretation

### 2. Biological Accuracy
- PSI reflects true proportional usage
- Novel alternatives correctly identified as rare
- Avoids 100% PSI for singleton clusters

### 3. Statistical Robustness
- Prevents spurious significance from low-count singleton clusters
- Reduces false positives
- Maintains sensitivity for true splicing changes

### 4. Reduced Multiple Testing
- Each intron tested once (not twice)
- Better FDR control
- Increased power

### 5. Computational Efficiency
- Single analysis instead of two
- Faster execution (~50% reduction in runtime)
- Simpler output structure

## Comparison with Other Tools

### Annotation-Based Tools (rMATS, MAJIQ, LeafCutter)
These tools avoid the singleton problem by:
- Pre-defining all alternative exons/junctions from annotation
- Always having multiple alternatives in the denominator

**Trade-off:**
- ✅ No singleton artifacts
- ❌ Miss truly novel events not in annotation
- ❌ Require high-quality annotation

### Diff-Splice-Finder Approach
Uses shared offsets to handle discovered events:
- ✅ Discovers novel events without annotation
- ✅ Prevents singleton artifacts
- ✅ Works with partial/incomplete annotation
- ✅ Appropriate normalization for any novel event

## When Does This Matter?

Shared offsets are **critical** for:
1. **Novel splice sites** (disease-specific, species without good annotation)
2. **Rare isoforms** that appear only in some conditions
3. **Long-read data** where novel events are common
4. **Non-model organisms** with incomplete annotation

Shared offsets are **less critical** but still beneficial for:
1. Well-annotated exon skipping events
2. Known alternative 5'/3' splice sites
3. Balanced alternative splicing (all alternatives ~equal usage)

## Filtering Strategy with Shared Offsets

### OR Logic for Cluster Filtering

Since we use `max(donor_total, acceptor_total)` as the offset, we only need **one** of the two clusters to be well-supported, not both. The pipeline uses **OR logic** when filtering:

```
Keep intron if: (donor_cluster passes) OR (acceptor_cluster passes)
```

**Why OR not AND?**
- Shared offset = max(donor, acceptor)
- If donor cluster is well-supported, we have a reliable offset
- If acceptor cluster is well-supported, we have a reliable offset
- We only need ONE to be reliable!

### Three Categories of Introns

**1. Both clusters pass** (most common)
- Donor cluster: ≥20 reads in ≥3 samples
- Acceptor cluster: ≥20 reads in ≥3 samples
- **Example**: Known exon skipping events
- **Offset**: Well-supported from either side

**2. Only donor cluster passes** (novel acceptors)
- Donor cluster: Well-supported
- Acceptor cluster: Low coverage (singleton or rare)
- **Example**: Novel disease-specific acceptor site
- **Offset**: Uses donor cluster total (appropriate!)

**3. Only acceptor cluster passes** (novel donors)
- Acceptor cluster: Well-supported
- Donor cluster: Low coverage (singleton or rare)
- **Example**: Novel disease-specific donor site
- **Offset**: Uses acceptor cluster total (appropriate!)

### Impact on Discovery

**With AND logic** (overly conservative):
```
Novel acceptor + common donor → REJECTED (acceptor cluster fails)
Common acceptor + novel donor → REJECTED (donor cluster fails)
```
Result: Miss many interesting novel splicing events!

**With OR logic** (appropriate):
```
Novel acceptor + common donor → KEPT (donor cluster passes)
Common acceptor + novel donor → KEPT (acceptor cluster passes)
```
Result: Discover novel events with appropriate normalization!

### Default Filtering Thresholds

```bash
--min_cluster_count 20     # Minimum reads per sample in cluster
--min_cluster_samples 3    # Samples that must meet threshold
```

An intron passes if **either**:
- Donor cluster has ≥20 reads in ≥3 samples, **OR**
- Acceptor cluster has ≥20 reads in ≥3 samples

This ensures the offset (whichever is used) is reliable.

## Summary

**Key Points:**
- Shared offsets = `max(donor_total, acceptor_total)` for each intron
- Prevents singleton cluster artifacts (PSI=1.0 for rare events)
- Same denominator used for both edgeR and PSI calculation
- Biologically accurate proportions for all events
- Each intron tested once with comprehensive information

**Bottom Line:**
Shared offsets ensure that PSI and statistical tests reflect **biological reality** rather than **technical artifacts** of cluster definitions.
