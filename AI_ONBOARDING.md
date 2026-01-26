# Intron-Focused Differential Splicing Analysis

## AI Onboarding / Design Context

### Purpose

We are building a **bulk RNA-seq splicing analysis pipeline** that works **consistently for both short-read Illumina RNA-seq and long-read PacBio Kinnex RNA-seq**, using **intron-level evidence** and **edgeR** for differential testing.

The goal is **not isoform switching per se**, but robust detection of **differential intron usage** (splicing changes) while avoiding confounding by gene expression.

This document explains:

* what the features are,
* how counts are constructed,
* how splicing is modeled statistically,
* and what assumptions are intentional.

---

## Core Design Principles

1. **Splicing is compositional**

   * Splicing reflects *choices among introns* within a local locus.
   * Absolute read depth is not the signal of interest.
   * We test **relative intron usage within clusters**, not global expression.

2. **Counts remain counts**

   * We use **integer intron-support counts**.
   * No TPM/PSI transformation before modeling.
   * All normalization happens via **offsets in the GLM**.

3. **Same abstraction for short + long reads**

   * Both technologies are reduced to the same representation:

     ```
     intron × sample → integer counts
     ```
   * Long reads may contribute to multiple introns per read; this is intentional.

4. **Annotation-light**

   * Introns are defined from observed splice junctions.
   * Strand is inferred from splice-site dinucleotides.
   * No dependence on transcript models.

---

## Input Data Assumptions

### Samples

* Bulk RNA-seq
* Multiple biological replicates per sample type
* No UMIs, no deduplication
* Optional batch covariates

### Technologies

* Short reads: Illumina RNA-seq
* Long reads: PacBio **bulk Kinnex**

---

## Feature Definition

### Intron (junction) features

Each feature corresponds to a spliced intron defined as:

```
(chr, donor_position, acceptor_position, strand)
```

* Strand is inferred from splice-site motifs (GT–AG, GC–AG, AT–AC).
* Counts represent **number of reads supporting that intron**.

### Counting rules

* **Short reads**: each junction-spanning read increments one intron.
* **Long reads**: a read spanning *k* introns increments **each of the k introns by +1**.

  * This intentionally measures *intron support*, not molecules.
  * Long isoforms therefore contribute more intron evidence.

---

## Splicing Clusters (Denominators)

To model intron *usage*, introns are grouped into **local splicing clusters**, which define the denominator for proportions.

### Default cluster types (recommended)

We typically run analyses using **both**:

1. **Donor clusters**

   * All introns sharing the same 5′ splice site
2. **Acceptor clusters**

   * All introns sharing the same 3′ splice site

These clusters:

* Stay small (important for long-read sparsity)
* Capture classic alt-5′ / alt-3′ and many cassette-like events
* Are stable across technologies

(Connected-component clustering à la LeafCutter is optional but not default.)

---

## Statistical Model (edgeR)

### What we test

For each intron *i* in cluster *C*:

> Does the **fraction of intron usage within its cluster** differ between sample types?

Not:

* gene expression,
* total junction abundance,
* isoform expression.

### Model structure

We fit a **negative binomial GLM** using edgeR:

[
\log(\mu_{i,s}) =
X_s \beta_i + \log(T_{C,s})
]

Where:

* ( \mu_{i,s} ) = expected intron count
* ( X_s ) = design matrix (e.g. sample type, batch)
* ( T_{C,s} ) = total intron-support reads in cluster *C* for sample *s*
* ( \log(T_{C,s}) ) is supplied as an **offset**

### Key consequence

The coefficient for sample type represents a **log fold-change in intron usage**, not expression.

### edgeR specifics

* We **do not apply library-size normalization** (`norm.factors = 1`)
* All normalization is via the cluster-total offset
* Use QL GLMs (`glmQLFit`, `glmQLFTest`) with robust dispersion

---

## Filtering Strategy

Because junction data (especially long reads) are sparse:

### Cluster-level filters

* Require cluster total ≥ *N* reads in ≥ *K* samples
  (e.g. ≥20 reads in ≥3 samples)

### Intron-level filters

* Require:

  * total count ≥10 across all samples
  * nonzero counts in ≥2–3 samples

These thresholds are tunable and dataset-dependent.

---

## Outputs

### Intron-level results

For each intron:

* logFC = log change in **usage**
* p-value, FDR
* cluster membership

### Cluster-level results (recommended)

* Aggregate intron p-values per cluster (e.g. ACAT)
* Report:

  * cluster-level FDR
  * contributing introns and directionality

This yields interpretable “splicing events” similar in spirit to LeafCutter, but using edgeR.

---

## Long-Read–Specific Notes

* Counting multiple introns per read is intentional.
* This slightly upweights long isoforms; acceptable for intron-support analyses.
* No isoform reconstruction or switching analysis is performed here.
* Future extensions may add intron-chain or isoform-fraction analyses, but those are out of scope for now.

---

## Non-Goals (Explicit)

* No transcript assembly
* No isoform EM / abundance estimation
* No PSI from exon inclusion
* No single-cell modeling
* No UMI deduplication

---

## Mental Model Summary

> “We treat each intron as a feature, but we test it **relative to its local splicing alternatives**, using edgeR with offsets so that **expression changes cancel out and splicing choices remain**.”

This is conceptually similar to LeafCutter, but:

* simpler,
* more flexible,
* and directly applicable to long reads.

---

If you want, next I can:

* write this as a **system prompt** instead of markdown,
* generate a **reference implementation skeleton** (`counts → clusters → offsets → edgeR → results`),
* or create a **checklist / unit-test spec** for validating the pipeline.

