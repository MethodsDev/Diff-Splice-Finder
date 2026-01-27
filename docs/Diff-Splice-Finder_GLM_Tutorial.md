# Diff-Splice-Finder  
## GLM Tutorial: Intron Usage Analysis with edgeR Offsets

This document explains **how generalized linear models (GLMs)** are used in Diff-Splice-Finder to detect **differential intron usage** in bulk RNA-seq data, for both short-read (Illumina) and long-read (PacBio Kinnex) technologies.

The goal is conceptual understanding rather than mathematical rigor.

---

## 1. What problem are we solving?

We want to detect **splicing changes**, not expression changes.

Specifically, we want to answer:

> *Within a local splicing choice (a donor or acceptor site), did the **relative usage** of an intron change between conditions?*

Raw intron counts alone cannot answer this, because:
- gene expression may increase or decrease, affecting all introns
- sequencing depth varies across samples

Splicing is **compositional**: introns compete for usage within a cluster.

---

## 2. What is being modeled?

### Observations

Each observation corresponds to:

- **intron** `i` (a splice junction)
- **sample** `s` (a bulk RNA-seq replicate)

We observe:

- `y[i, s]` = number of reads supporting intron *i* in sample *s*

Introns are grouped into **clusters**:
- donor clusters (shared 5′ splice site)
- acceptor clusters (shared 3′ splice site)

For each cluster `C` and sample `s`, we also compute:

- `T[C, s]` = total intron-support reads in that cluster and sample  
  (the sum of counts across introns in the cluster)

---

## 3. Why splicing is compositional (the pizza analogy)

Think of a splice site as a pizza:

- the **pizza size** = total splicing evidence at that site (`T[C, s]`)
- each **slice** = one intron’s usage

If one intron gets more slices, others must get fewer.

---

## 4. What is a GLM (in plain English)?

A **generalized linear model (GLM)** predicts an expected mean count and accounts for variability.

edgeR uses:
- a **negative binomial** distribution
- a **log link function**

---

## 5. Why β₁ is a log fold-change in intron *usage* (explicit derivation)

Model:

```
log(μ[i, s]) = β0[i] + β1[i] × Group[s] + log(T[C, s])
```

**Group = 0 (reference):**

```
μ₀ = T × exp(β0)
```

**Group = 1 (comparison):**

```
μ₁ = T × exp(β0 + β1)
```

Taking the ratio:

```
μ₁ / μ₀ = exp(β1)
```

So:

```
β1 = log( usage proportion in Group 1 / usage proportion in Group 0 )
```

This holds because the cluster total `T` cancels due to the offset.

---

## 6. Key takeaways

- Splicing is compositional
- Offsets convert count models into proportion tests
- `logFC` represents change in intron usage, not expression
