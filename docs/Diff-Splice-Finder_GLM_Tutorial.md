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

### Switching without expression change

| Condition | Intron A | Intron B | Intron C | Cluster Total |
|---------|----------|----------|----------|----------------|
| Group 1 | 80 | 15 | 5 | 100 |
| Group 2 | 40 | 10 | 50 | 100 |

The pizza stayed the same size, but the slices changed.
This is a **splicing change**.

### Expression change without switching

| Condition | Intron A | Intron B | Intron C | Cluster Total |
|---------|----------|----------|----------|----------------|
| Group 1 | 80 | 15 | 5 | 100 |
| Group 2 | 160 | 30 | 10 | 200 |

The pizza doubled, but slice proportions stayed the same.
This is **not** a splicing change.

---

## 4. What is a GLM (in plain English)?

A **generalized linear model (GLM)** predicts an expected mean count and accounts for variability.

edgeR uses:
- a **negative binomial** distribution (handles biological variability)
- a **log link function**

This means it models:

```
log(expected count) = linear combination of predictors
```

---

## 5. What does the index `i` mean?

When you see:

```
μ[i, s]
```

read it as:

> the expected count for **intron i** in **sample s**

- `i` indexes **introns**
- `s` indexes **samples**

Each intron gets its **own GLM**.

---

## 6. The naive (wrong) model

```
log(μ[i, s]) = β0[i] + β1[i] × Group[s]
```

This tests raw intron counts and is confounded by expression changes.

---

## 7. The key idea: offsets convert counts into proportions

```
log(μ[i, s]) = β0[i] + β1[i] × Group[s] + log(T[C, s])
```

Where `T[C, s]` is the cluster total.

After exponentiating:

```
μ[i, s] = T[C, s] × exp(β0[i] + β1[i] × Group[s])
```

The group coefficient becomes a log fold-change in **intron usage proportion**.

---

## 8. What does “Group = 0 or 1” mean?

- reference group → 0
- comparison group → 1

In R, factors are automatically encoded this way via `model.matrix()`.

---

## 9. Interpreting logFC

- `logFC = 1` → 2× higher intron usage
- `logFC = -1` → 2× lower intron usage

This refers to **usage proportion**, not expression.

---

## 10. Why library-size normalization is disabled

Normalization is done via cluster-total offsets.
Library-size normalization would distort usage estimates.

---

## 11. Why negative binomial?

It models biological variability better than Poisson models.

---

## 12. Why this works for long reads

Long reads support multiple introns per read.
Because both numerator and denominator are intron-support evidence,
relative usage remains meaningful.

---

## 13. Mental model summary

Observed counts ≈ (cluster total) × (intron usage share) + noise.

The GLM tests whether usage share changes between conditions.

---

## 14. Key takeaways

- Splicing is compositional
- Offsets isolate splicing from expression
- edgeR GLMs test intron usage proportions
- logFC reflects splicing changes, not abundance
