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

Each intron gets its **own GLM**, with its own coefficients.

---

## 6. The naive (wrong) model

```
log(μ[i, s]) = β0[i] + β1[i] × Group[s]
```

This treats an intron like a gene and tests raw intron counts.

**Problem:**  
If gene expression changes, all introns shift together, producing false splicing signals.

---

## 7. The key idea: offsets convert counts into proportions

To model **intron usage**, we include a fixed **offset**:

```
log(μ[i, s]) = β0[i] + β1[i] × Group[s] + log(T[C, s])
```

Where:
- `T[C, s]` is the total intron-support evidence in intron *i*’s cluster.

After exponentiating:

```
μ[i, s] = T[C, s] × exp(β0[i] + β1[i] × Group[s])
```

Interpretation:

- `T[C, s]` = total splicing evidence (pizza size)
- `exp(β0 + β1 × Group)` = intron’s **usage share**

The model now tests **relative usage**, not absolute abundance.

---

## 8. Conceptually: Group = 0 or 1

In the simplest two-group comparison:

- **Group = 0** → reference group (e.g. `control`)
- **Group = 1** → comparison group (e.g. `TDP43`)

The model:

```
log(μ[i, s]) = β0[i] + β1[i] × Group[s] + log(T[C, s])
```

has the following interpretations:

### When Group = 0 (reference)

```
log(μ) = β0[i] + log(T)
```

Expected count:
```
μ = T × exp(β0[i])
```

This defines the **baseline usage proportion** of intron *i*.

### When Group = 1 (comparison)

```
log(μ) = β0[i] + β1[i] + log(T)
```

Expected count:
```
μ = T × exp(β0[i] + β1[i])
```

### What β1 means

Subtracting the two cases:

```
β1[i] = log( usage in Group 1 / usage in Group 0 )
```

So **β1 is exactly the log fold-change in intron usage proportion**.

This is why `logFC` from edgeR should be interpreted as:

> “How much more (or less) of the local splicing pool this intron receives in one group vs the other.”

---

## 9. What does this look like in R?

In practice, you supply a factor, and R encodes it as 0/1:

```r
group <- factor(
  c("control","control","control",
    "TDP43","TDP43","TDP43"),
  levels = c("control","TDP43")
)

design <- model.matrix(~ group)
```

This produces a column `groupTDP43` with values:

- `0` for control samples
- `1` for TDP43 samples

The coefficient for `groupTDP43` is **β1**.

---

## 10. Interpreting logFC

Because the model is on the log scale:

- `logFC = 1` → intron usage is **2× higher**
- `logFC = -1` → intron usage is **2× lower**

Important:
- This refers to **usage proportion**, not expression
- A gene can go up in expression while an intron goes *down* in usage

---

## 11. Why library-size normalization is disabled

edgeR normally normalizes by library size.

In Diff-Splice-Finder:
- normalization is done via **cluster-total offsets**
- applying library-size normalization would distort usage estimates

Therefore:
- `norm.factors = 1`
- offsets are the *only* normalization

---

## 12. Why negative binomial?

RNA-seq counts show more variability than Poisson models allow.

The negative binomial includes a **dispersion** parameter that captures:
- biological variability
- technical variability

edgeR estimates dispersion from replicates and uses it for inference.

---

## 13. Why this works for long reads

In long-read RNA-seq:
- a single read may support multiple introns
- each intron receives +1 evidence

Because:
- the denominator (`T[C, s]`) is built from the same evidence
- the test is relative within clusters

The offset-based GLM still correctly tests **relative intron usage**.

---

## 14. Mental model summary

For each intron *i* and sample *s*:

1. The cluster has `T[C, s]` total splicing evidence
2. The intron has a usage proportion `p[i, s]`
3. Observed counts are noisy realizations of  
   `T[C, s] × p[i, s]`
4. The GLM tests whether `p[i, s]` changes with condition

---

## 15. Key takeaways

- Splicing is compositional
- Raw intron counts are misleading
- Offsets convert count models into proportion tests
- edgeR GLMs isolate splicing from expression
- `logFC` represents **change in intron usage proportion**, not abundance

For implementation details, see:
- `README.md`
- `AI_ONBOARDING.md`
