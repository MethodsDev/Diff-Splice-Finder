# Understanding Offsets in GLMs
## Why β₁ Means Different Things With and Without an Offset

This document builds intuition around **what changes — and what does not — when an offset is included in a generalized linear model (GLM)** for RNA-seq analysis.

The focus is on understanding:
- why the base GLM is a straight line with a slope,
- why β₁ is *always* a log fold-change mathematically,
- and how an offset fundamentally changes the **biological interpretation** of that fold-change.

This tutorial complements the main Diff-Splice-Finder GLM tutorial.

---

## 1. A two-point line

With a two-group comparison, the GLM

```
log(μ[i, s]) = β0[i] + β1[i] × Group[s]
```

is just a **line with two x-values**:

- `Group = 0` → reference group
- `Group = 1` → comparison group

The expected means are:
- μ₀ for Group 0
- μ₁ for Group 1

The model fits the straight line connecting:

```
(0, log(μ₀))   and   (1, log(μ₁))
```

The **slope of this line is β₁**.

---

## 2. Why this is a linear model

A linear model has the form:

```
y = a + b × x
```

Mapping to RNA-seq:

| Linear model | RNA-seq GLM |
|-------------|-------------|
| y | log(μ) |
| a | β0 |
| b | β1 |
| x | Group |

So β₁ is literally the change in log(expected count) per unit change in Group.

---

## 3. β₁ is always a log fold-change

From the definition of the slope:

```
β1 = log(μ₁) − log(μ₀)
```

which implies:

```
β1 = log( μ₁ / μ₀ )
```

This is true **with or without an offset**.

What changes is **what μ represents**.

---

## 4. No offset: abundance model

Model:

```
log(μ) = β0 + β1 × Group
```

Here:

```
μ = expected raw count
```

So:

```
β1 = log( raw counts in Group 1 / raw counts in Group 0 )
```

### Interpretation

- β₁ measures **abundance change**
- It mixes expression and sequencing depth
- This is appropriate for **gene-level differential expression**

---

## 5. Why this fails for splicing

If gene expression doubles but splicing does not change:

- all introns double in count
- β₁ ≈ log(2) for every intron

This produces false splicing signals.

---

## 6. With an offset: conditional model

Model:

```
log(μ) = β0 + β1 × Group + log(T)
```

Exponentiating:

```
μ = T × exp(β0 + β1 × Group)
```

Here:

```
μ = (available evidence) × (usage proportion)
```

---

## 7. What the offset does

The offset tells the model:

> “These counts are conditional on a shared total.”

It does **not** change the slope mathematically.
It changes what the slope refers to biologically.

---

## 8. Why β₁ becomes a usage fold-change

With the offset:

```
μ₀ = T × p₀
μ₁ = T × p₁
```

So:

```
μ₁ / μ₀ = p₁ / p₀
```

And therefore:

```
β1 = log( p₁ / p₀ )
```

The cluster total cancels exactly.

---

## 9. Same math, different question

| Model | Question |
|-----|----------|
| No offset | “Is there more RNA?” |
| With offset | “Is usage redistributed?” |

---

## 10. One-sentence takeaway

> Offsets do not change the mathematics of the GLM; they change the biological meaning of μ, and therefore the interpretation of β₁.
