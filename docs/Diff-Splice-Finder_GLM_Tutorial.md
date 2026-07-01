# Diff-Splice-Finder
## GLM Tutorial: Intron Usage Analysis with edgeR Offsets

This document explains how Diff-Splice-Finder uses edgeR GLMs to test
differential intron usage. The key idea is that raw junction counts are modeled
with an offset so the group coefficient reflects usage relative to a selected
local denominator, not raw abundance alone.

## 1. What Problem Are We Solving?

We want to detect splicing changes, not expression changes.

Raw intron counts alone cannot answer that because gene expression, local
coverage, and sequencing depth can change all junction counts at a locus. DSF
therefore models an intron's count relative to a denominator representing local
available evidence.

## 2. What Is Modeled?

Each observation is:

- intron `i`
- sample `s`
- observed numerator count `y[i, s]`
- selected denominator `D[i, s]`

The default DSF denominator is:

```text
D[i, s] = max(donor_site_window_depth, acceptor_site_window_depth)
```

Strict modes can instead use:

```text
D[i, s] = strict_local_depth
```

where `strict_local_depth` is the maximum of donor and acceptor
splice-decision depths. In `gene_median_strict_depth` mode, PSI still uses
`strict_local_depth`, but edgeR uses the per-gene median strict depth as `D`.

## 3. The Naive Model

Without an offset:

```text
log(mu[i, s]) = beta0[i] + beta1[i] * Group[s]
```

Here `mu` is the expected raw count. If a locus doubles in expression while
splicing proportions stay unchanged, this model can produce a group effect even
though splicing did not change.

## 4. The Offset Model

DSF supplies the denominator as a fixed offset:

```text
log(mu[i, s]) = beta0[i] + beta1[i] * Group[s] + log(D[i, s])
```

Equivalently:

```text
mu[i, s] = D[i, s] * exp(beta0[i] + beta1[i] * Group[s])
```

The coefficient is therefore fitted on the count scale while accounting for the
available local evidence.

## 5. Interpreting beta1

For a two-group comparison:

```text
Group = 0  reference
Group = 1  comparison
```

The fitted means are:

```text
mu0 = D * exp(beta0)
mu1 = D * exp(beta0 + beta1)
```

If the same denominator definition is used for PSI and edgeR, the denominator
cancels in the group ratio:

```text
mu1 / mu0 = [D * exp(beta0 + beta1)] / [D * exp(beta0)]
          = exp(beta1)
```

So `beta1` is the log fold-change in numerator usage relative to the selected
denominator. edgeR reports this on the log2 scale as `logFC`.

In `gene_median_strict_depth` mode, this interpretation is slightly different:
`logFC` is still the model coefficient under the supplied gene-median exposure,
but it is not the exact log-ratio of the reported strict-local mean PSI values.

## 6. Why Library-Size Normalization Is Disabled

edgeR normally estimates library-size normalization factors. DSF sets
`norm.factors = 1` because normalization is supplied through the denominator
offset. Applying both would mix global library normalization into a local usage
model.

## 7. What delta_PSI Adds

`logFC` is a model-based relative effect. `delta_PSI` is an observed absolute
difference:

```text
delta_PSI = mean_PSI_group1 - mean_PSI_group2
PSI = numerator_count / selected_denominator
```

Both are useful. `logFC` is sensitive to relative changes, especially at low
baseline PSI. `delta_PSI` reports how much of the local splicing output moved.

## 8. Filtering and FDR

DSF filters before edgeR:

- canonical splice motif unless `--keep_noncanonical` is used
- minimum numerator count
- minimum samples with nonzero numerator count
- minimum selected-denominator depth
- optional `--min_delta_psi`

edgeR FDR is computed over the introns that pass these prefilters. The current
pipeline does not recompute FDR after edgeR.

## 9. Strandedness

In default DSF mode, `--site_depth_strand_mode` filters the site-depth
denominator by transcript strand. Junction discovery itself is not explicitly
orientation-filtered; canonical junction strand is inferred from splice motifs.

Strict-depth counting currently falls back to genomic/unstranded behavior if a
stranded mode is supplied.

## 10. Summary

- DSF tests each intron once.
- Counts remain counts; PSI is not modeled directly.
- The selected denominator enters edgeR as a fixed log offset.
- Default denominator: max donor/acceptor site-depth window.
- Strict denominator: max donor/acceptor splice-decision depth.
- `logFC` is an offset-adjusted model coefficient.
- `delta_PSI` is an observed absolute usage difference.
