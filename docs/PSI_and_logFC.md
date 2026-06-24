# PSI and `logFC`: how they relate, and can you derive one from the other?

Diff-Splice-Finder reports two effect-size measures for each intron:

- **`logFC`** — the edgeR group coefficient (log2 scale), estimated with the
  shared cluster-total offset.
- **`delta_PSI`** — the difference in mean intron usage proportion between groups.

They describe the **same** usage shift on **two different scales** (a log-ratio
vs an arithmetic difference). This note makes the relationship precise and
answers the common question: *can I get one from the other?*

Short answer: **only if you also know the baseline PSI.** Neither quantity can be
recovered from the other in isolation.

---

## 1. Definitions

For an intron in group `g`:

```
PSI_g = intron_count_g / shared_cluster_total_g          # usage proportion, in [0, 1]
```

(`shared_cluster_total = max(donor_total, acceptor_total)`, the same denominator
edgeR uses as its offset — see [SHARED_OFFSETS.md](SHARED_OFFSETS.md).)

```
delta_PSI = PSI_A - PSI_B                                 # difference, in [-1, 1]
logFC     = log2( PSI_A / PSI_B )   (model-estimated)      # log-ratio, in (-inf, inf)
```

where `A` is the perturbation/test group and `B` the reference (control) group.

---

## 2. Why `logFC` is a PSI log-ratio (the exact identity)

Because PSI is `count / cluster_total`, the **log-ratio of PSI factorizes exactly**:

```
log2(PSI_A / PSI_B)
   = log2(count_A / count_B)  -  log2(clusterTotal_A / clusterTotal_B)
   = (raw intron-count logFC)  -  (cluster-total logFC)
```

The second term is exactly the **offset** edgeR puts in the GLM:

```
log(expected count) = group_effect + log(shared cluster total)
                                      └──────── offset ────────┘
```

Subtracting the offset is the same as dividing the count by the cluster total, so
the fitted **group effect is the log-fold-change of PSI**, not of the raw count.
That is why the reported `logFC` is an *offset-adjusted intron-usage* log
fold-change rather than an expression fold-change.

> **`logFC` is the model-based estimate of `log2(PSI_A / PSI_B)`.** It is not
> numerically identical to the log-ratio of the reported `mean_PSI` columns,
> because edgeR uses GLM-fitted, dispersion-weighted, replicate-aware values
> (a ratio of fitted means), whereas `mean_PSI` is the mean of per-sample
> ratios. The two diverge most when counts are low or a denominator PSI is near
> zero (mean-of-ratios ≠ ratio-of-means). Use `logFC` for the model estimate and
> `delta_PSI` for the observed magnitude.

---

## 3. Converting between `logFC` and `delta_PSI`

Treating `logFC` as `log2(PSI_A / PSI_B)`, and writing the baseline (control)
proportion as `PSI_B`:

```
PSI_A     = PSI_B · 2^logFC
delta_PSI = PSI_A - PSI_B = PSI_B · (2^logFC - 1)
logFC     = log2( 1 + delta_PSI / PSI_B )
```

So the two are interconvertible **only with a third quantity, the baseline
`PSI_B`** (`control_mean_PSI` in the results). The pipeline reports it, so you
can convert row by row — but you cannot go from one to the other on its own.

### Worked example: same `logFC`, very different `delta_PSI`

`logFC = +1` means usage doubles (2×). The absolute shift depends entirely on
where it started:

| `control_mean_PSI` (PSI_B) | `logFC` | implied `PSI_A` | `delta_PSI` |
|---|---|---|---|
| 0.02 | +1.0 | 0.04 | **+0.02** |
| 0.40 | +1.0 | 0.80 | **+0.40** |

Same fold-change, a 20× difference in absolute usage shift. The reverse holds
too — a fixed `delta_PSI = +0.10` is `logFC = log2(1 + 0.10/0.05) = +1.58` from a
5% baseline but only `log2(1 + 0.10/0.80) = +0.17` from an 80% baseline.

---

## 4. Why both are reported (practical implications)

- **`logFC` is a ratio** — unbounded and **hypersensitive at low PSI**. A rare
  intron going from 0.5% to 1% usage is `logFC = +1` (a "2-fold" change) but a
  biologically trivial `delta_PSI = +0.005`.
- **`delta_PSI` is a difference** — bounded to `[-1, 1]` and dominated by introns
  that carry a meaningful share of the cluster.

Because of this, significance requires **both**:

- `FDR < 0.05` (statistical reliability, from the count-based QL F-test), **and**
- `|delta_PSI| ≥ threshold` (a real change in usage proportion;
  `--min_delta_psi`, default 0.05).

When `--min_delta_psi` is enabled, the pipeline first filters to introns with
`|delta_PSI| >= threshold`, then recalculates Benjamini-Hochberg FDR from the
remaining `PValue` values. In the delta-PSI-filtered output, `FDR_original`
stores the FDR from the full tested set before this filtering, and `FDR` stores
the recomputed value used for final significance calls. `logFC` is not used to
choose the rows for this FDR recalculation.

The final `edgeR_results.significant_introns.tsv` file is therefore not just a
copy of every low-FDR row. With delta-PSI filtering enabled, reported rows must
come from the delta-PSI-filtered set, pass the configured FDR threshold using
the recomputed `FDR`, and pass the configured `--min_logFC` threshold.

The FDR/`logFC` side answers *whether* usage changed; `delta_PSI` answers *how
much* of the local splicing output actually moved. They are complementary, not
redundant.

---

## 5. Summary

| | `logFC` | `delta_PSI` |
|---|---|---|
| scale | log2 ratio of PSI | arithmetic difference of PSI |
| range | `(-inf, inf)` | `[-1, 1]` |
| source | edgeR GLM coefficient (offset-adjusted) | mean PSI difference |
| sensitive to | small switches at **low** PSI | switches among **abundant** introns |
| convert to the other | needs baseline `PSI_B` | needs baseline `PSI_B` |

`logFC` and `delta_PSI` are two scales of the same usage shift, linked by
`delta_PSI = PSI_B·(2^logFC − 1)`. You can derive either from the other **only**
with the control-group baseline `PSI_B`; from one alone, you cannot.
