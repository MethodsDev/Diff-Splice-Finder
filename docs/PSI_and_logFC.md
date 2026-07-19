# PSI and `logFC`: how they relate, and can you derive one from the other?

Diff-Splice-Finder reports two effect-size measures for each intron:

- **`logFC`** — the model-estimated effect size (log2 scale).
- **`delta_PSI`** — the difference in mean intron usage proportion between groups.

In **offset mode** with `splice_plus_retained`, `logFC` and `delta_PSI` describe
the same usage shift on two different scales (a log-ratio vs an arithmetic
difference). With `splice_vs_rest`, `logFC` uses a gene-wide model denominator
while reported `delta_PSI` remains SPR-based, so they are complementary but not
interconvertible scales of the same proportion.
In **interaction mode**, `logFC` is a log2 focal/other odds-ratio, not a PSI
log-ratio; `delta_PSI` remains the PSI-scale effect size and is unaffected by
`--stat_mode`.
This note covers the offset-mode relationship. The common question: *can I
get one from the other?*

Short answer: **only if you also know the baseline PSI.** Neither quantity can be
recovered from the other in isolation.

---

## 1. Definitions

For an intron in group `g`:

```
PSI_g = numerator_count_g / denominator_g               # usage proportion
```

In both offset modes, the reported PSI denominator is
`max_splice_plus_retained_depth`, the maximum of left/right boundary splice
depth plus intron-side retained depth. In `splice_vs_rest` mode, only the
statistical model switches to `D^svr`, the gene-total junction count plus the
local SPR remainder. This keeps PSI, delta PSI, and the tested intron set
comparable between SPR and SVR runs.

```
delta_PSI = PSI_A - PSI_B                                 # difference, in [-1, 1]
logFC     = log2( PSI_A / PSI_B )   (model-estimated, SPR offset mode only)
```

where `A` is the perturbation/test group and `B` the reference (control) group.

> **Interaction mode note**: in `--stat_mode interaction`, `logFC` is the log2
> focal/other odds-ratio change, not `log2(PSI_A/PSI_B)`. The remainder of this
> document describes the offset-mode relationship.

---

## 2. When `logFC` is a reported-PSI log-ratio

In `splice_plus_retained` mode, PSI and the edgeR offset use the same
denominator, so the **log-ratio of PSI factorizes exactly**:

```
log2(PSI_A / PSI_B)
   = log2(count_A / count_B)  -  log2(denominator_A / denominator_B)
   = (raw numerator-count logFC)  -  (denominator logFC)
```

The second term is exactly the **offset** edgeR puts in the GLM:

```
log(expected count) = group_effect + log(denominator)
                                      └──── offset ────┘
```

Subtracting the offset is the same as dividing the count by the denominator, so
the fitted **group effect is the log-fold-change of PSI**, not of the raw count.
That is why the reported `logFC` is an *offset-adjusted intron-usage* log
fold-change rather than an expression fold-change. In `splice_vs_rest` mode,
`logFC` instead describes usage relative to `D^svr` and is not directly
convertible to the reported SPR-based PSI or delta PSI.

> **In SPR offset mode, `logFC` estimates `log2(PSI_A / PSI_B)`.** It is not
> numerically identical to the log-ratio of the reported `mean_PSI` columns,
> because edgeR uses GLM-fitted, dispersion-weighted, replicate-aware values
> (a ratio of fitted means), whereas `mean_PSI` is the mean of per-sample
> ratios. The two diverge most when counts are low or a denominator PSI is near
> zero (mean-of-ratios ≠ ratio-of-means). Use `logFC` for the model estimate and
> `delta_PSI` for the observed magnitude.

---

## 3. Converting between `logFC` and `delta_PSI`

For SPR offset mode, treating `logFC` as `log2(PSI_A / PSI_B)`, and writing the
baseline (control) proportion as `PSI_B`:

```
PSI_A     = PSI_B · 2^logFC
delta_PSI = PSI_A - PSI_B = PSI_B · (2^logFC - 1)
logFC     = log2( 1 + delta_PSI / PSI_B )
```

So in SPR offset mode the two are interconvertible **only with a third quantity, the baseline
`PSI_B`** (`control_mean_PSI` in the results). The pipeline reports it, so you
can convert row by row — but you cannot go from one to the other on its own.
This conversion does not apply to SVR model `logFC`.

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
  that carry a meaningful share of the SPR denominator.

The FDR/`logFC` side answers *whether* usage changed; `delta_PSI` answers *how
much* of the local splicing output actually moved. They are complementary, not
redundant.

### Filtering and FDR in the current pipeline

`--min_delta_psi` is a **pre-test** filter in the current pipeline. Introns that
do not pass the configured absolute `delta_PSI` threshold are removed before the
selected model is fit, and Benjamini-Hochberg FDR is computed over the remaining
tested intron set. The pipeline does not recompute FDR after model fitting.

`edgeR_results.significant_introns.tsv` is produced from the tested result table
using edgeR's `FDR`, the configured `--fdr_threshold`, and the configured
`--min_logFC`. If `--min_delta_psi` was nonzero, all rows in the tested table
already passed that prefilter.

---

## 5. Summary

| | `logFC` | `delta_PSI` |
|---|---|---|
| scale | model-denominator-adjusted log2 ratio | arithmetic difference of SPR PSI |
| range | `(-inf, inf)` | `[-1, 1]` |
| source | edgeR GLM coefficient (offset-adjusted) | mean PSI difference |
| sensitive to | small switches at **low** PSI | switches among **abundant** introns |
| convert to the other | needs baseline `PSI_B`, and only for SPR offset mode | needs baseline `PSI_B`, and only for SPR offset mode |

In SPR offset mode, `logFC` and `delta_PSI` are linked by
`delta_PSI = PSI_B·(2^logFC − 1)` and require the control-group baseline
`PSI_B` for conversion. SVR-model and interaction `logFC` values are not
convertible to the reported SPR delta PSI with this identity.
