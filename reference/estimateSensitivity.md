# Layer 4 (sensitivity): minimum detectable effect and power

From the null-model variance fits
([`estimateVariance()`](https://smorabit.github.io/seny/reference/estimateVariance.md)),
computes the analytic standard error of the two-group log2 fold-change
contrast per gene (`SE^2 = 1/sum_A w + 1/sum_B w`,
`w = 1/(1/mu + phi_total)`), the minimum detectable effect at the target
power and FDR, and the power at the primary target effect. Summarises to
the sensitivity fraction: the share of testable genes detectable at
`targetLFC`. Uncalibrated in the MVP (`calibrate = FALSE`); the
empirical calibration from E2/E7 is applied in a later phase.

## Usage

``` r
estimateSensitivity(
  variance,
  targetLFC = 1,
  power = 0.8,
  fdr = 0.05,
  pi1 = 0.05,
  calibrate = FALSE
)
```

## Arguments

- variance:

  The list returned by
  [`estimateVariance()`](https://smorabit.github.io/seny/reference/estimateVariance.md).

- targetLFC:

  Primary target absolute log2 fold-change for the call. Default 1.

- power:

  Target power. Default 0.8.

- fdr:

  Target FDR. Default 0.05.

- pi1:

  Assumed fraction of true positives (maps FDR to a per-gene threshold,
  Jung 2005). Default 0.05.

- calibrate:

  Reserved; must be `FALSE` in the MVP.

## Value

A `data.frame`, one row per comparison, with `target_lfc`,
`sensitivity_fraction`, `median_mde` and `median_power`.
