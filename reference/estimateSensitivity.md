# Layer 4 (sensitivity): minimum detectable effect and power

From the variance fits
([`estimateVariance()`](https://smorabit.github.io/seny/reference/estimateVariance.md)),
computes the analytic standard error of the two-group log2 fold-change
contrast per gene (`SE^2 = 1/sum_A w + 1/sum_B w`, `w = 1/(1/mu + phi)`,
phi the dispersion with the condition in the design), the minimum
detectable effect at the target power and FDR, and the power at the
primary target effect. Summarises to the sensitivity fraction: the share
of testable genes detectable at `targetLFC`. Uncalibrated in the MVP
(`calibrate = FALSE`); the empirical calibration from E2/E7 is applied
in a later phase.

## Usage

``` r
estimateSensitivity(
  variance,
  targetLFC = 1,
  power = 0.8,
  fdr = 0.05,
  pi1 = "discoveries",
  alphaFixed = NULL,
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

  Fraction of true positives that maps the FDR to a per-gene threshold:
  `"discoveries"` (default) estimates it per comparison as the fraction
  of genes the comparison's own DE test calls at `fdr`; a number fixes
  it for every comparison (falls back to 0.05 for a comparison whose fit
  has no p-values).

- alphaFixed:

  Optional fixed per-gene significance level. `NULL` (default) uses the
  `pi1`-mapped level alone; a number makes the call use the stricter of
  the two (the larger MDE).

- calibrate:

  Reserved; must be `FALSE` in the MVP.

## Value

A `data.frame`, one row per comparison, with `target_lfc`, the `pi1` and
per-gene `alpha` used, `sensitivity_fraction`, `median_mde`,
`median_mde_fdr`, `median_mde_fixed` (`NA` without `alphaFixed`) and
`median_power`.

## Details

The FDR target is mapped to a per-gene significance level through the
fraction of true positives `pi1` (Jung 2005), capped at `fdr` (BH never
tests a gene at more than the FDR level). BH's realised per-gene
threshold scales with the fraction of genes it calls, so no fixed `pi1`
suits every contrast: a sparse contrast (male vs female in blood,
realised threshold ~1e-4) makes a fixed 5% optimistic, a dense one makes
it conservative. By default (`pi1 = "discoveries"`) each comparison uses
its own fraction of BH discoveries, `R / m` from the DE test of the
chosen method at `fdr` (floored at `1 / m`), so power is evaluated at
the threshold BH actually uses there. On OneK1K with 0.5-20\\ of genes
injected as DE and on male vs female it was the best calibrated of the
rules tested (E4b).
