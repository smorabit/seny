# Assemble the comparison table and make the adequacy call

Merges the layer metrics into one tidy table, one row per comparison,
and assigns each comparison a call: **Invalid** (the contrast is not
estimable, or too few samples per group), **Caution** (estimable but the
condition is confounded with a covariate, or cell-count imbalance tracks
the condition), **Underpowered** (valid but the expected recall at the
target effect is below `thresholds$sensitivity`, or too few testable
genes), or **Adequate**. Severity order is Invalid \> Caution \>
Underpowered \> Adequate; the `reason` column lists every flag so the
single call never hides a second problem.

## Usage

``` r
callComparisons(
  design,
  pseudobulk,
  variance,
  sensitivity,
  thresholds = defaultThresholds()
)
```

## Arguments

- design, pseudobulk, variance, sensitivity:

  The outputs of
  [`assessDesign()`](https://smorabit.github.io/seny/reference/assessDesign.md),
  [`assessPseudobulk()`](https://smorabit.github.io/seny/reference/assessPseudobulk.md),
  [`estimateVariance()`](https://smorabit.github.io/seny/reference/estimateVariance.md)
  (its `$summary`) and
  [`estimateSensitivity()`](https://smorabit.github.io/seny/reference/estimateSensitivity.md).

- thresholds:

  A list from
  [`defaultThresholds()`](https://smorabit.github.io/seny/reference/defaultThresholds.md).

## Value

A
[S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html),
one row per comparison, with all metrics, the `call` and a `reason`
string.
