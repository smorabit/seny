# Variance attribution bars

Stacked depth / cell / donor variance shares per comparison: where each
comparison's per-sample noise comes from. A donor-dominated bar means
more cells or reads will not help, only more donors. Comparisons fit
with the voom method have no negative-binomial decomposition and are
dropped.

## Usage

``` r
plotVariance(res)
```

## Arguments

- res:

  The result of
  [`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md)
  or
  [`callComparisons()`](https://smorabit.github.io/seny/reference/callComparisons.md).

## Value

A `ggplot`, faceted by annotation set when there is more than one.
