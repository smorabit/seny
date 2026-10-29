# Variance regime ternary

The depth / cell / donor decomposition as a single simplex: each
comparison is a point among the three variance terms, coloured by its
call and sized by its median detectable effect. Position shows where the
noise is, size shows how hard detection is. The axes are gridded so the
shares can be read off directly. Comparisons fit with the voom method
have no decomposition and are dropped.

## Usage

``` r
plotRegime(res)
```

## Arguments

- res:

  The result of
  [`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md)
  or
  [`callComparisons()`](https://smorabit.github.io/seny/reference/callComparisons.md).
  Comparisons from all annotation sets are pooled into one triangle;
  subset `res` first for a single annotation.

## Value

A `ggplot`.
