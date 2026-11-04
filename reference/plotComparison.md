# Comparison card: a per-comparison deep dive

Assembles one comparison into a single card: the call and its reason,
the variance attribution as a donut, the sensitivity curve, and cells
per sample. This is the per-comparison unit behind the eventual report.

## Usage

``` r
plotComparison(res, comparison)
```

## Arguments

- res:

  The result of
  [`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md)
  (its `pb` and per-gene `pergene` are required).

- comparison:

  A comparison key from `res$comparison`.

## Value

A `patchwork` composition of `ggplot`s.

## Examples

``` r
sce <- readRDS(system.file("extdata", "example_sce.rds", package = "seny"))
res <- assessAll(sce, sample = "donor", group = "condition",
                 annotations = "cell_type", targetLFC = 1, thetaReps = 0)
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
plotComparison(res, "cell_type:CD4T:B_vs_A")
```
