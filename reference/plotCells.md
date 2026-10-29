# Cells per sample by condition

The cells contributed per sample, split by condition, as a beeswarm over
a boxplot and faceted by cluster. When cell counts track the condition
the pseudobulks are heteroscedastic, which is the risk behind a
cell-count imbalance Caution.

## Usage

``` r
plotCells(res)
```

## Arguments

- res:

  The result of
  [`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md)
  (its `pb` setup object is required).

## Value

A `ggplot`, faceted by cluster (and annotation set when there is more
than one).
