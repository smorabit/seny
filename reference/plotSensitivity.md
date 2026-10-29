# Detectable effect vs expression

For every gene, the minimum detectable log2 fold-change against its
expression, with the target effect drawn in and one panel per cluster.
The share of genes below the line is the sensitivity fraction; genes
above it are out of reach at the target.

## Usage

``` r
plotSensitivity(res)
```

## Arguments

- res:

  The result of
  [`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md)
  (its per-gene `pergene` table is required).

## Value

A `ggplot`, faceted by cluster (and annotation set when there is more
than one).
