# Power vs effect size

Mean power across a sweep of target log2 fold-changes, one line per
cluster, with the 0.8 reference. Shows where the call would change at a
smaller or larger target effect.

## Usage

``` r
plotPower(res, targets = seq(0.25, 2, by = 0.05))
```

## Arguments

- res:

  The result of
  [`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md)
  (its per-gene `pergene` table is required).

- targets:

  Numeric vector of target absolute log2 fold-changes to sweep.

## Value

A `ggplot`, faceted by annotation set when there is more than one.
