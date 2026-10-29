# Traffic-light matrix of adequacy calls

Renders the comparison table from
[`callComparisons()`](https://smorabit.github.io/seny/reference/callComparisons.md)
(or
[`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md))
as a clusters-by-contrasts grid, one plot per annotation set, each tile
coloured by the comparison's call (Adequate / Underpowered / Caution /
Invalid). The call reasons stay in the table; the matrix is the
at-a-glance view.

## Usage

``` r
plotAdequacy(res)
```

## Arguments

- res:

  The
  [S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  returned by
  [`callComparisons()`](https://smorabit.github.io/seny/reference/callComparisons.md)
  or
  [`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md),
  one row per comparison, with `annotation`, `cluster`, `A`, `B` and
  `call`.

## Value

A named list of `ggplot` objects, one per annotation set (a single-level
input returns a one-element list).
