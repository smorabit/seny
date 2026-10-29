# Layer 2: pseudobulk quality

For every comparison, summarises the pseudobulks the test will run on:
the cells contributed per sample and their inequality, whether cell
counts or library sizes track the condition (which breaks the
equal-variance assumption), and how many genes are testable.
Specifically, per comparison: cells per sample (`min_cells`,
`median_cells`, `gini_cells`, and `frac_below_min`, the fraction of the
cluster's samples dropped for having fewer than `minCells` cells); group
imbalance (`cells_ratio`, the smaller group median over the larger, and
`cells_cond_cor`, the correlation of log cells with condition);
`libsize_cond_cor` (correlation of log library size with condition); and
testable genes (`n_testable` passing
[edgeR::filterByExpr](https://rdrr.io/pkg/edgeR/man/filterByExpr.html)
under the comparison's design, and `frac_testable`). The testable-gene
count is the denominator for the sensitivity fraction in
[`estimateSensitivity()`](https://smorabit.github.io/seny/reference/estimateSensitivity.md).
Outlier-sample detection is part of layer 3, not here.

## Usage

``` r
assessPseudobulk(x)
```

## Arguments

- x:

  A
  [S4Vectors::SimpleList](https://rdrr.io/pkg/S4Vectors/man/SimpleList-class.html)
  from
  [`pbSetup()`](https://smorabit.github.io/seny/reference/pbSetup.md).

## Value

A `data.frame`, one row per comparison, with `min_cells`,
`median_cells`, `gini_cells`, `frac_below_min`, `cells_ratio`,
`cells_cond_cor`, `libsize_cond_cor`, `n_testable` and `frac_testable`.
