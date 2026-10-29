# Layer 4 (variance): null-model dispersion and its decomposition

For every comparison, fits the null model (design without the condition
term) with edgeR to get the pseudobulk negative-binomial dispersion
`phi_total`, then splits it into a cell-heterogeneity term (`theta`/n,
with the cell-level `theta` estimated at
[`pbSetup()`](https://smorabit.github.io/seny/reference/pbSetup.md)) and
a between-donor term `phi_donor` via the moment estimate
`phi_donor = max(0, phi_total - theta * mean(1/n))`. Reports the
per-gene fits needed by
[`estimateSensitivity()`](https://smorabit.github.io/seny/reference/estimateSensitivity.md)
and a per-comparison summary of the variance shares (depth, cell, donor)
with a regime label.

## Usage

``` r
estimateVariance(x)
```

## Arguments

- x:

  A
  [S4Vectors::SimpleList](https://rdrr.io/pkg/S4Vectors/man/SimpleList-class.html)
  from
  [`pbSetup()`](https://smorabit.github.io/seny/reference/pbSetup.md).

## Value

A list with `summary` (a `data.frame`, one row per comparison) and
`pergene` (a named list of per-comparison fit objects for sensitivity).
