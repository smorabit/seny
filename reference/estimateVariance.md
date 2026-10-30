# Layer 4 (variance): null-model dispersion, decomposition and per-gene SE

For every comparison, fits the chosen DE method and returns the per-gene
standard error of the two-group log2 fold-change contrast (used by
[`estimateSensitivity()`](https://smorabit.github.io/seny/reference/estimateSensitivity.md)).
The SE uses the dispersion estimated with the condition in the design,
as the user's test sees it, evaluated at the condition-blind
(null-model) fitted means. Power is always computed at the target
effect, never the observed one. A null-model dispersion would absorb
real effects as noise and understate power exactly where DE is strong.
The variance decomposition stays condition-blind. For the count methods
(`edgeR`, `deseq2`) it also splits the pseudobulk negative-binomial
dispersion `phi_total` into a cell-heterogeneity term (`theta`/n, with
the cell-level `theta` from
[`pbSetup()`](https://smorabit.github.io/seny/reference/pbSetup.md)) and
a between-donor term `phi_donor` (moment estimate
`max(0, phi_total - theta*mean(1/n))`), and reports the depth/cell/donor
variance shares with a regime label. The `voom` method uses limma's
precision-weight variance model, so it has no NB decomposition (those
columns are `NA`).

## Usage

``` r
estimateVariance(x, method = c("edgeR", "deseq2", "voom", "dreamlet"))
```

## Arguments

- x:

  A
  [S4Vectors::SimpleList](https://rdrr.io/pkg/S4Vectors/man/SimpleList-class.html)
  from
  [`pbSetup()`](https://smorabit.github.io/seny/reference/pbSetup.md).

- method:

  One of `"edgeR"`, `"deseq2"`, `"voom"` (or the reserved `"dreamlet"`).
  Default `"edgeR"`.

## Value

A list with `summary` (a `data.frame`, one row per comparison, with the
`method`, the NB decomposition and regime) and `pergene` (a named list
of per-comparison fits carrying `se_log2` and `df` for sensitivity).

## Details

The method should match the differential-expression method the user
intends to run. `edgeR` and `deseq2` require sum aggregation (integer
counts); `voom` works on either sum (voom) or mean (limma-trend)
pseudobulks. `dreamlet` is reserved for mixed designs and is not yet
implemented.
