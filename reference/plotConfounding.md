# Confounding strip: condition-covariate association per comparison

A Layer-1 lollipop of how strongly the condition is tangled with its
covariates, with the Caution threshold drawn in and points coloured by
the call. A comparison past the line cannot cleanly separate the
condition from the covariate.

## Usage

``` r
plotConfounding(res, metric = c("cramerV", "vif"))
```

## Arguments

- res:

  The result of
  [`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md)
  or
  [`callComparisons()`](https://smorabit.github.io/seny/reference/callComparisons.md).

- metric:

  `"cramerV"` (bias-corrected Cramer's V, categorical covariates) or
  `"vif"` (variance inflation factor). Default `"cramerV"`.

## Value

A `ggplot`, faceted by annotation set when there is more than one.
