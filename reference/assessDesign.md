# Layer 1: design validity and confounding

For every comparison, measures whether the two-group contrast is
estimable (full design rank with residual degrees of freedom) and how
strongly the condition is confounded with each candidate covariate: the
variance inflation factor, bias-corrected Cramer's V for categorical
covariates, and the standardised mean difference for numeric ones. A
rank-deficient or fully nested design is flagged non-estimable, which
downstream becomes an Invalid call.

## Usage

``` r
assessDesign(x)
```

## Arguments

- x:

  A
  [S4Vectors::SimpleList](https://rdrr.io/pkg/S4Vectors/man/SimpleList-class.html)
  from
  [`pbSetup()`](https://smorabit.github.io/seny/reference/pbSetup.md).

## Value

A `data.frame`, one row per comparison, with design and confounding
metrics (`estimable`, `resid_df`, `nA`, `nB`, `max_vif`, `max_cramerV`,
`max_smd`).
