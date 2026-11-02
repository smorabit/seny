# Composition of one comparison: subtype shares and composition-driven DE

Two panels for a comparison whose cluster has a finer annotation nested
in it (see
[`assessComposition()`](https://smorabit.github.io/seny/reference/assessComposition.md)).
Left: each subtype's share of the cluster's cells per sample, by
condition; a shift here makes the subtypes' markers look DE in the
cluster's pseudobulk. Right: per gene, the observed log2 fold change of
the cluster's own DE test against the change the subtype shift alone
predicts. Genes the test calls DE are coloured by their attribution;
genes on the diagonal are explained by composition.

## Usage

``` r
plotComposition(res, comparison, nLabel = 10)
```

## Arguments

- res:

  The result of
  [`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md)
  run with two nested annotation levels.

- comparison:

  A comparison key from `res$comparison`.

- nLabel:

  Number of composition-driven genes to label, largest observed change
  first. Default 10.

## Value

A `patchwork` of two `ggplot`s.
