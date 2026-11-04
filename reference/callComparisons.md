# Assemble the comparison table and make the adequacy call

Merges the layer metrics into one tidy table, one row per comparison,
and assigns each comparison a call: **Invalid** (the contrast is not
estimable, or too few samples per group), **Caution** (estimable but the
condition is confounded with a covariate, cell-count imbalance tracks
the condition under DESeq2 (whose FDR control it breaks; with edgeR or
voom the imbalance is reported in the reason but does not change the
call), or a subtype shift within the cluster would by itself make enough
genes look DE; see
[`assessComposition()`](https://smorabit.github.io/seny/reference/assessComposition.md)),
**Underpowered** (valid but the expected recall at the target effect is
below `thresholds$sensitivity`, or too few testable genes), or
**Adequate**. Severity order is Invalid \> Caution \> Underpowered \>
Adequate; the `reason` column lists every flag so the single call never
hides a second problem.

## Usage

``` r
callComparisons(
  design,
  pseudobulk,
  variance,
  sensitivity,
  thresholds = defaultThresholds(),
  composition = NULL
)
```

## Arguments

- design, pseudobulk, variance, sensitivity:

  The outputs of
  [`assessDesign()`](https://smorabit.github.io/seny/reference/assessDesign.md),
  [`assessPseudobulk()`](https://smorabit.github.io/seny/reference/assessPseudobulk.md),
  [`estimateVariance()`](https://smorabit.github.io/seny/reference/estimateVariance.md)
  (its `$summary`) and
  [`estimateSensitivity()`](https://smorabit.github.io/seny/reference/estimateSensitivity.md).

- thresholds:

  A list from
  [`defaultThresholds()`](https://smorabit.github.io/seny/reference/defaultThresholds.md).

- composition:

  Optional `$summary` of
  [`assessComposition()`](https://smorabit.github.io/seny/reference/assessComposition.md).
  `NULL` (default) skips the composition columns and flag.

## Value

A
[S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html),
one row per comparison, with all metrics, the `call`, a `reason` string
and `flags`, the same flags as short codes separated by `;`
(`confounded`, `imbalance`, `composition`, `few_testable`,
`underpowered`, `too_few_samples`, `not_estimable`; empty when
adequate).

## Examples

``` r
sce <- readRDS(system.file("extdata", "example_sce.rds", package = "seny"))
x <- pbSetup(sce, sample = "donor", group = "condition",
             annotations = "cell_type", thetaReps = 0)
variance <- estimateVariance(x)
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
sens <- estimateSensitivity(variance, targetLFC = 1)
res <- callComparisons(assessDesign(x), assessPseudobulk(x), variance$summary, sens)
res[, c("comparison", "call", "reason")]
#> DataFrame with 3 rows and 3 columns
#>              comparison         call                 reason
#>             <character>     <factor>            <character>
#> 1 cell_type:CD4T:B_vs_A Adequate                   adequate
#> 2 cell_type:CD8T:B_vs_A Adequate     cell-count imbalance..
#> 3   cell_type:NK:B_vs_A Underpowered underpowered (expect..
```
