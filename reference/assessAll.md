# Assess every pseudobulk comparison end to end

The main entry point. Runs the pipeline on a single-cell dataset:
[`pbSetup()`](https://smorabit.github.io/seny/reference/pbSetup.md) -\>
[`assessDesign()`](https://smorabit.github.io/seny/reference/assessDesign.md)
-\>
[`assessPseudobulk()`](https://smorabit.github.io/seny/reference/assessPseudobulk.md)
-\>
[`estimateVariance()`](https://smorabit.github.io/seny/reference/estimateVariance.md)
-\>
[`estimateSensitivity()`](https://smorabit.github.io/seny/reference/estimateSensitivity.md)
-\>
[`callComparisons()`](https://smorabit.github.io/seny/reference/callComparisons.md),
and returns the comparison table with an Adequate / Underpowered /
Caution / Invalid call per cluster-by-contrast comparison.

## Usage

``` r
assessAll(
  sce,
  sample,
  group,
  annotations,
  covariates = character(0),
  contrasts = NULL,
  minCells = 10,
  aggregation = c("sum", "mean"),
  assayName = "counts",
  targetLFC = 1,
  power = 0.8,
  fdr = 0.05,
  pi1 = 0.05,
  thresholds = defaultThresholds(),
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- sce:

  A
  [SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  or
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  with a raw-count assay (non-negative integers).

- sample:

  Character scalar; `colData` column giving the biological sample
  (donor) each cell belongs to.

- group:

  Character scalar; `colData` column giving the condition compared.

- annotations:

  Character vector of `colData` columns, each a clustering (annotation
  level) to pseudobulk, e.g. `c("cell_type_broad")`.

- covariates:

  Character vector of `colData` columns to carry as candidate
  confounders. Default none.

- contrasts:

  Optional list of length-2 character vectors `c(A, B)` naming `group`
  levels to compare (B vs A). Default `NULL`: all non-reference levels
  versus the first level.

- minCells:

  Integer; pseudobulks with fewer contributing cells are dropped.
  Default 10.

- aggregation:

  Either `"sum"` (counts, for the negative-binomial model) or `"mean"`
  (reserved for the limma path). Default `"sum"`.

- assayName:

  Name of the raw-count assay in `sce`. Default `"counts"`.

- targetLFC:

  Primary target absolute log2 fold-change driving the call.

- power, fdr, pi1:

  Sensitivity settings passed to
  [`estimateSensitivity()`](https://smorabit.github.io/seny/reference/estimateSensitivity.md).

- thresholds:

  A list from
  [`defaultThresholds()`](https://smorabit.github.io/seny/reference/defaultThresholds.md).

- BPPARAM:

  A `BiocParallelParam` object for the per-cluster theta estimation.
  Default
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html).

## Value

A
[S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html),
one row per comparison, with the metrics, the `call` and a `reason`. The
[`pbSetup()`](https://smorabit.github.io/seny/reference/pbSetup.md)
object is attached as `attr(, "pb")` for inspection.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- assessAll(sce, sample = "donor", group = "condition",
                 annotations = "cell_type")
table(res$call)
} # }
```
