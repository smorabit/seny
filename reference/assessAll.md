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
[`assessComposition()`](https://smorabit.github.io/seny/reference/assessComposition.md)
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
  method = c("edgeR", "deseq2", "voom"),
  targetLFC = 1,
  power = 0.8,
  fdr = 0.05,
  pi1 = "discoveries",
  alphaFixed = NULL,
  thresholds = defaultThresholds(),
  thetaReps = 0,
  compReference = c("reference", "pooled"),
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- sce:

  A
  [SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  or
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  with a raw-count assay (non-negative integers). The assay may be a
  base or sparse matrix, or an out-of-memory matrix such as a
  `DelayedMatrix` or an HDF5-backed matrix (for atlases beyond one
  sparse matrix's 2^31 nonzeros): it is then checked and summed a block
  of cells at a time (`options(seny.block_cells = ...)`, default 1e5).

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

- method:

  DE method whose variance model the call should match, passed to
  [`estimateVariance()`](https://smorabit.github.io/seny/reference/estimateVariance.md):
  `"edgeR"`, `"deseq2"` or `"voom"`. `"edgeR"`/`"deseq2"` require
  `aggregation = "sum"`; `aggregation = "mean"` requires `"voom"`.

- targetLFC:

  Primary target absolute log2 fold-change driving the call.

- power, fdr, pi1, alphaFixed:

  Sensitivity settings passed to
  [`estimateSensitivity()`](https://smorabit.github.io/seny/reference/estimateSensitivity.md).

- thresholds:

  A list from
  [`defaultThresholds()`](https://smorabit.github.io/seny/reference/defaultThresholds.md).

- thetaReps:

  Subsample repeats for the cell-level theta fit, passed to
  [`pbSetup()`](https://smorabit.github.io/seny/reference/pbSetup.md);
  `0` skips it (the call is unchanged, the variance decomposition is
  `NA`).

- compReference:

  Source of the per-cell subtype expression in
  [`assessComposition()`](https://smorabit.github.io/seny/reference/assessComposition.md):
  `"reference"` (default) or `"pooled"`.

- BPPARAM:

  A `BiocParallelParam` object for the per-cluster theta estimation.
  Default
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html).

## Value

A
[S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html),
one row per comparison, with the metrics, the `call`, a `reason` and
`flags`. The
[`pbSetup()`](https://smorabit.github.io/seny/reference/pbSetup.md)
object is attached as `attr(, "pb")`, the per-gene table as
`attr(, "pergene")`, the full
[`assessComposition()`](https://smorabit.github.io/seny/reference/assessComposition.md)
result as `attr(, "composition")` and the settings of the run (method,
target, thresholds, seny version) as `attr(, "settings")`, which
[`reportSeny()`](https://smorabit.github.io/seny/reference/reportSeny.md)
reads. Composition beyond the cluster share needs two nested annotation
levels in `annotations` (e.g. broad and fine).

## Examples

``` r
sce <- readRDS(system.file("extdata", "example_sce.rds", package = "seny"))
#> Loading required namespace: SingleCellExperiment
res <- assessAll(sce, sample = "donor", group = "condition",
                 annotations = "cell_type", targetLFC = 1, thetaReps = 0)
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
table(res$call)
#> 
#>     Adequate Underpowered      Caution      Invalid 
#>            2            1            0            0 
res[, c("cluster", "call", "reason")]
#> DataFrame with 3 rows and 3 columns
#>       cluster         call                 reason
#>   <character>     <factor>            <character>
#> 1        CD4T Adequate                   adequate
#> 2        CD8T Adequate     cell-count imbalance..
#> 3          NK Underpowered underpowered (expect..
```
