#' Assess every pseudobulk comparison end to end
#'
#' The main entry point. Runs the pipeline on a single-cell dataset:
#' [pbSetup()] -> [assessDesign()] -> [assessPseudobulk()] ->
#' [estimateVariance()] -> [estimateSensitivity()] -> [callComparisons()], and
#' returns the comparison table with an Adequate / Underpowered / Caution /
#' Invalid call per cluster-by-contrast comparison.
#'
#' @inheritParams pbSetup
#' @param targetLFC Primary target absolute log2 fold-change driving the call.
#' @param power,fdr,pi1 Sensitivity settings passed to [estimateSensitivity()].
#' @param thresholds A list from [defaultThresholds()].
#' @return A [S4Vectors::DataFrame], one row per comparison, with the metrics,
#'   the `call` and a `reason`. The [pbSetup()] object is attached as
#'   `attr(, "pb")` for inspection.
#' @export
#' @examples
#' \dontrun{
#' res <- assessAll(sce, sample = "donor", group = "condition",
#'                  annotations = "cell_type")
#' table(res$call)
#' }
assessAll <- function(sce, sample, group, annotations, covariates = character(0),
                      contrasts = NULL, minCells = 10, aggregation = c('sum', 'mean'),
                      assayName = 'counts', targetLFC = 1, power = 0.8, fdr = 0.05,
                      pi1 = 0.05, thresholds = defaultThresholds(),
                      BPPARAM = BiocParallel::SerialParam()) {
    x <- pbSetup(sce, sample = sample, group = group, annotations = annotations,
                 covariates = covariates, contrasts = contrasts, minCells = minCells,
                 aggregation = aggregation, assayName = assayName, BPPARAM = BPPARAM)
    design <- assessDesign(x)
    pseudobulk <- assessPseudobulk(x)
    variance <- estimateVariance(x)
    sensitivity <- estimateSensitivity(variance, targetLFC = targetLFC, power = power,
                                       fdr = fdr, pi1 = pi1)
    res <- callComparisons(design, pseudobulk, variance$summary, sensitivity, thresholds)
    attr(res, 'pb') <- x
    res
}
