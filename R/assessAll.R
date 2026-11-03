#' Assess every pseudobulk comparison end to end
#'
#' The main entry point. Runs the pipeline on a single-cell dataset:
#' [pbSetup()] -> [assessDesign()] -> [assessPseudobulk()] ->
#' [estimateVariance()] -> [estimateSensitivity()] -> [assessComposition()] ->
#' [callComparisons()], and
#' returns the comparison table with an Adequate / Underpowered / Caution /
#' Invalid call per cluster-by-contrast comparison.
#'
#' @inheritParams pbSetup
#' @param method DE method whose variance model the call should match, passed to
#'   [estimateVariance()]: `"edgeR"`, `"deseq2"` or `"voom"`. `"edgeR"`/`"deseq2"`
#'   require `aggregation = "sum"`; `aggregation = "mean"` requires `"voom"`.
#' @param targetLFC Primary target absolute log2 fold-change driving the call.
#' @param power,fdr,pi1,alphaFixed Sensitivity settings passed to
#'   [estimateSensitivity()].
#' @param thresholds A list from [defaultThresholds()].
#' @param compReference Source of the per-cell subtype expression in
#'   [assessComposition()]: `"reference"` (default) or `"pooled"`.
#' @param thetaReps Subsample repeats for the cell-level theta fit, passed to
#'   [pbSetup()]; `0` skips it (the call is unchanged, the variance
#'   decomposition is `NA`).
#' @return A [S4Vectors::DataFrame], one row per comparison, with the metrics,
#'   the `call` and a `reason`. The [pbSetup()] object is attached as
#'   `attr(, "pb")`, the per-gene table as `attr(, "pergene")` and the full
#'   [assessComposition()] result as `attr(, "composition")`. Composition beyond
#'   the cluster share needs two nested annotation levels in `annotations`
#'   (e.g. broad and fine).
#' @export
#' @examples
#' \dontrun{
#' res <- assessAll(sce, sample = "donor", group = "condition",
#'                  annotations = "cell_type")
#' table(res$call)
#' }
assessAll <- function(sce, sample, group, annotations, covariates = character(0),
                      contrasts = NULL, minCells = 10, aggregation = c('sum', 'mean'),
                      assayName = 'counts', method = c('edgeR', 'deseq2', 'voom'),
                      targetLFC = 1, power = 0.8, fdr = 0.05,
                      pi1 = 'discoveries', alphaFixed = NULL, thresholds = defaultThresholds(), thetaReps = 3,
                      compReference = c('reference', 'pooled'),
                      BPPARAM = BiocParallel::SerialParam()) {
    method <- match.arg(method)
    x <- pbSetup(sce, sample = sample, group = group, annotations = annotations,
                 covariates = covariates, contrasts = contrasts, minCells = minCells,
                 aggregation = aggregation, assayName = assayName, thetaReps = thetaReps,
                 BPPARAM = BPPARAM)
    design <- assessDesign(x)
    pseudobulk <- assessPseudobulk(x)
    variance <- estimateVariance(x, method = method)
    sensitivity <- estimateSensitivity(variance, targetLFC = targetLFC, power = power,
                                       fdr = fdr, pi1 = pi1, alphaFixed = alphaFixed)
    composition <- assessComposition(x, variance, sensitivity, power = power, fdr = fdr,
                                     reference = match.arg(compReference), thresholds = thresholds)
    res <- callComparisons(design, pseudobulk, variance$summary, sensitivity, thresholds,
                           composition = composition$summary)
    attr(res, 'pb') <- x
    attr(res, 'composition') <- composition
    attr(res, 'pergene') <- collect_pergene(variance$pergene)
    res
}

#' Flatten the per-comparison fits into one per-gene table for the plots
#' @noRd
collect_pergene <- function(pergene) {
    rows <- lapply(names(pergene), function(k) {
        f <- pergene[[k]]
        data.frame(comparison = k, gene = f$genes, acpm = f$acpm,
                   se_log2 = f$se_log2, df = f$df, stringsAsFactors = FALSE)
    })
    if (length(rows) == 0) return(NULL)
    do.call(rbind, rows)
}
