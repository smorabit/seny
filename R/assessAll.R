#' Assess every pseudobulk comparison end to end
#'
#' The main entry point. Runs the pipeline on a single-cell dataset:
#' [pbSetup()] -> [assessDesign()] -> [assessPseudobulk()] ->
#' [estimateVariance()] -> [estimateSensitivity()] -> [callComparisons()], and
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
                      assayName = 'counts', method = c('edgeR', 'deseq2', 'voom'),
                      targetLFC = 1, power = 0.8, fdr = 0.05,
                      pi1 = 0.05, alphaFixed = 0.001, thresholds = defaultThresholds(), thetaReps = 3,
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
    res <- callComparisons(design, pseudobulk, variance$summary, sensitivity, thresholds)
    attr(res, 'pb') <- x
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
