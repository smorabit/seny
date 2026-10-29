#' seny: sanity checks for differential expression in noisy data
#'
#' For every cluster-by-contrast comparison in a multi-sample single-cell
#' dataset, seny decides whether a pseudobulk differential expression test is
#' Adequate, Underpowered, Caution or Invalid, and explains why. The diagnostics
#' are analytic and run on the user's own data: design validity and confounding
#' (layer 1), pseudobulk quality (layer 2), and variance structure plus
#' sensitivity (layer 4), the last giving a per-gene minimum detectable effect
#' from a negative-binomial variance model. No per-dataset simulation.
#'
#' The name is a backronym -- Sanity checks for differential gene Expression in
#' NoisY data -- and also the Catalan word for common sense.
#'
#' The main entry point is [assessAll()], which chains [pbSetup()],
#' [assessDesign()], [assessPseudobulk()], [estimateVariance()],
#' [estimateSensitivity()] and [callComparisons()].
#'
#' @keywords internal
'_PACKAGE'
