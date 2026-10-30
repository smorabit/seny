#' Default decision thresholds for the adequacy call
#'
#' Provisional values used by [callComparisons()]. The final numbers are to be
#' derived from the validation experiments (E1-E8); until then these are sensible
#' placeholders, documented so a user can override them.
#'
#' @return A named list: `cramerV` (condition-covariate association above which a
#'   comparison is flagged Caution), `vif` (variance inflation factor threshold),
#'   `sensitivity` (sensitivity fraction at or above which a comparison is
#'   adequately powered), `minSamples` (fewer samples per group is Invalid),
#'   `cellsCondCor` (absolute correlation of log cells per sample with condition
#'   above which cell-count imbalance is flagged Caution), `cellsCondAlpha`
#'   (the correlation must also exceed the critical value at this level for the
#'   comparison's sample size, so chance correlations in small designs are not
#'   flagged: the cutoff is `max(cellsCondCor, r_crit(n))`, about 0.81 at 6
#'   samples and `cellsCondCor` from ~16 samples), and `minTestable` (fewer
#'   testable genes is Underpowered). The layer-2 values, like the others, are
#'   provisional placeholders until E1-E8.
#' @export
defaultThresholds <- function() {
    list(cramerV = 0.3, vif = 5, sensitivity = 0.5, minSamples = 2,
         cellsCondCor = 0.5, cellsCondAlpha = 0.05, minTestable = 200)
}

#' Cutoff for the cell-count-vs-condition correlation at n samples: the larger of
#' the fixed floor and the two-sided critical Pearson r at level alpha
#' @noRd
cells_cond_cutoff <- function(n, floor, alpha = NULL) {
    if (is.null(alpha) || is.na(n)) return(floor)
    if (n < 3) return(Inf)
    tc <- stats::qt(1 - alpha / 2, n - 2)
    max(floor, tc / sqrt(n - 2 + tc^2))
}
