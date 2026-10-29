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
#'   above which cell-count imbalance is flagged Caution), and `minTestable`
#'   (fewer testable genes is Underpowered). The last two gate the layer-2
#'   metrics and, like the others, are provisional placeholders until E1-E8.
#' @export
defaultThresholds <- function() {
    list(cramerV = 0.3, vif = 5, sensitivity = 0.5, minSamples = 2,
         cellsCondCor = 0.5, minTestable = 200)
}
