#' Default decision thresholds for the adequacy call
#'
#' Provisional values used by [callComparisons()]. The final numbers are to be
#' derived from the validation experiments (E1-E8); until then these are sensible
#' placeholders, documented so a user can override them.
#'
#' @return A named list: `cramerV` (condition-covariate association above which a
#'   comparison is flagged Caution), `vif` (variance inflation factor threshold),
#'   `sensitivity` (sensitivity fraction at or above which a comparison is
#'   adequately powered), and `minSamples` (fewer samples per group is Invalid).
#' @export
defaultThresholds <- function() {
    list(cramerV = 0.3, vif = 5, sensitivity = 0.5, minSamples = 2)
}
