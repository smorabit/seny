#' Assemble the comparison table and make the adequacy call
#'
#' Merges the layer metrics into one tidy table, one row per comparison, and
#' assigns each comparison a call: **Invalid** (the contrast is not estimable, or
#' too few samples per group), **Caution** (estimable but the condition is
#' confounded with a covariate), **Underpowered** (valid but the sensitivity
#' fraction is below target), or **Adequate**. Severity order is
#' Invalid > Caution > Underpowered > Adequate; the `reason` column lists every
#' flag so the single call never hides a second problem.
#'
#' @param design,pseudobulk,variance,sensitivity The outputs of [assessDesign()],
#'   [assessPseudobulk()], [estimateVariance()] (its `$summary`) and
#'   [estimateSensitivity()].
#' @param thresholds A list from [defaultThresholds()].
#' @return A [S4Vectors::DataFrame], one row per comparison, with all metrics,
#'   the `call` and a `reason` string.
#' @export
#' @importFrom S4Vectors DataFrame
callComparisons <- function(design, pseudobulk, variance, sensitivity,
                            thresholds = defaultThresholds()) {
    tab <- Reduce(function(a, b) merge(a, b, by = 'comparison', all = TRUE),
                  list(design, pseudobulk, variance, sensitivity))
    th <- thresholds
    call <- character(nrow(tab)); reason <- character(nrow(tab))

    for (i in seq_len(nrow(tab))) {
        r <- tab[i, ]
        flags <- character(0)
        confounded <- (!is.na(r$max_cramerV) && r$max_cramerV >= th$cramerV) ||
                      (!is.na(r$max_vif) && r$max_vif >= th$vif)
        underpowered <- !is.na(r$sensitivity_fraction) && r$sensitivity_fraction < th$sensitivity
        too_few <- is.na(r$estimable) || !r$estimable ||
                   min(r$nA, r$nB) < th$minSamples

        if (confounded) flags <- c(flags, sprintf('confounded (Cramer V %.2f, VIF %.1f)',
                                                  r$max_cramerV, r$max_vif))
        if (underpowered) flags <- c(flags, sprintf('underpowered (%.0f%% of genes detectable at log2FC %.2g)',
                                                    100 * r$sensitivity_fraction, r$target_lfc))
        if (too_few) flags <- c(flags, if (isTRUE(r$estimable)) 'too few samples per group'
                                       else 'design not estimable (rank-deficient or nested)')

        call[i] <- if (too_few) 'Invalid'
                   else if (confounded) 'Caution'
                   else if (underpowered) 'Underpowered'
                   else 'Adequate'
        reason[i] <- if (length(flags)) paste(flags, collapse = '; ') else 'adequate'
    }
    tab$call <- factor(call, levels = c('Adequate', 'Underpowered', 'Caution', 'Invalid'))
    tab$reason <- reason
    S4Vectors::DataFrame(tab)
}
