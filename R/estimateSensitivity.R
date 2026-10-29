#' Layer 4 (sensitivity): minimum detectable effect and power
#'
#' From the null-model variance fits ([estimateVariance()]), computes the
#' analytic standard error of the two-group log2 fold-change contrast per gene
#' (`SE^2 = 1/sum_A w + 1/sum_B w`, `w = 1/(1/mu + phi_total)`), the minimum
#' detectable effect at the target power and FDR, and the power at the primary
#' target effect. Summarises to the sensitivity fraction: the share of testable
#' genes detectable at `targetLFC`. Uncalibrated in the MVP (`calibrate = FALSE`);
#' the empirical calibration from E2/E7 is applied in a later phase.
#'
#' @param variance The list returned by [estimateVariance()].
#' @param targetLFC Primary target absolute log2 fold-change for the call.
#'   Default 1.
#' @param power Target power. Default 0.8.
#' @param fdr Target FDR. Default 0.05.
#' @param pi1 Assumed fraction of true positives (maps FDR to a per-gene
#'   threshold, Jung 2005). Default 0.05.
#' @param calibrate Reserved; must be `FALSE` in the MVP.
#' @return A `data.frame`, one row per comparison, with `target_lfc`,
#'   `sensitivity_fraction`, `median_mde` and `median_power`.
#' @export
estimateSensitivity <- function(variance, targetLFC = 1, power = 0.8, fdr = 0.05,
                                pi1 = 0.05, calibrate = FALSE) {
    if (isTRUE(calibrate)) stop('power calibration is not available in the MVP (calibrate must be FALSE)')
    rows <- lapply(names(variance$pergene), function(key) {
        ap <- analytic_power(variance$pergene[[key]], targetLFC, power, fdr, pi1)
        data.frame(comparison = key, target_lfc = targetLFC,
                   sensitivity_fraction = mean(ap$mde <= targetLFC),
                   median_mde = stats::median(ap$mde),
                   median_power = stats::median(ap$power), stringsAsFactors = FALSE)
    })
    # comparisons with no fit (non-estimable) get an empty-sensitivity row
    seen <- vapply(rows, function(r) r$comparison, character(1))
    missing <- setdiff(variance$summary$comparison, seen)
    for (key in missing) rows[[length(rows) + 1]] <- data.frame(comparison = key,
        target_lfc = targetLFC, sensitivity_fraction = NA_real_,
        median_mde = NA_real_, median_power = NA_real_, stringsAsFactors = FALSE)
    do.call(rbind, rows)
}

#' Per-gene MDE and power from the per-gene SE (method-agnostic)
#'
#' [estimateVariance()] computes `se_log2` and `df` per gene under the chosen DE
#' method (NB weights for edgeR/deseq2, precision weights for voom), so this turns
#' them into power and the minimum detectable effect the same way for every method.
#' @noRd
analytic_power <- function(fit, targetLFC, power, fdr, pi1) {
    se_log2 <- fit$se_log2
    df <- fit$df
    beta <- 1 - power
    alpha_star <- fdr * (1 - beta) * pi1 / ((1 - fdr) * (1 - pi1))
    ncp <- targetLFC / se_log2
    tcrit <- stats::qt(1 - alpha_star / 2, df)
    pwr <- stats::pt(tcrit, df, ncp, lower.tail = FALSE) +
        stats::pt(-tcrit, df, ncp, lower.tail = TRUE)
    mde <- (stats::qt(1 - alpha_star / 2, df) + stats::qt(power, df)) * se_log2
    data.frame(gene = fit$genes, se_log2 = se_log2, power = pwr, mde = mde,
               stringsAsFactors = FALSE)
}
