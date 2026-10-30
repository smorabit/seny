#' Layer 4 (sensitivity): minimum detectable effect and power
#'
#' From the variance fits ([estimateVariance()]), computes the
#' analytic standard error of the two-group log2 fold-change contrast per gene
#' (`SE^2 = 1/sum_A w + 1/sum_B w`, `w = 1/(1/mu + phi)`, phi the dispersion with the
#' condition in the design), the minimum
#' detectable effect at the target power and FDR, and the power at the primary
#' target effect. Summarises to the sensitivity fraction: the share of testable
#' genes detectable at `targetLFC`. Uncalibrated in the MVP (`calibrate = FALSE`);
#' the empirical calibration from E2/E7 is applied in a later phase.
#'
#' The FDR target is mapped to a per-gene significance level through the
#' fraction of true positives `pi1` (Jung 2005), capped at `fdr` (BH never tests
#' a gene at more than the FDR level). BH's realised per-gene threshold scales
#' with the fraction of genes it calls, so no fixed `pi1` suits every contrast: a
#' sparse contrast (male vs female in blood, realised threshold ~1e-4) makes a
#' fixed 5% optimistic, a dense one makes it conservative. By default
#' (`pi1 = "discoveries"`) each comparison uses its own fraction of BH discoveries,
#' `R / m` from the DE test of the chosen method at `fdr` (floored at `1 / m`), so
#' power is evaluated at the threshold BH actually uses there. On OneK1K with 0.5-20\%
#' of genes injected as DE and on male vs female it was the best calibrated of the
#' rules tested (E4b).
#'
#' @param variance The list returned by [estimateVariance()].
#' @param targetLFC Primary target absolute log2 fold-change for the call.
#'   Default 1.
#' @param power Target power. Default 0.8.
#' @param fdr Target FDR. Default 0.05.
#' @param pi1 Fraction of true positives that maps the FDR to a per-gene
#'   threshold: `"discoveries"` (default) estimates it per comparison as the
#'   fraction of genes the comparison's own DE test calls at `fdr`; a number
#'   fixes it for every comparison (falls back to 0.05 for a comparison whose
#'   fit has no p-values).
#' @param alphaFixed Optional fixed per-gene significance level. `NULL` (default)
#'   uses the `pi1`-mapped level alone; a number makes the call use the stricter
#'   of the two (the larger MDE).
#' @param calibrate Reserved; must be `FALSE` in the MVP.
#' @return A `data.frame`, one row per comparison, with `target_lfc`, the `pi1`
#'   and per-gene `alpha` used, `sensitivity_fraction`, `median_mde`,
#'   `median_mde_fdr`, `median_mde_fixed` (`NA` without `alphaFixed`) and
#'   `median_power`.
#' @export
estimateSensitivity <- function(variance, targetLFC = 1, power = 0.8, fdr = 0.05,
                                pi1 = 'discoveries', alphaFixed = NULL, calibrate = FALSE) {
    if (isTRUE(calibrate)) stop('power calibration is not available in the MVP (calibrate must be FALSE)')
    if (is.character(pi1) && !identical(pi1, 'discoveries')) {
        stop('pi1 must be "discoveries" or a number in (0, 1)')
    }
    rows <- lapply(names(variance$pergene), function(key) {
        fit <- variance$pergene[[key]]
        cur_pi1 <- resolve_pi1(pi1, fit$pvalue, fdr)
        ap <- analytic_power(fit, targetLFC, power, fdr, cur_pi1, alphaFixed)
        data.frame(comparison = key, target_lfc = targetLFC, pi1 = cur_pi1,
                   alpha = sensitivity_alpha(power, fdr, cur_pi1, alphaFixed)[['call']],
                   sensitivity_fraction = mean(ap$mde <= targetLFC),
                   median_mde = stats::median(ap$mde),
                   median_mde_fdr = stats::median(ap$mde_fdr),
                   median_mde_fixed = stats::median(ap$mde_fixed),
                   median_power = stats::median(ap$power), stringsAsFactors = FALSE)
    })
    # comparisons with no fit (non-estimable) get an empty-sensitivity row
    seen <- vapply(rows, function(r) r$comparison, character(1))
    missing <- setdiff(variance$summary$comparison, seen)
    for (key in missing) rows[[length(rows) + 1]] <- data.frame(comparison = key,
        target_lfc = targetLFC, pi1 = NA_real_, alpha = NA_real_, sensitivity_fraction = NA_real_,
        median_mde = NA_real_, median_mde_fdr = NA_real_, median_mde_fixed = NA_real_,
        median_power = NA_real_, stringsAsFactors = FALSE)
    do.call(rbind, rows)
}

#' The comparison's pi1: a fixed number, or its fraction of BH discoveries at
#' fdr (floored at 1/m, so the mapped level stays near BH's first step)
#' @noRd
resolve_pi1 <- function(pi1, pvalue, fdr) {
    if (is.numeric(pi1)) return(pi1)
    p <- pvalue[!is.na(pvalue)]
    if (length(p) == 0) return(0.05)
    m <- length(p)
    max(sum(stats::p.adjust(p, 'BH') < fdr) / m, 1 / m)
}

#' Per-gene significance levels: the FDR mapped through pi1 (Jung 2005, capped
#' at fdr) and the optional fixed alpha. The call tests at the stricter of the two.
#' @noRd
sensitivity_alpha <- function(power, fdr, pi1, alphaFixed = NULL) {
    alpha_fdr <- min(fdr * power * pi1 / ((1 - fdr) * (1 - pi1)), fdr)
    fixed <- if (is.null(alphaFixed)) NA_real_ else alphaFixed
    c(fdr = alpha_fdr, fixed = fixed, call = min(alpha_fdr, fixed, na.rm = TRUE))
}

#' Per-gene MDE and power from the per-gene SE (method-agnostic)
#'
#' [estimateVariance()] computes `se_log2` and `df` per gene under the chosen DE
#' method (NB weights for edgeR/deseq2, precision weights for voom), so this turns
#' them into power and the minimum detectable effect the same way for every
#' method. `mde` is the larger of the FDR-based and fixed-alpha MDEs (the
#' FDR-based one without `alphaFixed`) and `power` the power at the stricter
#' threshold.
#' @noRd
analytic_power <- function(fit, targetLFC, power, fdr, pi1, alphaFixed = NULL) {
    a <- sensitivity_alpha(power, fdr, pi1, alphaFixed)
    se_log2 <- fit$se_log2
    df <- fit$df
    mde_at <- function(alpha) (stats::qt(1 - alpha / 2, df) + stats::qt(power, df)) * se_log2
    ncp <- targetLFC / se_log2
    tcrit <- stats::qt(1 - a[['call']] / 2, df)
    pwr <- stats::pt(tcrit, df, ncp, lower.tail = FALSE) +
        stats::pt(-tcrit, df, ncp, lower.tail = TRUE)
    mde_fdr <- mde_at(a[['fdr']])
    mde_fixed <- if (is.na(a[['fixed']])) rep(NA_real_, length(se_log2)) else mde_at(a[['fixed']])
    data.frame(gene = fit$genes, se_log2 = se_log2, power = pwr,
               mde = pmax(mde_fdr, mde_fixed, na.rm = TRUE), mde_fdr = mde_fdr,
               mde_fixed = mde_fixed, stringsAsFactors = FALSE)
}
