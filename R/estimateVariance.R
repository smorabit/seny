#' Layer 4 (variance): null-model dispersion and its decomposition
#'
#' For every comparison, fits the null model (design without the condition term)
#' with edgeR to get the pseudobulk negative-binomial dispersion `phi_total`,
#' then splits it into a cell-heterogeneity term (`theta`/n, with the cell-level
#' `theta` estimated at [pbSetup()]) and a between-donor term `phi_donor` via the
#' moment estimate `phi_donor = max(0, phi_total - theta * mean(1/n))`. Reports
#' the per-gene fits needed by [estimateSensitivity()] and a per-comparison
#' summary of the variance shares (depth, cell, donor) with a regime label.
#'
#' @param x A [S4Vectors::SimpleList] from [pbSetup()].
#' @return A list with `summary` (a `data.frame`, one row per comparison) and
#'   `pergene` (a named list of per-comparison fit objects for sensitivity).
#' @export
#' @importFrom edgeR DGEList filterByExpr calcNormFactors estimateDisp glmQLFit aveLogCPM
estimateVariance <- function(x) {
    summ <- list(); pergene <- list()
    for (se in as.list(x)) {
        cmp <- enumerate_comparisons(se)
        for (i in seq_len(nrow(cmp))) {
            d <- comparison_data(se, cmp$cluster[i], cmp$A[i], cmp$B[i])
            key <- cmp$comparison[i]
            fit <- if (!is.null(d) && length(d$group) >= 2) fit_null(d) else NULL
            if (is.null(fit)) {
                summ[[key]] <- data.frame(comparison = key, median_phi_total = NA_real_,
                    median_phi_donor = NA_real_, median_theta = NA_real_,
                    share_depth = NA_real_, share_cell = NA_real_, share_donor = NA_real_,
                    regime = NA_character_, stringsAsFactors = FALSE)
                next
            }
            inv_n <- mean(1 / fit$n_i)
            theta_g <- fit$theta
            phi_donor <- pmax(0, fit$phi_total - theta_g * inv_n)
            v_depth <- 1 / fit$mean_count
            v_cell <- theta_g * inv_n
            v_donor <- phi_donor
            v_tot <- v_depth + v_cell + v_donor
            sh <- cbind(depth = v_depth / v_tot, cell = v_cell / v_tot, donor = v_donor / v_tot)
            regime <- c('depth', 'cell', 'donor')[max.col(sh, ties.method = 'first')]

            fit$phi_donor <- phi_donor
            pergene[[key]] <- fit
            summ[[key]] <- data.frame(comparison = key,
                median_phi_total = stats::median(fit$phi_total),
                median_phi_donor = stats::median(phi_donor),
                median_theta = stats::median(theta_g),
                share_depth = mean(sh[, 'depth']), share_cell = mean(sh[, 'cell']),
                share_donor = mean(sh[, 'donor']),
                regime = names(which.max(table(regime))), stringsAsFactors = FALSE)
        }
    }
    list(summary = do.call(rbind, summ), pergene = pergene)
}

#' Null-model edgeR fit for one comparison; aligns stored theta to the kept genes
#' @noRd
fit_null <- function(d) {
    design <- build_null_design(d$covariates, length(d$group))
    y <- edgeR::DGEList(d$counts)
    keep <- edgeR::filterByExpr(y, design)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    y <- edgeR::estimateDisp(y, design, robust = TRUE)
    qfit <- edgeR::glmQLFit(y, design, robust = TRUE)
    g <- rownames(y)
    theta_map <- stats::setNames(d$theta$theta, d$theta$gene)
    theta_g <- unname(theta_map[g]); theta_g[is.na(theta_g)] <- 0
    list(genes = g, mu = qfit$fitted.values, phi_total = y$tagwise.dispersion,
         df_prior = stats::median(qfit$df.prior), group = d$group,
         n_i = d$nCells, theta = theta_g,
         mean_count = rowMeans(d$counts[g, , drop = FALSE]),
         acpm = edgeR::aveLogCPM(y), resid_df = ncol(design) )
}

#' Null design: covariates only, no condition term
#' @noRd
build_null_design <- function(covariates, n) {
    if (is.null(covariates) || ncol(covariates) == 0) matrix(1, n, 1)
    else stats::model.matrix(~ ., data = as.data.frame(covariates))
}
