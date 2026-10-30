#' Layer 4 (variance): null-model dispersion, decomposition and per-gene SE
#'
#' For every comparison, fits the chosen DE method and returns the per-gene
#' standard error of the two-group log2 fold-change contrast (used by
#' [estimateSensitivity()]). The SE uses the dispersion estimated with the
#' condition in the design, as the user's test sees it, evaluated at the
#' condition-blind (null-model) fitted means. Power is always computed at the
#' target effect, never the observed one. A null-model dispersion would absorb
#' real effects as noise and understate power exactly where DE is strong. The
#' variance decomposition stays condition-blind. For the count methods (`edgeR`, `deseq2`) it also
#' splits the pseudobulk negative-binomial dispersion `phi_total` into a
#' cell-heterogeneity term (`theta`/n, with the cell-level `theta` from [pbSetup()])
#' and a between-donor term `phi_donor` (moment estimate
#' `max(0, phi_total - theta*mean(1/n))`), and reports the depth/cell/donor variance
#' shares with a regime label. The `voom` method uses limma's precision-weight
#' variance model, so it has no NB decomposition (those columns are `NA`).
#'
#' The method should match the differential-expression method the user intends to
#' run. `edgeR` and `deseq2` require sum aggregation (integer counts); `voom`
#' works on either sum (voom) or mean (limma-trend) pseudobulks. `dreamlet` is
#' reserved for mixed designs and is not yet implemented.
#'
#' @param x A [S4Vectors::SimpleList] from [pbSetup()].
#' @param method One of `"edgeR"`, `"deseq2"`, `"voom"` (or the reserved
#'   `"dreamlet"`). Default `"edgeR"`.
#' @return A list with `summary` (a `data.frame`, one row per comparison, with the
#'   `method`, the NB decomposition and regime) and `pergene` (a named list of
#'   per-comparison fits carrying `se_log2` and `df` for sensitivity, and the
#'   method's own p-values for the condition, `pvalue`, behind the
#'   `pi1 = "discoveries"` mapping in [estimateSensitivity()]).
#' @export
#' @importFrom edgeR DGEList filterByExpr calcNormFactors estimateDisp glmQLFit aveLogCPM cpm
estimateVariance <- function(x, method = c('edgeR', 'deseq2', 'voom', 'dreamlet')) {
    method <- match.arg(method)
    summ <- list(); pergene <- list()
    for (se in as.list(x)) {
        aggregation <- S4Vectors::metadata(se)$aggregation
        check_method_aggregation(method, aggregation)
        cmp <- enumerate_comparisons(se)
        for (i in seq_len(nrow(cmp))) {
            d <- comparison_data(se, cmp$cluster[i], cmp$A[i], cmp$B[i])
            key <- cmp$comparison[i]
            fit <- if (!is.null(d) && length(d$group) >= 2) fit_null(d, method, aggregation) else NULL
            if (is.null(fit)) {
                summ[[key]] <- variance_summary(key, method, NULL)
                next
            }
            pergene[[key]] <- fit
            summ[[key]] <- variance_summary(key, method, fit)
        }
    }
    list(summary = do.call(rbind, summ), pergene = pergene)
}

#' Guard: count methods need integer sums; mean pseudobulks need limma (voom)
#' @noRd
check_method_aggregation <- function(method, aggregation) {
    if (method %in% c('edgeR', 'deseq2') && identical(aggregation, 'mean')) {
        stop('method "', method, '" needs sum aggregation (integer counts); ',
             'use aggregation = "sum" or method = "voom"')
    }
    if (identical(aggregation, 'mean') && method != 'voom') {
        stop('mean aggregation pairs with limma; use method = "voom"')
    }
}

#' Per-comparison summary row: NB decomposition for count methods, NA for voom
#' @noRd
variance_summary <- function(key, method, fit) {
    na_row <- data.frame(comparison = key, method = method,
        median_phi_total = NA_real_, median_phi_donor = NA_real_, median_theta = NA_real_,
        share_depth = NA_real_, share_cell = NA_real_, share_donor = NA_real_,
        regime = NA_character_, stringsAsFactors = FALSE)
    if (is.null(fit) || !isTRUE(fit$nb)) return(na_row)
    # theta skipped (thetaReps = 0): total dispersion only, no decomposition
    if (all(is.na(fit$theta))) {
        na_row$median_phi_total <- stats::median(fit$phi_total)
        return(na_row)
    }

    inv_n <- mean(1 / fit$n_i)
    phi_donor <- pmax(0, fit$phi_total - fit$theta * inv_n)
    v_depth <- 1 / fit$mean_count
    v_cell <- fit$theta * inv_n
    v_donor <- phi_donor
    v_tot <- v_depth + v_cell + v_donor
    sh <- cbind(depth = v_depth / v_tot, cell = v_cell / v_tot, donor = v_donor / v_tot)
    regime <- c('depth', 'cell', 'donor')[max.col(sh, ties.method = 'first')]
    data.frame(comparison = key, method = method,
        median_phi_total = stats::median(fit$phi_total),
        median_phi_donor = stats::median(phi_donor),
        median_theta = stats::median(fit$theta),
        share_depth = mean(sh[, 'depth']), share_cell = mean(sh[, 'cell']),
        share_donor = mean(sh[, 'donor']),
        regime = names(which.max(table(regime))), stringsAsFactors = FALSE)
}

#' Dispatch the null-model fit to the chosen method
#' @noRd
fit_null <- function(d, method, aggregation) {
    switch(method,
        edgeR = fit_edger(d),
        deseq2 = fit_deseq2(d),
        voom = fit_voom(d, aggregation),
        dreamlet = stop('method "dreamlet" is post-MVP (mixed / repeated-measures designs) ',
                        'and not yet implemented'))
}

#' SE of the two-group log2FC contrast from an NB mean/dispersion fit
#' @noRd
nb_se <- function(mu, phi, group) {
    is_a <- group == 'A'
    w <- 1 / (1 / mu + phi)
    sum_a <- rowSums(w[, is_a, drop = FALSE])
    sum_b <- rowSums(w[, !is_a, drop = FALSE])
    sqrt(1 / sum_a + 1 / sum_b) / log(2)
}

#' Align the stored cell-level theta to a set of genes (missing -> 0; all NA
#' when theta was skipped with thetaReps = 0)
#' @noRd
align_theta <- function(theta_df, genes) {
    if (is.null(theta_df)) return(rep(NA_real_, length(genes)))
    m <- stats::setNames(theta_df$theta, theta_df$gene)
    t <- unname(m[genes]); t[is.na(t)] <- 0
    t
}

#' edgeR-QL null fit: tagwise NB dispersion, fitted means, QL prior df
#' @noRd
fit_edger <- function(d) {
    design <- build_null_design(d$covariates, length(d$group))
    y <- edgeR::DGEList(d$counts)
    # filter by condition group, as the user's DE (~ group) and layer 2's n_testable
    # do; the null design alone demands expression in nearly every sample and drops
    # genes expressed mainly in one condition, the strongest DE genes
    keep <- edgeR::filterByExpr(y, group = d$group)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    # null fit (condition-blind): the fitted means for the SE and phi_total for the
    # decomposition diagnostics
    y_null <- edgeR::estimateDisp(y, design, robust = TRUE)
    qfit <- edgeR::glmQLFit(y_null, design, robust = TRUE)
    g <- rownames(y)
    phi <- y_null$tagwise.dispersion
    mu <- qfit$fitted.values
    # full fit (condition in the design): the dispersion and QL prior df behind the
    # SE, as the user's test sees them. under the null fit a large true effect is
    # absorbed as noise and inflates the SE of exactly the strongly DE genes
    # (E4: ceiling genes predicted 0.22 vs detected 0.77 at 3 donors per group).
    # power stays evaluated at the target effect, never the observed one.
    full <- disp_full(y, d)
    phi_se <- if (is.null(full)) phi else full$phi
    df_prior <- if (is.null(full)) stats::median(qfit$df.prior) else full$df_prior
    df <- (length(d$group) - 2) + df_prior
    list(genes = g, se_log2 = nb_se(mu, phi_se, d$group), df = rep(df, length(g)),
         pvalue = full$pvalue,
         group = d$group, n_i = d$nCells, nb = TRUE, phi_total = phi, phi_se = phi_se,
         theta = align_theta(d$theta, g), mean_count = rowMeans(d$counts[g, , drop = FALSE]),
         acpm = edgeR::aveLogCPM(y))
}

#' Tagwise dispersion, QL prior df and the condition's QL-test p-values with the
#' condition in the design; NULL when the full design is not estimable (the
#' caller falls back to the null fit)
#' @noRd
disp_full <- function(y, d) {
    bd <- build_design(d$group, d$covariates)
    if (!bd$estimable) return(NULL)
    tryCatch({
        yf <- edgeR::estimateDisp(y, bd$design, robust = TRUE)
        qf <- edgeR::glmQLFit(yf, bd$design, robust = TRUE)
        qt <- edgeR::glmQLFTest(qf, coef = bd$coef)
        list(phi = yf$tagwise.dispersion, df_prior = stats::median(qf$df.prior),
             pvalue = stats::setNames(qt$table$PValue, rownames(qt)))
    }, error = function(e) NULL)
}

#' DESeq2 null fit: MAP-shrunk NB dispersion and fitted means (residual df)
#' @noRd
#' @importFrom stats as.formula
fit_deseq2 <- function(d) {
    counts <- round(as.matrix(d$counts)); mode(counts) <- 'integer'
    cd <- data.frame(group = d$group, row.names = colnames(counts))
    rhs <- 'group'
    if (!is.null(d$covariates) && ncol(d$covariates) > 0) {
        cd <- cbind(cd, d$covariates)
        rhs <- paste(c(colnames(d$covariates), 'group'), collapse = ' + ')
    }
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, cd, stats::as.formula(paste('~', rhs)))
    dds <- DESeq2::estimateSizeFactors(dds, type = 'poscounts')
    dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
    wald <- DESeq2::results(dds, name = utils::tail(DESeq2::resultsNames(dds), 1))
    phi <- DESeq2::dispersions(dds)
    mu <- SummarizedExperiment::assays(dds)[['mu']]
    ok <- is.finite(phi) & rowSums(is.na(mu)) == 0
    g <- rownames(dds)[ok]
    phi <- phi[ok]; mu <- mu[ok, , drop = FALSE]
    n_cov <- if (is.null(d$covariates)) 0 else ncol(d$covariates)
    df <- length(d$group) - 2 - n_cov
    list(genes = g, se_log2 = nb_se(mu, phi, d$group), df = rep(df, length(g)),
         pvalue = stats::setNames(wald$pvalue[ok], g),
         group = d$group, n_i = d$nCells, nb = TRUE, phi_total = phi,
         theta = align_theta(d$theta, g), mean_count = rowMeans(as.matrix(d$counts)[g, , drop = FALSE]),
         acpm = edgeR::aveLogCPM(edgeR::DGEList(as.matrix(d$counts)[g, , drop = FALSE])))
}

#' voom / limma-trend fit: precision-weight SE of the group coefficient
#' @noRd
#' @importFrom limma voom lmFit eBayes
fit_voom <- function(d, aggregation) {
    bd <- build_design(d$group, d$covariates)
    design <- bd$design
    coef <- bd$coef
    if (identical(aggregation, 'mean')) {
        lib <- colSums(d$counts)
        keep <- edgeR::filterByExpr(d$counts, design, lib.size = lib)
        logcpm <- edgeR::cpm(d$counts[keep, , drop = FALSE], lib.size = lib,
                             log = TRUE, prior.count = 1)
        fit <- limma::eBayes(limma::lmFit(logcpm, design), trend = TRUE)
    } else {
        y <- edgeR::DGEList(d$counts)
        keep <- edgeR::filterByExpr(y, design)
        y <- y[keep, , keep.lib.sizes = FALSE]
        y <- edgeR::calcNormFactors(y)
        v <- limma::voom(y, design)
        fit <- limma::eBayes(limma::lmFit(v, design))
    }
    g <- rownames(fit$coefficients)
    se_log2 <- fit$stdev.unscaled[, coef] * sqrt(fit$s2.post)
    list(genes = g, se_log2 = se_log2, df = fit$df.total, pvalue = fit$p.value[, coef],
         group = d$group, n_i = d$nCells, nb = FALSE, acpm = fit$Amean)
}

#' Null design: covariates only, no condition term
#' @noRd
build_null_design <- function(covariates, n) {
    if (is.null(covariates) || ncol(covariates) == 0) matrix(1, n, 1)
    else stats::model.matrix(~ ., data = as.data.frame(covariates))
}
