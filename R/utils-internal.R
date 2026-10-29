#---------------------------------------------------------#
# Internal helpers (not exported)
#---------------------------------------------------------#

#' Sum a genes x cells matrix into genes x groups
#' @param counts genes x cells (sparse or dense)
#' @param grp grouping factor aligned to columns
#' @return dense genes x groups matrix, columns named by group level
#' @noRd
pb_sum <- function(counts, grp) {
    g <- factor(as.character(grp))
    design <- Matrix::sparse.model.matrix(~ 0 + g)
    pb <- counts %*% design
    colnames(pb) <- levels(g)
    as.matrix(pb)
}

#' Subsample at most max_cells column indices per group
#' @noRd
cap_cells <- function(grp, max_cells) {
    idx <- split(seq_along(grp), as.character(grp))
    unlist(lapply(idx, function(ii) {
        if (length(ii) > max_cells) sample(ii, max_cells) else ii
    }), use.names = FALSE)
}

#' 'normed_sum' size factors: column totals scaled to geometric mean 1
#' @noRd
normed_sf <- function(mat) {
    tot <- Matrix::colSums(mat)
    tot / exp(mean(log(tot[tot > 0])))
}

#' Gini coefficient of a non-negative vector (0 = equal, 1 = maximally unequal)
#' @noRd
gini <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 2 || sum(x) == 0) return(0)
    x <- sort(x)
    2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1) / n
}

#' Standardised mean difference between two numeric vectors
#' @noRd
smd <- function(x, y) {
    s <- sqrt((stats::var(x) + stats::var(y)) / 2)
    if (is.na(s) || s == 0) return(0)
    (mean(x) - mean(y)) / s
}

#' Bias-corrected Cramer's V (Bergsma 2013) between two categorical vectors
#' @noRd
cramers_v <- function(a, b) {
    tab <- table(a, b)
    tab <- tab[, colSums(tab) > 0, drop = FALSE]
    tab <- tab[rowSums(tab) > 0, , drop = FALSE]
    n <- sum(tab)
    if (n == 0 || nrow(tab) < 2 || ncol(tab) < 2) return(0)
    chi2 <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$statistic)
    phi2 <- as.numeric(chi2) / n
    r <- nrow(tab); k <- ncol(tab)
    phi2_c <- max(0, phi2 - (k - 1) * (r - 1) / (n - 1))
    r_c <- r - (r - 1)^2 / (n - 1)
    k_c <- k - (k - 1)^2 / (n - 1)
    denom <- min(k_c - 1, r_c - 1)
    if (denom <= 0) return(NA_real_)
    sqrt(phi2_c / denom)
}

#' Variance inflation factor of the condition indicator against one covariate
#' @noRd
vif_cov <- function(cond, cov) {
    if (!is.numeric(cov)) {
        cov <- droplevels(factor(cov))
        if (nlevels(cov) < 2) return(1)
    } else if (stats::var(cov) == 0) {
        return(1)
    }
    r2 <- suppressWarnings(summary(stats::lm(cond ~ cov))$r.squared)
    if (r2 >= 1 - 1e-9) Inf else 1 / (1 - r2)
}

#' Build the two-group design and report whether the contrast is estimable
#' @noRd
build_design <- function(group, covariates = NULL) {
    if (is.null(covariates) || ncol(covariates) == 0) {
        design <- stats::model.matrix(~ group)
    } else {
        design <- stats::model.matrix(~ ., data = cbind(covariates, group = group))
    }
    coef_idx <- grep('group', colnames(design))
    coef_idx <- coef_idx[length(coef_idx)]
    rank_ok <- qr(design)$rank == ncol(design)
    resid_df <- nrow(design) - qr(design)$rank
    list(design = design, coef = coef_idx, estimable = rank_ok && resid_df >= 1,
         resid_df = resid_df)
}
