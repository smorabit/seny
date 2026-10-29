#---------------------------------------------------------#
# Cell-level overdispersion theta_g (internal)
#
# theta_g is the within-donor cell-to-cell NB overdispersion. It is estimated
# with a cell-level glmGamPoi fit using donor as the only covariate, so the
# overdispersion is the residual cell heterogeneity after removing donor means.
# glmGamPoi has no fast path for a many-donor design and is memory-bandwidth
# bound, so a narrow subsample is fit serially and averaged over a few draws
# (theta_g is shared across donors, so the subsample estimates it just as well).
#---------------------------------------------------------#

#' Draw one small (donor-subsampled, cell-capped) matrix + donor factor
#' @noRd
draw_theta_sample <- function(counts, donor, max_cells = 50, max_donors = 50,
                              min_cells_gene = 20, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    donor <- as.character(donor)
    ud <- unique(donor)
    if (length(ud) > max_donors) {
        sel <- donor %in% sample(ud, max_donors)
        counts <- counts[, sel]
        donor <- donor[sel]
    }
    sub <- cap_cells(donor, max_cells)
    mat <- counts[, sub]
    keep <- Matrix::rowSums(mat > 0) >= min_cells_gene
    list(mat = mat[keep, ], d = droplevels(factor(donor[sub])))
}

#' Fit theta_g (unshrunk overdispersion) on one prepared sample
#' @noRd
fit_theta <- function(sample) {
    mat <- sample$mat
    sf <- normed_sf(mat)
    fit <- glmGamPoi::glm_gp(mat, design = stats::model.matrix(~ sample$d),
                             size_factors = sf, overdispersion = TRUE,
                             overdispersion_shrinkage = FALSE)
    stats::setNames(fit$overdispersions, rownames(mat))
}

#' Ensemble theta_g over n_rep random (donor, cell) draws
#'
#' Returns a data.frame(gene, theta, theta_sd). Fits run serially: a single
#' glm_gp fit is memory-bandwidth bound, so running several at once just
#' saturates bandwidth and each one crawls.
#' @noRd
estimate_cell_theta <- function(counts, donor, n_rep = 3, max_cells = 50,
                                max_donors = 50, min_cells_gene = 20, seed = 1) {
    reps <- lapply(seq_len(n_rep), function(r) {
        s <- draw_theta_sample(counts, donor, max_cells, max_donors,
                               min_cells_gene, seed = seed + r)
        fit_theta(s)
    })
    gg <- Reduce(intersect, lapply(reps, names))
    tmat <- vapply(reps, function(v) v[gg], numeric(length(gg)))
    data.frame(gene = gg, theta = rowMeans(tmat),
               theta_sd = if (n_rep > 1) apply(tmat, 1, stats::sd) else NA_real_,
               stringsAsFactors = FALSE)
}
