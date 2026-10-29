#' Layer 1: design validity and confounding
#'
#' For every comparison, measures whether the two-group contrast is estimable
#' (full design rank with residual degrees of freedom) and how strongly the
#' condition is confounded with each candidate covariate: the variance inflation
#' factor, bias-corrected Cramer's V for categorical covariates, and the
#' standardised mean difference for numeric ones. A rank-deficient or fully
#' nested design is flagged non-estimable, which downstream becomes an Invalid
#' call.
#'
#' @param x A [S4Vectors::SimpleList] from [pbSetup()].
#' @return A `data.frame`, one row per comparison, with design and confounding
#'   metrics (`estimable`, `resid_df`, `nA`, `nB`, `max_vif`, `max_cramerV`,
#'   `max_smd`).
#' @export
assessDesign <- function(x) {
    do.call(rbind, lapply(x, function(se) {
        cmp <- enumerate_comparisons(se)
        do.call(rbind, lapply(seq_len(nrow(cmp)), function(i) {
            d <- comparison_data(se, cmp$cluster[i], cmp$A[i], cmp$B[i])
            row <- cmp[i, ]
            if (is.null(d) || length(d$group) < 2) {
                return(cbind(row, nA = 0L, nB = 0L, estimable = FALSE,
                             resid_df = NA_integer_, max_vif = NA_real_,
                             max_cramerV = NA_real_, max_smd = NA_real_))
            }
            cond <- as.integer(d$group == 'B')
            bd <- build_design(d$group, d$covariates)
            vif <- cram <- smds <- numeric(0)
            if (!is.null(d$covariates)) {
                for (cov in d$covariates) {
                    vif <- c(vif, vif_cov(cond, cov))
                    if (is.numeric(cov)) smds <- c(smds, abs(smd(cov[cond == 0], cov[cond == 1])))
                    else cram <- c(cram, cramers_v(d$group, cov))
                }
            }
            cbind(row, nA = sum(d$group == 'A'), nB = sum(d$group == 'B'),
                  estimable = bd$estimable, resid_df = bd$resid_df,
                  max_vif = if (length(vif)) max(vif) else 1,
                  max_cramerV = if (length(cram)) max(cram) else 0,
                  max_smd = if (length(smds)) max(smds) else 0)
        }))
    }))
}
