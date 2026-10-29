#' Confounding strip: condition-covariate association per comparison
#'
#' A Layer-1 lollipop of how strongly the condition is tangled with its
#' covariates, with the Caution threshold drawn in and points coloured by the
#' call. A comparison past the line cannot cleanly separate the condition from the
#' covariate.
#'
#' @param res The result of [assessAll()] or [callComparisons()].
#' @param metric `"cramerV"` (bias-corrected Cramer's V, categorical covariates)
#'   or `"vif"` (variance inflation factor). Default `"cramerV"`.
#' @return A `ggplot`, faceted by annotation set when there is more than one.
#' @export
#' @importFrom ggplot2 .data
plotConfounding <- function(res, metric = c('cramerV', 'vif')) {
    metric <- match.arg(metric)
    df <- as.data.frame(res)
    df$call <- factor(as.character(df$call), levels = names(seny_call_colors))
    df$cluster <- factor(df$cluster, levels = rev(sort(unique(df$cluster))))
    th <- defaultThresholds()
    if (metric == 'cramerV') {
        df$val <- df$max_cramerV; thr <- th$cramerV
        xlab <- "condition-covariate association (Cramer's V)"
    } else {
        df$val <- df$max_vif; thr <- th$vif
        xlab <- 'variance inflation factor'
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$val, y = .data$cluster, color = .data$call)) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$val, yend = .data$cluster),
                              color = 'grey70', linewidth = 1) +
        ggplot2::geom_point(size = 4) +
        ggplot2::geom_vline(xintercept = thr, linetype = 'dashed', color = '#D73027') +
        ggplot2::scale_color_manual(values = seny_call_colors, drop = FALSE, name = 'call') +
        ggplot2::labs(x = xlab, y = NULL, title = 'Confounding (Layer 1)') +
        seny_theme()
    facet_by_annotation(p, df)
}
