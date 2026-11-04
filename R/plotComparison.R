#' Comparison card: a per-comparison deep dive
#'
#' Assembles one comparison into a single card: the call and its reason, the
#' variance attribution as a donut, the sensitivity curve, and cells per sample.
#' This is the per-comparison unit behind the eventual report.
#'
#' @param res The result of [assessAll()] (its `pb` and per-gene `pergene` are
#'   required).
#' @param comparison A comparison key from `res$comparison`.
#' @return A `patchwork` composition of `ggplot`s.
#' @examples
#' sce <- readRDS(system.file("extdata", "example_sce.rds", package = "seny"))
#' res <- assessAll(sce, sample = "donor", group = "condition",
#'                  annotations = "cell_type", targetLFC = 1, thetaReps = 0)
#' plotComparison(res, "cell_type:CD4T:B_vs_A")
#' @export
#' @importFrom ggplot2 .data
plotComparison <- function(res, comparison) {
    df <- as.data.frame(res)
    if (!comparison %in% df$comparison) {
        stop('unknown comparison: "', comparison, '"; see res$comparison', call. = FALSE)
    }
    row <- df[df$comparison == comparison, ]
    hdr <- ggplot2::ggplot() + ggplot2::theme_void() +
        ggplot2::annotate('text', x = 0, y = 1, hjust = 0, size = 5.3, fontface = 'bold',
            label = paste0(row$cluster, ' : ', row$B, ' vs ', row$A, '  \u2014  ', as.character(row$call))) +
        ggplot2::annotate('text', x = 0, y = 0.55, hjust = 0, size = 3.5, colour = 'grey30',
            label = paste0('reason: ', row$reason)) +
        ggplot2::annotate('text', x = 0, y = 0.15, hjust = 0, size = 3.2, colour = 'grey40',
            label = composition_text(row)) +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1.3)
    # without theta (thetaReps = 0) there is no decomposition, so no donut
    parts <- list(sens_one(res, comparison), cells_one(res, row))
    if (!is.na(row$share_donor)) parts <- c(list(comparison_donut(row)), parts)
    panels <- patchwork::wrap_plots(parts, nrow = 1)
    patchwork::wrap_plots(hdr, panels, ncol = 1, heights = c(0.28, 1))
}

#' Variance donut for one comparison
#' @noRd
comparison_donut <- function(row) {
    sh <- data.frame(term = factor(c('depth', 'cell', 'donor'), levels = c('depth', 'cell', 'donor')),
                     share = c(row$share_depth, row$share_cell, row$share_donor))
    dom <- sh[which.max(sh$share), ]
    center <- paste0(dom$term, '\n', round(100 * dom$share), '%')
    ggplot2::ggplot(sh, ggplot2::aes(x = 2, y = .data$share, fill = .data$term)) +
        ggplot2::geom_col(width = 1, color = 'white', linewidth = 0.6) +
        ggplot2::coord_polar(theta = 'y') + ggplot2::xlim(0, 2.5) +
        ggplot2::annotate('text', x = 0, y = 0, label = center, size = 4.2, fontface = 'bold', lineheight = 0.9) +
        ggplot2::scale_fill_manual(values = seny_share_colors, name = NULL) +
        ggplot2::labs(title = 'attribution') +
        ggplot2::theme_void(base_size = 12) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = 'bold', hjust = 0.5),
                       legend.position = 'bottom')
}

#' Per-gene MDE vs expression for one comparison
#' @noRd
sens_one <- function(res, comparison) {
    pg <- require_attr(res, 'pergene')
    pg <- pg[pg$comparison == comparison, , drop = FALSE]
    df <- as.data.frame(res)
    target <- df$target_lfc[df$comparison == comparison][1]
    pg$mde <- mde_from_se(pg$se_log2, pg$df, df$alpha[df$comparison == comparison][1])
    ymax <- stats::quantile(pg$mde, 0.98, na.rm = TRUE)
    ggplot2::ggplot(pg, ggplot2::aes(x = .data$acpm, y = .data$mde)) +
        ggplot2::geom_point(size = 0.5, alpha = 0.3, colour = '#4C72B0') +
        ggplot2::geom_hline(yintercept = target, linetype = 'dashed', colour = '#D73027') +
        ggplot2::coord_cartesian(ylim = c(0, ymax)) +
        ggplot2::labs(x = 'aveLogCPM', y = 'MDE', title = 'sensitivity') + seny_theme()
}

#' Cells per sample for one comparison's cluster and groups
#' @noRd
cells_one <- function(res, row) {
    pb <- require_attr(res, 'pb')
    se <- pb[[row$annotation]]
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    cd <- cd[cd$cluster == row$cluster & cd$group %in% c(row$A, row$B), , drop = FALSE]
    cd$group <- factor(cd$group, levels = c(row$A, row$B))
    ggplot2::ggplot(cd, ggplot2::aes(x = .data$group, y = .data$nCells, color = .data$group)) +
        ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5) +
        ggbeeswarm::geom_quasirandom(width = 0.25, size = 1.2, alpha = 0.6) +
        ggplot2::scale_color_manual(values = seny_cond_colors, guide = 'none') +
        ggplot2::labs(x = 'condition', y = 'cells', title = 'cells per sample') + seny_theme()
}

#' One line on the comparison's composition for the card header
#' @noRd
composition_text <- function(row) {
    if (is.null(row$composition)) return('')
    if (row$composition != 'assessed') {
        return(paste0('composition: ', row$composition,
                      if (row$composition == 'not assessed') ' (needs a finer annotation nested in this one)' else ''))
    }
    sprintf('composition: most shifted subtype %s %.0f%% -> %.0f%% (FDR %.2g); ~%.0f genes expected to look DE from the shift alone; %d of %d parent DEGs composition-driven',
            row$comp_child, 100 * row$comp_share_A, 100 * row$comp_share_B, row$comp_min_fdr,
            row$exp_comp_deg, as.integer(row$n_comp_driven), as.integer(row$n_parent_deg))
}
