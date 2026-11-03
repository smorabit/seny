#' Composition: cluster shares, subtype shifts and composition-driven DE
#'
#' A broad cluster can hold subtypes whose shares shift between the groups.
#' Pseudobulk DE on the broad cluster then reports the subtypes' marker genes as
#' differentially expressed although no gene changes within any subtype. This
#' step measures that risk for every comparison, in four parts:
#'
#' 1. **Cluster share** (every annotation level, reported only): the cluster's
#'    share of each sample's cells, A vs B, tested with `speckle::propeller`
#'    (arcsin-square-root proportions, limma). A cluster whose share shifts is
#'    not a problem for DE within it, so this never changes the call.
#' 2. **Subtype shift within the cluster** (needs a finer annotation nested in
#'    this one, see [pbSetup()]): propeller on the children's shares of the
#'    parent's cells, over the comparison's samples.
#' 3. **Composition-sensitive genes** (prospective): per gene, the log2 fold
#'    change the observed subtype shift alone would produce in the parent's
#'    library-size-normalised pseudobulk,
#'    `log2[(sum_k pB_k c_gk / sum_k pB_k L_k) / (sum_k pA_k c_gk / sum_k pA_k L_k)]`,
#'    with `p` the donor-averaged child shares, `c_gk` the per-cell mean count of
#'    gene g in child k and `L_k` the child's per-cell library size. Each gene's
#'    power to be called DE at that shift-only change (layer 4, at the
#'    comparison's per-gene alpha) is `comp_power`; their sum, `exp_comp_deg`, is
#'    the number of DEGs the shift alone is expected to produce, and
#'    `frac_comp_expected` its share of the testable genes. [callComparisons()]
#'    flags the comparison **Caution** when a child's shift is significant and
#'    `frac_comp_expected` reaches `thresholds$compFrac`. Genes whose shift-only
#'    change exceeds their own MDE (`sensitive`, counted in `n_comp_sensitive`)
#'    are reported too, but many genes each at modest power add up to many
#'    artefacts without any of them clearing the MDE, so that count understates
#'    the risk (E6).
#' 4. **Attribution of the observed parent DEGs** (retrospective, reported
#'    only): each gene the method's own test calls DE in the parent at `fdr` is
#'    labelled *composition-driven* (DE in no child comparison, and the
#'    shift-only change has the same sign and explains at least
#'    `thresholds$compExplained` of the observed change), *within-child DE* (DE
#'    in at least one child), *mixed* (both) or *unresolved* (neither, or no
#'    child comparison could be fitted).
#'
#' Clusters without a finer annotation nested in them get
#' `composition = "not assessed"`; parts 2-4 need two annotation levels.
#'
#' @param x A [S4Vectors::SimpleList] from [pbSetup()].
#' @param variance The list returned by [estimateVariance()].
#' @param sensitivity The table returned by [estimateSensitivity()] (its
#'   per-gene `alpha` sets the MDE).
#' @param power Target power for the MDE. Default 0.8.
#' @param fdr FDR at which parent and child genes count as DE for the
#'   attribution. Default 0.05.
#' @param reference Where the per-cell child expression `c_gk` comes from:
#'   `"reference"` (the contrast's reference group A only; default, which avoids
#'   attributing real within-child DE in B to composition) or `"pooled"` (all
#'   samples of the contrast's groups).
#' @param thresholds A list from [defaultThresholds()]; uses `compExplained`.
#' @return A list with `summary` (a `data.frame`, one row per comparison: the
#'   cluster share `share_A`, `share_B`, `share_fdr`; the `composition` status;
#'   the most shifted child `comp_child` with `comp_share_A`, `comp_share_B` and
#'   `comp_min_fdr`; `exp_comp_deg`, `frac_comp_expected`; `n_comp_sensitive`,
#'   `frac_comp_sensitive`; `n_parent_deg`,
#'   `n_comp_driven`, `frac_comp_driven`), `genes` (per gene and comparison:
#'   `lfc_comp`, `comp_power`, `mde`, `sensitive`, the observed `logfc`, `padj`, `parent_de`,
#'   `child_de` and the attribution `label`) and `proportions` (per sample and
#'   child, the child's share of the parent).
#' @export
assessComposition <- function(x, variance, sensitivity, power = 0.8, fdr = 0.05,
                              reference = c('reference', 'pooled'),
                              thresholds = defaultThresholds()) {
    reference <- match.arg(reference)
    th <- utils::modifyList(defaultThresholds(), thresholds)
    hier <- S4Vectors::metadata(x)$hierarchy
    sens <- as.data.frame(sensitivity)
    alpha_of <- stats::setNames(sens$alpha, sens$comparison)
    rows <- list(); genes <- list(); props <- list()

    for (se in as.list(x)) {
        md <- S4Vectors::metadata(se)
        ann <- md$annotation
        cmp <- enumerate_comparisons(se)
        child_ann <- immediate_child(hier, ann)
        share_cache <- list()
        for (i in seq_len(nrow(cmp))) {
            key <- cmp$comparison[i]; cl <- cmp$cluster[i]; A <- cmp$A[i]; B <- cmp$B[i]
            ck <- paste(A, B)
            if (is.null(share_cache[[ck]])) share_cache[[ck]] <- cluster_share(md$cell_counts, A, B)
            sh <- share_cache[[ck]]
            row <- empty_composition_row(key)
            if (!is.null(sh) && cl %in% sh$cluster) {
                j <- match(cl, sh$cluster)
                row$share_A <- sh$share_A[j]; row$share_B <- sh$share_B[j]; row$share_fdr <- sh$fdr[j]
            }
            fit <- variance$pergene[[key]]
            if (is.null(child_ann)) {
                rows[[key]] <- row; next
            }
            if (is.null(fit)) {
                row$composition <- 'no fit'; rows[[key]] <- row; next
            }
            res <- parent_composition(se, x[[child_ann]], hier, ann, child_ann, cl, A, B,
                                      key, fit, alpha_of[[key]], variance$pergene, power,
                                      fdr, reference, th)
            row$composition <- res$status
            if (!is.null(res$summary)) row[names(res$summary)] <- res$summary
            rows[[key]] <- row
            genes[[key]] <- res$genes
            props[[key]] <- res$proportions
        }
    }
    list(summary = do.call(rbind, rows), genes = do.call(rbind, genes),
         proportions = do.call(rbind, props))
}

#' Composition row with nothing assessed yet
#' @noRd
empty_composition_row <- function(key) {
    data.frame(comparison = key, share_A = NA_real_, share_B = NA_real_, share_fdr = NA_real_,
               composition = 'not assessed', comp_child = NA_character_,
               comp_share_A = NA_real_, comp_share_B = NA_real_, comp_min_fdr = NA_real_,
               exp_comp_deg = NA_real_, frac_comp_expected = NA_real_,
               n_comp_sensitive = NA_integer_, frac_comp_sensitive = NA_real_,
               n_parent_deg = NA_integer_, n_comp_driven = NA_integer_,
               frac_comp_driven = NA_real_, stringsAsFactors = FALSE)
}

#' The finest annotation level directly nested in `ann` (fewest clusters among
#' its children), or NULL
#' @noRd
immediate_child <- function(hier, ann) {
    if (is.null(hier)) return(NULL)
    m <- hier$map[hier$map$parent_annotation == ann, ]
    if (nrow(m) == 0) return(NULL)
    n_child <- tapply(m$child, m$child_annotation, function(z) length(unique(z)))
    names(n_child)[which.min(n_child)]
}

#' propeller on a count matrix (clusters x samples) for group B vs A; NULL when
#' it cannot be tested (fewer than two clusters or three samples, or a group
#' without samples)
#' @noRd
propeller_test <- function(counts, grp) {
    if (nrow(counts) < 2 || ncol(counts) < 3 || length(unique(grp)) < 2) return(NULL)
    prop <- sweep(counts, 2, colSums(counts), '/')
    prop_list <- list(Counts = counts, TransformedProps = asin(sqrt(prop)), Proportions = prop)
    design <- stats::model.matrix(~ 0 + factor(grp, levels = c('A', 'B')))
    out <- tryCatch(suppressMessages(suppressWarnings(speckle::propeller.ttest(
        prop_list, design, contrasts = c(-1, 1), robust = TRUE, trend = FALSE, sort = FALSE))),
        error = function(e) NULL)
    if (is.null(out)) return(NULL)
    data.frame(cluster = rownames(counts),
               share_A = rowMeans(prop[, grp == 'A', drop = FALSE]),
               share_B = rowMeans(prop[, grp == 'B', drop = FALSE]),
               fdr = out$FDR, stringsAsFactors = FALSE)
}

#' Every cluster's share of each sample's cells, A vs B (part 1)
#' @noRd
cluster_share <- function(cell_counts, A, B) {
    cc <- cell_counts[cell_counts$group %in% c(A, B), ]
    if (nrow(cc) == 0) return(NULL)
    m <- stats::xtabs(nCells ~ cluster + sample, data = cc)
    m <- matrix(m, nrow(m), dimnames = dimnames(m))
    s_grp <- tapply(cc$group, cc$sample, function(z) z[1])[colnames(m)]
    propeller_test(m, ifelse(s_grp == B, 'B', 'A'))
}

#' Parts 2-4 for one parent comparison
#' @noRd
parent_composition <- function(se, child_se, hier, ann, child_ann, cl, A, B, key, fit,
                               alpha, pergene, power, fdr, reference, th) {
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    sel <- cd$cluster == cl & cd$group %in% c(A, B)
    s_grp <- stats::setNames(ifelse(cd$group[sel] == B, 'B', 'A'), cd$sample[sel])
    hc <- hier$cells[hier$cells$parent_annotation == ann & hier$cells$child_annotation == child_ann &
                     hier$cells$parent == cl & hier$cells$sample %in% names(s_grp), ]
    if (length(unique(hc$child)) < 2) return(list(status = 'single child'))

    m <- stats::xtabs(nCells ~ child + sample, data = hc)
    m <- matrix(m, nrow(m), dimnames = dimnames(m))
    grp <- unname(s_grp[colnames(m)])
    shift <- propeller_test(m, grp)
    prop <- sweep(m, 2, colSums(m), '/')
    p_a <- rowMeans(prop[, grp == 'A', drop = FALSE])
    p_b <- rowMeans(prop[, grp == 'B', drop = FALSE])
    proportions <- data.frame(comparison = key, sample = rep(colnames(m), each = nrow(m)),
                              group = rep(ifelse(grp == 'B', B, A), each = nrow(m)),
                              child = rep(rownames(m), ncol(m)), nCells = as.vector(m),
                              share = as.vector(prop), stringsAsFactors = FALSE)

    # part 3: the shift-only log2 fold change, each gene's power to be called DE at
    # it, and the genes whose MDE it clears
    g <- fit$genes
    expr <- child_expression(S4Vectors::metadata(child_se)$cluster_totals, rownames(m), g,
                             A, reference)
    mix <- function(p) as.vector(expr$c %*% p[colnames(expr$c)]) / sum(p[colnames(expr$c)] * expr$L)
    lfc_comp <- log2(mix(p_b)) - log2(mix(p_a))
    lfc_comp[!is.finite(lfc_comp)] <- NA
    mde <- if (is.na(alpha)) rep(NA_real_, length(g)) else mde_from_se(fit$se_log2, fit$df, alpha, power)
    sensitive <- !is.na(lfc_comp) & !is.na(mde) & abs(lfc_comp) > mde
    comp_power <- if (is.na(alpha)) rep(NA_real_, length(g))
                  else power_from_se(fit$se_log2, fit$df, abs(ifelse(is.na(lfc_comp), 0, lfc_comp)), alpha)

    # part 4: attribution of the parent's own DEGs
    logfc <- if (is.null(fit$logfc)) rep(NA_real_, length(g)) else unname(fit$logfc[g])
    padj <- if (is.null(fit$pvalue)) rep(NA_real_, length(g)) else stats::p.adjust(unname(fit$pvalue[g]), 'BH')
    parent_de <- !is.na(padj) & padj < fdr
    child_keys <- paste0(child_ann, ':', rownames(m), ':', B, '_vs_', A)
    child_fits <- Filter(function(f) !is.null(f) && !is.null(f$pvalue), pergene[child_keys])
    child_de <- rep(FALSE, length(g))
    for (f in child_fits) {
        de_k <- stats::p.adjust(f$pvalue, 'BH') < fdr
        hit <- de_k[g]
        child_de <- child_de | (!is.na(hit) & hit)
    }
    explained <- !is.na(lfc_comp) & !is.na(logfc) & sign(lfc_comp) == sign(logfc) &
                 abs(lfc_comp) >= th$compExplained * abs(logfc)
    label <- if (length(child_fits) == 0) rep('unresolved', length(g))
             else ifelse(child_de & explained, 'mixed',
                  ifelse(child_de, 'within-child DE',
                  ifelse(explained, 'composition-driven', 'unresolved')))
    label[!parent_de] <- NA_character_

    top <- if (is.null(shift)) NULL else shift[which.min(shift$fdr), ]
    n_deg <- sum(parent_de)
    n_comp <- sum(label == 'composition-driven', na.rm = TRUE)
    summary <- list(
        comp_child = if (is.null(top)) NA_character_ else top$cluster,
        comp_share_A = if (is.null(top)) NA_real_ else top$share_A,
        comp_share_B = if (is.null(top)) NA_real_ else top$share_B,
        comp_min_fdr = if (is.null(top)) NA_real_ else top$fdr,
        exp_comp_deg = sum(comp_power), frac_comp_expected = mean(comp_power),
        n_comp_sensitive = sum(sensitive), frac_comp_sensitive = mean(sensitive),
        n_parent_deg = n_deg, n_comp_driven = n_comp,
        frac_comp_driven = if (n_deg > 0) n_comp / n_deg else NA_real_)
    genes <- data.frame(comparison = key, gene = g, lfc_comp = lfc_comp, comp_power = comp_power, mde = mde,
                        sensitive = sensitive, logfc = logfc, padj = padj,
                        parent_de = parent_de, child_de = child_de, label = label,
                        stringsAsFactors = FALSE)
    list(status = 'assessed', summary = summary, genes = genes, proportions = proportions)
}

#' Per-cell mean count of each gene in each child (`c`, genes x children) and the
#' children's per-cell library sizes (`L`), pooled or from the reference group
#' (pooled for a child without reference-group cells)
#' @noRd
child_expression <- function(totals, children, genes, A, reference) {
    children <- intersect(children, colnames(totals$counts))
    cnt <- totals$counts[, children, drop = FALSE]
    cells <- totals$cells[children]
    if (reference == 'reference' && !is.null(totals$by_group[[A]])) {
        ref <- totals$by_group[[A]]
        ref_cells <- ref$cells[children]
        has <- children[!is.na(ref_cells) & ref_cells > 0]
        cnt[, has] <- ref$counts[, has, drop = FALSE]
        cells[has] <- ref$cells[has]
    }
    per_cell <- sweep(cnt, 2, cells, '/')
    list(c = per_cell[genes, , drop = FALSE], L = colSums(per_cell))
}
