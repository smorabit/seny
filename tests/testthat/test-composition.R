comp_run <- function(shift, seed = 21) {
    sce <- sim_sce(n_per_group = 12, phi_donor = 0.01, theta = 0.2, n_genes = 300, n_de = 0,
                   cells_per = 120, subtype_shift = shift, seed = seed)
    assessAll(sce, sample = 'donor', group = 'condition',
              annotations = c('cell_type', 'cell_type_fine'), targetLFC = 1, thetaReps = 0)
}

test_that('the hierarchy is inferred from the cells, coarser level as parent', {
    sce <- sim_sce(n_per_group = 4, n_genes = 50, cells_per = 20, n_de = 0,
                   subtype_shift = 0.2, seed = 2)
    x <- pbSetup(sce, sample = 'donor', group = 'condition',
                 annotations = c('cell_type_fine', 'cell_type'), thetaReps = 0)
    h <- S4Vectors::metadata(x)$hierarchy
    expect_equal(unique(h$map$parent_annotation), 'cell_type')
    expect_equal(h$map$parent[match(c('T1', 'T2', 'B'), h$map$child)], c('T', 'T', 'B'))
    expect_true(all(h$map$purity == 1))
    expect_null(S4Vectors::metadata(pbSetup(sce, sample = 'donor', group = 'condition',
                                            annotations = 'cell_type', thetaReps = 0))$hierarchy)
})

test_that('a subtype shift without within-subtype DE is flagged and attributed', {
    res <- comp_run(0.3)
    t_row <- res[res$comparison == 'cell_type:T:B_vs_A', ]
    expect_equal(as.character(t_row$call), 'Caution')
    expect_match(t_row$reason, 'composition shift')
    expect_lt(t_row$comp_min_fdr, 0.05)
    g <- attr(res, 'composition')$genes
    g <- g[g$comparison == 'cell_type:T:B_vs_A', ]
    markers <- paste0('g', 1:40)
    # the shift-only change of the markers clears their MDE, and has the right sign
    expect_gt(mean(g$sensitive[g$gene %in% markers]), 0.8)
    expect_true(all(g$lfc_comp[g$gene %in% paste0('g', 1:20)] > 0))
    expect_true(all(g$lfc_comp[g$gene %in% paste0('g', 21:40)] < 0))
    expect_lt(mean(g$sensitive[!g$gene %in% markers]), 0.05)
    # the expected composition DEGs are carried by the markers
    expect_gt(sum(g$comp_power[g$gene %in% markers]), 30)
    expect_lt(mean(g$comp_power[!g$gene %in% markers]), 0.05)
    expect_equal(t_row$exp_comp_deg, sum(g$comp_power))
    # the parent's DEGs are the markers, attributed to composition
    deg <- g[g$parent_de, ]
    expect_gt(nrow(deg), 10)
    expect_gt(mean(deg$label[deg$gene %in% markers] == 'composition-driven'), 0.8)
    expect_true(all(is.na(g$label[!g$parent_de])))
    # the subtypes themselves, B (single subtype) and the fine level are not flagged
    others <- res[res$comparison != 'cell_type:T:B_vs_A', ]
    expect_false(any(grepl('composition shift', others$reason)))
    expect_equal(res$composition[res$comparison == 'cell_type:B:B_vs_A'], 'single child')
    expect_true(all(res$composition[res$annotation == 'cell_type_fine'] == 'not assessed'))
})

test_that('no subtype shift, no composition flag', {
    res <- comp_run(0)
    expect_false(any(grepl('composition shift', res$reason)))
    t_row <- res[res$comparison == 'cell_type:T:B_vs_A', ]
    expect_equal(t_row$composition, 'assessed')
    expect_gt(t_row$comp_min_fdr, 0.05)
})

test_that('a single annotation is not assessed and calls are unchanged', {
    sce <- sim_sce(n_per_group = 10, phi_donor = 0.01, theta = 0.2, n_genes = 200,
                   cells_per = 60, seed = 8)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type', thetaReps = 0)
    expect_true(all(res$composition == 'not assessed'))
    x <- attr(res, 'pb')
    v <- estimateVariance(x)
    s <- estimateSensitivity(v)
    base <- callComparisons(assessDesign(x), assessPseudobulk(x), v$summary, s)
    expect_identical(as.character(res$call), as.character(base$call))
    expect_identical(res$reason, base$reason)
    # the cluster share is still reported
    expect_true(all(!is.na(res$share_A)))
})

test_that('the per-gene fit carries the observed log2FC of the edgeR test', {
    sce <- sim_sce(n_per_group = 6, n_genes = 150, cells_per = 40, seed = 9)
    x <- pbSetup(sce, sample = 'donor', group = 'condition', annotations = 'cell_type',
                 thetaReps = 0)
    f <- estimateVariance(x)$pergene[['cell_type:T:B_vs_A']]
    d <- comparison_data(x[[1]], 'T', 'A', 'B')
    y <- edgeR::DGEList(d$counts)
    design <- stats::model.matrix(~ d$group)
    y <- y[edgeR::filterByExpr(y, group = d$group), , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    qf <- edgeR::glmQLFit(edgeR::estimateDisp(y, design, robust = TRUE), design, robust = TRUE)
    ref <- edgeR::glmQLFTest(qf, coef = 2)$table$logFC
    expect_equal(unname(f$logfc[rownames(y)]), ref, tolerance = 1e-8)
})

test_that('real within-subtype DE is not attributed to composition', {
    sce <- sim_sce(n_per_group = 12, phi_donor = 0.01, theta = 0.2, n_genes = 300, n_de = 15,
                   log2fc = 1.5, cells_per = 120, subtype_shift = 0.3, seed = 22)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = c('cell_type', 'cell_type_fine'), thetaReps = 0)
    g <- attr(res, 'composition')$genes
    g <- g[g$comparison == 'cell_type:T:B_vs_A' & g$parent_de, ]
    real <- g[g$gene %in% paste0('g', 1:15), ]
    markers <- g[g$gene %in% paste0('g', 16:55), ]
    expect_gt(mean(real$label == 'within-child DE'), 0.8)
    expect_gt(mean(markers$label == 'composition-driven'), 0.8)
})
