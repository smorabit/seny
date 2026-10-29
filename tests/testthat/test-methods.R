test_that('deseq2 and voom run end to end and carry the method', {
    sce <- sim_sce(n_per_group = 12, phi_donor = 0.02, theta = 0.3,
                   n_genes = 200, cells_per = 60, seed = 2)
    for (mth in c('edgeR', 'deseq2', 'voom')) {
        res <- assessAll(sce, sample = 'donor', group = 'condition',
                         annotations = 'cell_type', method = mth, targetLFC = 1)
        d <- as.data.frame(res)
        expect_true(all(d$method == mth))
        expect_true(all(is.finite(d$sensitivity_fraction)))
        expect_true(all(d$call %in% c('Adequate', 'Underpowered', 'Caution', 'Invalid')))
    }
})

test_that('deseq2 dispersion is in the same ballpark as edgeR', {
    # sim_sce plants one dispersion for all genes, so there is no gene-to-gene
    # variation to correlate; check the typical magnitude agrees instead.
    sce <- sim_sce(n_per_group = 15, phi_donor = 0.03, theta = 0.3,
                   n_genes = 250, cells_per = 60, seed = 3)
    x <- pbSetup(sce, sample = 'donor', group = 'condition', annotations = 'cell_type')
    fe <- estimateVariance(x, 'edgeR')$pergene[[1]]
    fd <- estimateVariance(x, 'deseq2')$pergene[[1]]
    expect_true(length(intersect(fe$genes, fd$genes)) > 50)
    ratio <- stats::median(fd$phi_total) / stats::median(fe$phi_total)
    expect_true(ratio > 0.3 && ratio < 3)
})

test_that('voom has no NB decomposition but a finite SE and df', {
    sce <- sim_sce(n_per_group = 12, n_genes = 200, cells_per = 60, seed = 4)
    x <- pbSetup(sce, sample = 'donor', group = 'condition', annotations = 'cell_type')
    v <- estimateVariance(x, 'voom')
    expect_true(all(is.na(v$summary$median_phi_total)))
    fit <- v$pergene[[1]]
    expect_false(isTRUE(fit$nb))
    expect_true(all(is.finite(fit$se_log2)))
    expect_true(all(is.finite(fit$df)))
})

test_that('an underpowered sim is Underpowered under every method', {
    sce <- sim_sce(n_per_group = 3, phi_donor = 0.25, theta = 0.5,
                   n_genes = 200, cells_per = 30, seed = 5)
    for (mth in c('edgeR', 'deseq2', 'voom')) {
        res <- assessAll(sce, sample = 'donor', group = 'condition',
                         annotations = 'cell_type', method = mth, targetLFC = 0.5)
        expect_true(any(as.character(res$call) == 'Underpowered'))
    }
})

test_that('mean aggregation works with voom and the method/aggregation guards fire', {
    sce <- sim_sce(n_per_group = 12, n_genes = 200, cells_per = 60, seed = 6)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type', aggregation = 'mean', method = 'voom')
    expect_true(all(is.finite(as.data.frame(res)$sensitivity_fraction)))

    expect_error(assessAll(sce, sample = 'donor', group = 'condition',
                           annotations = 'cell_type', aggregation = 'mean', method = 'edgeR'),
                 'sum aggregation|pairs with limma')

    x <- pbSetup(sce, sample = 'donor', group = 'condition', annotations = 'cell_type')
    expect_error(estimateVariance(x, 'dreamlet'), 'post-MVP')
})
