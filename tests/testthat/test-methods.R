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

test_that('the edgeR null fit keeps genes expressed in only one condition', {
    sce <- sim_sce(n_per_group = 6, n_genes = 200, cells_per = 40, seed = 10)
    cnt <- SummarizedExperiment::assay(sce, 'counts')
    # a gene silent in condition A and strongly expressed in B (a Y-gene analogue)
    on <- sce$condition == 'B'
    cnt[1, ] <- 0
    cnt[1, on] <- stats::rpois(sum(on), 5)
    SummarizedExperiment::assay(sce, 'counts') <- cnt
    x <- pbSetup(sce, sample = 'donor', group = 'condition', annotations = 'cell_type',
                 thetaReps = 0)
    v <- estimateVariance(x, method = 'edgeR')
    expect_true(all(vapply(v$pergene, function(f) rownames(cnt)[1] %in% f$genes, logical(1))))
})

test_that('the edgeR SE uses full-model dispersion; the decomposition stays null-model', {
    fit_at <- function(fc) {
        sce <- sim_sce(n_per_group = 4, n_genes = 200, n_de = 15, log2fc = fc,
                       cells_per = 40, seed = 11)
        x <- pbSetup(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type', thetaReps = 0)
        estimateVariance(x, method = 'edgeR')$pergene[[1]]
    }
    f <- fit_at(3)      # 15 strong DE genes (g1..g15, log2FC 3)
    f0 <- fit_at(0)     # same data without the effect: the true noise level
    de <- f$genes %in% paste0('g', 1:15)
    null0 <- stats::median(f0$phi_se[!(f0$genes %in% paste0('g', 1:15))])
    # the null fit absorbs the effect as noise, for the DE genes and, through the
    # shared dispersion trend, for the rest; the full fit recovers the true level
    expect_lt(stats::median(f$phi_se[de]), stats::median(f$phi_total[de]) / 5)
    expect_gt(stats::median(f$phi_total[!de]), 1.5 * null0)
    expect_lt(abs(log(stats::median(f$phi_se[!de]) / null0)), 0.2)
    # without DE the two fits agree
    expect_lt(abs(log(stats::median(f0$phi_se) / stats::median(f0$phi_total))), 0.05)
})
