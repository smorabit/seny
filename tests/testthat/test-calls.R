test_that('a well-powered clean design is Adequate', {
    sce <- sim_sce(n_per_group = 20, phi_donor = 0.01, theta = 0.2,
                   n_genes = 200, cells_per = 80, seed = 3)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type', targetLFC = 1)
    expect_true(all(as.character(res$call) == 'Adequate'))
})

test_that('few donors with high donor variance is Underpowered', {
    sce <- sim_sce(n_per_group = 3, phi_donor = 0.25, theta = 0.5,
                   n_genes = 200, cells_per = 30, seed = 4)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type', targetLFC = 0.5)
    expect_true(any(as.character(res$call) == 'Underpowered'))
    expect_false(any(as.character(res$call) == 'Adequate'))
})

test_that('a batch nested in condition is Invalid', {
    sce <- sim_sce(n_per_group = 10, n_genes = 200, cells_per = 50,
                   batch = 'nested', seed = 5)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type', covariates = 'batch')
    expect_true(all(as.character(res$call) == 'Invalid'))
})

test_that('a partially confounded batch is Caution', {
    sce <- sim_sce(n_per_group = 15, phi_donor = 0.01, n_genes = 200,
                   cells_per = 60, batch = 'partial', seed = 6)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type', covariates = 'batch', targetLFC = 1)
    expect_true(all(as.character(res$call) == 'Caution'))
})

test_that('cell-count imbalance correlated with condition is Caution', {
    sce <- sim_sce(n_per_group = 15, phi_donor = 0.01, theta = 0.2, n_genes = 200,
                   cells_per = 80, cells_imbalance = 0.2, seed = 7)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type', targetLFC = 1)
    expect_true(all(as.character(res$call) == 'Caution'))
    expect_true(all(grepl('imbalance', res$reason)))
})

test_that('too few testable genes is Underpowered', {
    cmp <- 'a:c:B_vs_A'
    design <- data.frame(comparison = cmp, estimable = TRUE, nA = 10, nB = 10,
                         max_cramerV = 0, max_vif = 1, stringsAsFactors = FALSE)
    pseudobulk <- data.frame(comparison = cmp, cells_cond_cor = 0, n_testable = 50,
                             stringsAsFactors = FALSE)
    variance <- data.frame(comparison = cmp, stringsAsFactors = FALSE)
    sensitivity <- data.frame(comparison = cmp, sensitivity_fraction = 0.9,
                              target_lfc = 1, stringsAsFactors = FALSE)
    out <- callComparisons(design, pseudobulk, variance, sensitivity)
    expect_equal(as.character(out$call), 'Underpowered')
    expect_true(grepl('testable', out$reason))
})

test_that('a group emptied by the minCells filter is Invalid, not an error', {
    # B donors have 6 cells per cell type, all below minCells = 10
    sce <- sim_sce(n_per_group = 6, n_genes = 200, cells_per = 60,
                   cells_imbalance = 0.1, seed = 7)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type')
    expect_true(all(as.character(res$call) == 'Invalid'))
    expect_true(all(res$nB == 0))
    expect_true(all(is.na(res$cells_ratio)))
})

test_that('a cluster with no pseudobulk above minCells is Invalid, not dropped', {
    # every donor has 6 cells per cell type, all below minCells = 10
    sce <- sim_sce(n_per_group = 4, n_genes = 200, cells_per = 6, seed = 8)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type')
    expect_equal(nrow(res), 2)
    expect_true(all(as.character(res$call) == 'Invalid'))
})

test_that('thetaReps = 0 skips theta without changing the call', {
    sce <- sim_sce(n_per_group = 6, n_genes = 200, cells_per = 40, seed = 9)
    with_theta <- assessAll(sce, sample = 'donor', group = 'condition',
                            annotations = 'cell_type', thetaReps = 3)
    no_theta <- assessAll(sce, sample = 'donor', group = 'condition',
                          annotations = 'cell_type', thetaReps = 0)
    expect_identical(as.character(no_theta$call), as.character(with_theta$call))
    expect_equal(no_theta$sensitivity_fraction, with_theta$sensitivity_fraction)
    expect_equal(no_theta$median_phi_total, with_theta$median_phi_total)
    expect_true(all(is.na(no_theta$share_donor)))
    expect_true(all(is.na(no_theta$regime)))
    expect_error(pbSetup(sce, sample = 'donor', group = 'condition',
                         annotations = 'cell_type', thetaReps = -1), 'non-negative')
})

test_that('the cell-imbalance cutoff scales with sample size', {
    expect_equal(round(cells_cond_cutoff(6, 0.5, 0.05), 2), 0.81)
    expect_equal(cells_cond_cutoff(80, 0.5, 0.05), 0.5)
    expect_equal(cells_cond_cutoff(6, 0.5, NULL), 0.5)

    call_at <- function(n_per, r) {
        cmp <- 'a:c:B_vs_A'
        design <- data.frame(comparison = cmp, estimable = TRUE, nA = n_per, nB = n_per,
                             max_cramerV = 0, max_vif = 1, stringsAsFactors = FALSE)
        pseudobulk <- data.frame(comparison = cmp, cells_cond_cor = r, n_testable = 500,
                                 stringsAsFactors = FALSE)
        variance <- data.frame(comparison = cmp, stringsAsFactors = FALSE)
        sensitivity <- data.frame(comparison = cmp, sensitivity_fraction = 0.9,
                                  target_lfc = 1, stringsAsFactors = FALSE)
        as.character(callComparisons(design, pseudobulk, variance, sensitivity)$call)
    }
    # r = 0.6 is chance-level at 3 vs 3 but a real imbalance at 20 vs 20
    expect_equal(call_at(3, 0.6), 'Adequate')
    expect_equal(call_at(20, 0.6), 'Caution')
})

test_that('the call uses the larger of the FDR-based and fixed-alpha MDEs', {
    sce <- sim_sce(n_per_group = 6, n_genes = 200, cells_per = 40, seed = 12)
    x <- pbSetup(sce, sample = 'donor', group = 'condition', annotations = 'cell_type',
                 thetaReps = 0)
    v <- estimateVariance(x, method = 'edgeR')
    both <- estimateSensitivity(v, targetLFC = 1, pi1 = 0.05, alphaFixed = 0.001)
    fdr_only <- estimateSensitivity(v, targetLFC = 1, pi1 = 0.05)
    # alpha_fdr ~ 0.0022 > 0.001, so the fixed-alpha MDE is the larger one
    expect_equal(both$median_mde, both$median_mde_fixed)
    expect_true(all(both$median_mde_fixed > both$median_mde_fdr))
    expect_true(all(both$sensitivity_fraction <= fdr_only$sensitivity_fraction))
    expect_equal(fdr_only$median_mde, fdr_only$median_mde_fdr)
    expect_true(all(is.na(fdr_only$median_mde_fixed)))
    # a lenient fixed alpha leaves the FDR-based MDE in charge
    loose <- estimateSensitivity(v, targetLFC = 1, pi1 = 0.05, alphaFixed = 0.01)
    expect_equal(loose$median_mde, loose$median_mde_fdr)
})

test_that('pi1 = "discoveries" uses each comparison\'s fraction of BH discoveries', {
    sce <- sim_sce(n_per_group = 6, n_genes = 300, n_de = 60, cells_per = 40, seed = 12)
    x <- pbSetup(sce, sample = 'donor', group = 'condition', annotations = 'cell_type',
                 thetaReps = 0)
    for (cur_method in c('edgeR', 'deseq2', 'voom')) {
        v <- estimateVariance(x, method = cur_method)
        s <- estimateSensitivity(v, targetLFC = 1)
        for (key in names(v$pergene)) {
            p <- v$pergene[[key]]$pvalue
            p <- p[!is.na(p)]
            expect_gt(length(p), 0)
            expected <- max(sum(p.adjust(p, 'BH') < 0.05) / length(p), 1 / length(p))
            expect_equal(s$pi1[s$comparison == key], expected)
        }
        # alpha follows the Jung mapping at that pi1
        r <- s[1, ]
        expect_equal(r$alpha, 0.05 * 0.8 * r$pi1 / (0.95 * (1 - r$pi1)))
    }
})

test_that('more true DE genes raise pi1, alpha and sensitivity under "discoveries"', {
    fit_of <- function(n_de) {
        sce <- sim_sce(n_per_group = 6, n_genes = 300, n_de = n_de, cells_per = 40, seed = 3)
        x <- pbSetup(sce, sample = 'donor', group = 'condition', annotations = 'cell_type',
                     thetaReps = 0)
        estimateVariance(x, method = 'edgeR')
    }
    v_sparse <- fit_of(3)
    v_dense <- fit_of(90)
    sparse <- estimateSensitivity(v_sparse, targetLFC = 1)
    dense <- estimateSensitivity(v_dense, targetLFC = 1)
    expect_true(all(dense$pi1 > sparse$pi1))
    expect_true(all(dense$alpha > sparse$alpha))
    expect_true(all(dense$sensitivity_fraction >= sparse$sensitivity_fraction))
    # a fixed pi1 ignores the contrast
    expect_true(all(estimateSensitivity(v_dense, pi1 = 0.05)$pi1 == 0.05))
})

test_that('the pi1-mapped alpha never exceeds the FDR level', {
    a <- sensitivity_alpha(power = 0.8, fdr = 0.05, pi1 = 0.9)
    expect_equal(a[['fdr']], 0.05)
    expect_equal(a[['call']], 0.05)
    expect_true(is.na(a[['fixed']]))
    expect_error(estimateSensitivity(list(pergene = list(), summary = data.frame(comparison = character(0))),
                                     pi1 = 'storey'), 'discoveries')
})
