# small single-cell simulator for the tests: two cell types, two conditions,
# NB cells with a donor random effect (phi_donor) and within-donor overdispersion
# (theta), a handful of DE genes, and an optional confounding batch.

sim_sce <- function(n_per_group = 12, phi_donor = 0.03, theta = 0.3, n_genes = 300,
                    n_de = 15, log2fc = 1.5, cells_per = 60, cells_imbalance = 1,
                    batch = c('none', 'partial', 'nested'), seed = 1) {
    batch <- match.arg(batch)
    set.seed(seed)
    suppressMessages({
        requireNamespace('SingleCellExperiment')
        requireNamespace('Matrix')
    })
    n_donor <- 2 * n_per_group
    donors <- paste0('d', seq_len(n_donor))
    cond <- rep(c('A', 'B'), each = n_per_group)
    if (batch == 'nested') {
        bat <- ifelse(cond == 'A', 'b1', 'b2')
    } else if (batch == 'partial') {
        bat <- ifelse(cond == 'A', 'b1', 'b2')
        flip <- sample(n_donor, round(0.2 * n_donor))
        bat[flip] <- ifelse(bat[flip] == 'b1', 'b2', 'b1')
    } else {
        bat <- sample(c('b1', 'b2'), n_donor, replace = TRUE)
    }
    names(cond) <- names(bat) <- donors

    base <- stats::rnorm(n_genes, log(20), 1)             # per-gene log mean
    de_genes <- seq_len(n_de)
    de_sign <- sample(c(-1, 1), n_de, replace = TRUE)

    cts <- c('T', 'B')
    mats <- list(); meta <- list()
    for (ct in cts) {
        for (dn in donors) {
            n_cell <- if (cond[dn] == 'B') round(cells_per * cells_imbalance) else cells_per
            deff <- stats::rnorm(n_genes, 0, sqrt(phi_donor))          # donor random effect
            lm <- base + deff
            if (cond[dn] == 'B') lm[de_genes] <- lm[de_genes] + de_sign * log2fc * log(2)
            mu <- exp(lm)
            m <- matrix(stats::rnbinom(n_genes * n_cell, mu = mu, size = 1 / theta),
                        nrow = n_genes)
            mats[[length(mats) + 1]] <- m
            meta[[length(meta) + 1]] <- data.frame(
                donor = rep(dn, n_cell), condition = unname(cond[dn]),
                batch = unname(bat[dn]), cell_type = ct, stringsAsFactors = FALSE)
        }
    }
    counts <- do.call(cbind, mats)
    cdat <- do.call(rbind, meta)
    rownames(counts) <- paste0('g', seq_len(n_genes))
    colnames(counts) <- paste0('c', seq_len(ncol(counts)))
    rownames(cdat) <- colnames(counts)
    SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = as(counts, 'CsparseMatrix')),
        colData = cdat)
}
