# small single-cell simulator for the tests: two cell types, two conditions,
# NB cells with a donor random effect (phi_donor) and within-donor overdispersion
# (theta), a handful of DE genes, and an optional confounding batch.
# subtype_shift adds a fine annotation: T cells split into T1/T2, each with its
# own marker genes (no DE within either), T1's share 0.5 in A and
# 0.5 + subtype_shift in B; B cells stay one subtype.

sim_sce <- function(n_per_group = 12, phi_donor = 0.03, theta = 0.3, n_genes = 300,
                    n_de = 15, log2fc = 1.5, cells_per = 60, cells_imbalance = 1,
                    batch = c('none', 'partial', 'nested'), subtype_shift = NULL,
                    n_marker = 20, marker_fc = 3, seed = 1) {
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
    t1_markers <- n_de + seq_len(n_marker)
    t2_markers <- n_de + n_marker + seq_len(n_marker)

    cts <- c('T', 'B')
    mats <- list(); meta <- list()
    for (ct in cts) {
        for (dn in donors) {
            n_cell <- if (cond[dn] == 'B') round(cells_per * cells_imbalance) else cells_per
            deff <- stats::rnorm(n_genes, 0, sqrt(phi_donor))          # donor random effect
            lm <- base + deff
            if (cond[dn] == 'B') lm[de_genes] <- lm[de_genes] + de_sign * log2fc * log(2)
            fine <- rep(ct, n_cell)
            if (!is.null(subtype_shift) && ct == 'T') {
                share_t1 <- 0.5 + if (cond[dn] == 'B') subtype_shift else 0
                fine <- rep(c('T1', 'T2'), c(round(n_cell * share_t1), n_cell - round(n_cell * share_t1)))
            }
            m <- do.call(cbind, lapply(unique(fine), function(cur_fine) {
                cur_lm <- lm
                if (cur_fine == 'T1') cur_lm[t1_markers] <- cur_lm[t1_markers] + log(marker_fc)
                if (cur_fine == 'T2') cur_lm[t2_markers] <- cur_lm[t2_markers] + log(marker_fc)
                matrix(stats::rnbinom(n_genes * sum(fine == cur_fine), mu = exp(cur_lm), size = 1 / theta),
                       nrow = n_genes)
            }))
            mats[[length(mats) + 1]] <- m
            meta[[length(meta) + 1]] <- data.frame(
                donor = rep(dn, n_cell), condition = unname(cond[dn]),
                batch = unname(bat[dn]), cell_type = ct,
                cell_type_fine = fine, stringsAsFactors = FALSE)
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
