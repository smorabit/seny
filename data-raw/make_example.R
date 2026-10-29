# generate the small example dataset shipped with seny. three cell types in
# three regimes (adequate, cell-count imbalance correlated with condition, and
# underpowered) so the traffic-light matrix shows a mix of calls.
# run from the package root: pixi run Rscript data-raw/make_example.R

suppressMessages({
    library(SingleCellExperiment)
    library(Matrix)
})

set.seed(11)

n_per_group <- 12
n_genes <- 400
n_de <- 25
log2fc <- 1.2

donors <- paste0('d', seq_len(2 * n_per_group))
cond <- rep(c('A', 'B'), each = n_per_group)
names(cond) <- donors

base <- rnorm(n_genes, log(25), 1)                       # per-gene log mean
de_genes <- seq_len(n_de)
de_sign <- sample(c(-1, 1), n_de, replace = TRUE)

sim_celltype <- function(ct, cells_a, cells_b, phi_donor, theta) {
    mats <- list(); meta <- list()
    for (dn in donors) {
        n_cell <- if (cond[dn] == 'A') cells_a else cells_b
        deff <- rnorm(n_genes, 0, sqrt(phi_donor))       # donor random effect
        logmu <- base + deff
        if (cond[dn] == 'B') logmu[de_genes] <- logmu[de_genes] + de_sign * log2fc * log(2)
        m <- matrix(rnbinom(n_genes * n_cell, mu = exp(logmu), size = 1 / theta), nrow = n_genes)
        mats[[length(mats) + 1]] <- m
        meta[[length(meta) + 1]] <- data.frame(
            donor = rep(dn, n_cell), condition = unname(cond[dn]),
            cell_type = ct, stringsAsFactors = FALSE)
    }
    list(counts = do.call(cbind, mats), meta = do.call(rbind, meta))
}

regimes <- list(
    CD4T = list(cells_a = 90, cells_b = 90, phi_donor = 0.008, theta = 0.2),
    CD8T = list(cells_a = 70, cells_b = 15, phi_donor = 0.01, theta = 0.25),
    NK   = list(cells_a = 25, cells_b = 25, phi_donor = 0.20, theta = 0.5))

parts <- lapply(names(regimes), function(ct) {
    r <- regimes[[ct]]
    sim_celltype(ct, r$cells_a, r$cells_b, r$phi_donor, r$theta)
})

counts <- do.call(cbind, lapply(parts, `[[`, 'counts'))
cdat <- do.call(rbind, lapply(parts, `[[`, 'meta'))
rownames(counts) <- paste0('g', seq_len(n_genes))
colnames(counts) <- paste0('c', seq_len(ncol(counts)))
rownames(cdat) <- colnames(counts)

sce <- SingleCellExperiment(
    assays = list(counts = as(counts, 'CsparseMatrix')),
    colData = cdat)

dir.create('inst/extdata', recursive = TRUE, showWarnings = FALSE)
saveRDS(sce, 'inst/extdata/example_sce.rds')
cat('wrote inst/extdata/example_sce.rds:', ncol(sce), 'cells,', nrow(sce), 'genes\n')
