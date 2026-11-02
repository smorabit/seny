# generate the small example dataset shipped with seny. three cell types in
# three regimes (adequate, cell-count imbalance correlated with condition, and
# underpowered) so the traffic-light matrix shows a mix of calls. a finer
# annotation (cell_type_fine) splits CD4 T into naive and memory subtypes, each
# with its own markers, whose shares shift with condition (naive 60% in A, 35%
# in B): at the broad level the markers look DE, the composition flag's case.
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
n_marker <- 20
marker_fc <- 3

sim_celltype <- function(ct, cells_a, cells_b, phi_donor, theta, subtypes = NULL) {
    mats <- list(); meta <- list()
    for (dn in donors) {
        n_cell <- if (cond[dn] == 'A') cells_a else cells_b
        deff <- rnorm(n_genes, 0, sqrt(phi_donor))       # donor random effect
        logmu <- base + deff
        if (cond[dn] == 'B') logmu[de_genes] <- logmu[de_genes] + de_sign * log2fc * log(2)
        fine <- rep(ct, n_cell)
        if (!is.null(subtypes)) {
            share <- if (cond[dn] == 'A') subtypes$share_a else subtypes$share_b
            n_first <- round(n_cell * share)
            fine <- rep(subtypes$names, c(n_first, n_cell - n_first))
        }
        m <- do.call(cbind, lapply(unique(fine), function(cur_fine) {
            cur_logmu <- logmu
            k <- match(cur_fine, subtypes$names)
            if (!is.na(k)) {
                markers <- n_de + (k - 1) * n_marker + seq_len(n_marker)
                cur_logmu[markers] <- cur_logmu[markers] + log(marker_fc)
            }
            matrix(rnbinom(n_genes * sum(fine == cur_fine), mu = exp(cur_logmu), size = 1 / theta),
                   nrow = n_genes)
        }))
        mats[[length(mats) + 1]] <- m
        meta[[length(meta) + 1]] <- data.frame(
            donor = rep(dn, n_cell), condition = unname(cond[dn]),
            cell_type = ct, cell_type_fine = fine, stringsAsFactors = FALSE)
    }
    list(counts = do.call(cbind, mats), meta = do.call(rbind, meta))
}

regimes <- list(
    CD4T = list(cells_a = 90, cells_b = 90, phi_donor = 0.008, theta = 0.2,
                subtypes = list(names = c('CD4 naive', 'CD4 memory'), share_a = 0.6, share_b = 0.35)),
    CD8T = list(cells_a = 70, cells_b = 15, phi_donor = 0.01, theta = 0.25),
    NK   = list(cells_a = 20, cells_b = 20, phi_donor = 0.30, theta = 0.5))

parts <- lapply(names(regimes), function(ct) {
    r <- regimes[[ct]]
    sim_celltype(ct, r$cells_a, r$cells_b, r$phi_donor, r$theta, r$subtypes)
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
