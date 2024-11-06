### read data

morris_small <- readRDS("../data/morris_raw_small_screen.rds")

# set thresholds for nUMI, nGenes, %mito

mito_genes <- grepl("^MT-", rowData(morris_small)$Symbol)

morris_small$n_umi = colSums(counts(morris_small))
morris_small$n_genes = colSums(counts(morris_small) > 0)
morris_small$perc_mito = colSums(counts(morris_small[mito_genes, ])) / colSums(counts(morris_small))

colData(morris_small) %>% data.frame() %>%
  ggplot(aes(x = n_umi, y = perc_mito)) + geom_point(size = .1) + geom_vline(xintercept = 5000, col = "red") + geom_hline(yintercept = .1, col = "red")

colData(morris_small) %>%
  ggplot(aes(x = n_genes, y = perc_mito)) + geom_point(size = .1) + geom_vline(xintercept = 2000, col = "red") + geom_hline(yintercept = .1, col = "red")

thresh_n_umi <- 5000
thresh_n_genes <- 2000
thresh_perc_mito <- .1

morris_small <- morris_small[ , morris_small$n_umi > thresh_n_umi & morris_small$n_genes > thresh_n_genes & morris_small$perc_mito < thresh_perc_mito]

# assign gRNAs

hist(log10(colSums(counts(altExps(morris_small)[["cre_pert"]]))))
