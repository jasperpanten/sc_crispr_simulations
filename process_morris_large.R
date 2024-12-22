library(scater)
library(scran)
library(tidyverse)
library(DropletUtils)
library(cluster)

source("./differential_expression_fun.R")
source("./power_simulations_fun.R")

setwd(paste0(dirname(rstudioapi::getSourceEditorContext()$path)))

get_quantiles <- function(values, qlow = .1, qhigh = .9){
  return(c(quantile(values, qlow), quantile(values, qhigh)))
}

for (i in c(1, 2, 3, 4)){
  
  morris_dataset_1 <- readRDS(paste0("../data/morris_raw_small_screen_", i, ".rds"))
  
  # calculate qc metrics and set filters
  mito_genes <- rowData(morris_dataset_1) %>% data.frame() %>% dplyr::filter(grepl("^MT-", Symbol)) %>% pull(ID)
  colData(morris_dataset_1) <- cbind(colData(morris_dataset_1), perCellQCMetrics(morris_dataset_1, subset = list("mito" = mito_genes)))
  
  quantiles_sum <- get_quantiles(morris_dataset_1$sum, 0.10, 0.99)
  qplot(x = "1", y = morris_dataset_1$sum) + geom_violin() + geom_hline(yintercept = quantiles_sum, col = "red")
  
  quantiles_detected <- get_quantiles(morris_dataset_1$detected, 0.01, 0.99)
  qplot(x = "1", y = morris_dataset_1$detected) + geom_violin() + geom_hline(yintercept = quantiles_detected, col = "red")
  
  quantiles_mito <- get_quantiles(morris_dataset_1$subsets_mito_percent, 0.01, 0.90)
  qplot(x = "1", y = morris_dataset_1$subsets_mito_percent) + geom_violin() + geom_hline(yintercept = quantiles_mito, col = "red")
  
  cells_keep_quality <- (morris_dataset_1$sum > quantiles_sum[[1]] & morris_dataset_1$sum < quantiles_sum[[2]]) & 
    (morris_dataset_1$detected > quantiles_detected[[1]] & morris_dataset_1$detected < quantiles_detected[[2]]) & 
    (morris_dataset_1$subsets_mito_percent > quantiles_mito[[1]] & morris_dataset_1$subsets_mito_percent < quantiles_mito[[2]])
  
  # 
  
  center_log_transformation <- function(mat, pseudocount = 1){
    mat = log2(mat + pseudocount)
    apply(mat, 2, function(x){
      x - mean(x)
    })
  }
  
  logcounts(altExps(morris_dataset_1)[["hto"]]) <- center_log_transformation(counts(altExps(morris_dataset_1)[["hto"]]))
  
  MaxN <- function(x, N = 1){
    return(sort(x)[[N]])
  }
  
  # reimplementation of HTODemux for scater objects
  HTODemux_scater <- function(
    object,
    assay = "hto",
    count_transformation = "counts",
    positive.quantile = 0.99,
    init = NULL,
    nstarts = 100,
    kfunc = "clara",
    nsamples = 100,
    seed = 42,
    verbose = TRUE
  ) {
    if (!is.null(x = seed)) {
      set.seed(seed = seed)
    }
    #initial clustering
    # assay <- assay %||% DefaultAssay(object = object)
    # data <- GetAssayData(object = object, assay = assay)
    # counts <- GetAssayData(
    #   object = object,
    #   assay = assay,
    #   slot = 'counts'
    # )[, colnames(x = object)]
    # counts <- as.matrix(x = counts)
    
    data <-  altExps(object)[[assay]]
    # object <- altExps(object)[[assay]]
    counts <- assays(data)[[count_transformation]]
    
    # ncenters <- init %||% (nrow(x = data) + 1)
    ncenters = nrow(counts) + 1
    
    # switch(
    #   EXPR = kfunc,
    #   'kmeans' = {
    #     init.clusters <- kmeans(
    #       x = t(x = GetAssayData(object = object, assay = assay)),
    #       centers = ncenters,
    #       nstart = nstarts
    #     )
    #     #identify positive and negative signals for all HTO
    #     Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
    #   },
    #   'clara' = {
    #     #use fast k-medoid clustering
    #     init.clusters <- clara(
    #       x = t(x = GetAssayData(object = object, assay = assay)),
    #       k = ncenters,
    #       samples = nsamples
    #     )
    #     #identify positive and negative signals for all HTO
    #     Idents(object = object, cells = names(x = init.clusters$clustering), drop = TRUE) <- init.clusters$clustering
    #   },
    #   stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'")
    # )
    
    switch(
      EXPR = kfunc,
      'kmeans' = {
        init.clusters <- kmeans(
          x = t(x = counts),
          centers = ncenters,
          nstart = nstarts
        )
        #identify positive and negative signals for all HTO
        data$clustering <- init.clusters$cluster
      },
      'clara' = {
        #use fast k-medoid clustering
        init.clusters <- clara(
          x = t(x = counts),
          k = ncenters,
          samples = nsamples
        )
        #identify positive and negative signals for all HTO
        # Idents(object = object, cells = names(x = init.clusters$clustering), drop = TRUE) <- init.clusters$clustering
        data$clustering <- init.clusters$clustering
      },
      stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'")
    )
    
    #average hto signals per cluster
    #work around so we don't average all the RNA levels which takes time
    # average.expression <- suppressWarnings(
    #   AverageExpression(
    #     object = object,
    #     assays = assay,
    #     verbose = FALSE
    #   )[[assay]]
    # )
    
    average.expression <- assays(aggregateAcrossCells(data, data$clustering, statistics = "mean"))[[count_transformation]]
    
    print(average.expression)
    
    #checking for any cluster with all zero counts for any barcode
    if (any(colSums(average.expression == 0) > 0)) {
      stop("Cells with zero counts exist as a cluster.")
    }
    #create a matrix to store classification result
    discrete <- counts
    discrete[discrete > 0] <- 0
    # for each HTO, we will use the minimum cluster for fitting
    # for (iter in rownames(x = data)) {
    #   values <- counts[iter, colnames(data)]
    #   #commented out if we take all but the top cluster as background
    #   #values_negative=values[setdiff(object@cell.names,WhichCells(object,which.max(average.expression[iter,])))]
    #   values.use <- values[WhichCells(
    #     object = object,
    #     idents = levels(x = Idents(object = object))[[which.min(x = average.expression[iter, ])]]
    #   )]
    #   fit <- suppressWarnings(expr = fitdist(data = values.use, distr = "nbinom"))
    #   cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
    #   discrete[iter, names(x = which(x = values > cutoff))] <- 1
    #   if (verbose) {
    #     message(paste0("Cutoff for ", iter, " : ", cutoff, " reads"))
    #   }
    # }
    for (iter in rownames(x = data)) {
      values <- counts[iter, colnames(data)]
      #commented out if we take all but the top cluster as background
      #values_negative=values[setdiff(object@cell.names,WhichCells(object,which.max(average.expression[iter,])))]
      values.use <- values[data$clustering == which.min(x = average.expression[iter, ])]
      fit <- suppressWarnings(expr = fitdistrplus::fitdist(data = values.use, distr = "nbinom"))
      cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
      discrete[iter, which(x = values > cutoff)] <- 1
      if (verbose) {
        message(paste0("Cutoff for ", iter, " : ", cutoff, " reads"))
      }
    }
    # now assign cells to HTO based on discretized values
    npositive <- colSums(x = discrete)
    classification.global <- npositive
    classification.global[npositive == 0] <- "Negative"
    classification.global[npositive == 1] <- "Singlet"
    classification.global[npositive > 1] <- "Doublet"
    donor.id = rownames(x = data)
    hash.max <- apply(X = counts, MARGIN = 2, FUN = max)
    hash.maxID <- apply(X = counts, MARGIN = 2, FUN = which.max)
    hash.second <- apply(X = counts, MARGIN = 2, FUN = MaxN, N = 2)
    hash.maxID <- as.character(x = donor.id[sapply(
      X = 1:ncol(x = data),
      FUN = function(x) {
        return(which(x = counts[, x] == hash.max[x])[1])
      }
    )])
    hash.secondID <- as.character(x = donor.id[sapply(
      X = 1:ncol(x = data),
      FUN = function(x) {
        return(which(x = counts[, x] == hash.second[x])[1])
      }
    )])
    hash.margin <- hash.max - hash.second
    doublet_id <- sapply(
      X = 1:length(x = hash.maxID),
      FUN = function(x) {
        return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), collapse = "_"))
      }
    )
    # doublet_names <- names(x = table(doublet_id))[-1] # Not used
    classification <- classification.global
    classification[classification.global == "Negative"] <- "Negative"
    classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == "Singlet")]
    classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == "Doublet")]
    classification.metadata <- data.frame(
      hash.maxID,
      hash.secondID,
      hash.margin,
      classification,
      classification.global
    )
    colnames(x = classification.metadata) <- paste(
      assay,
      c('maxID', 'secondID', 'margin', 'classification', 'classification.global'),
      sep = '_'
    )
    colData(object) <- cbind(colData(object), classification.metadata)
    # Idents(object) <- paste0(assay, '_classification')
    # Idents(object, cells = rownames(object@meta.data[object@meta.data$classification.global == "Doublet", ])) <- "Doublet"
    # doublets <- rownames(x = object[[]])[which(object[[paste0(assay, "_classification.global")]] == "Doublet")]
    # Idents(object = object, cells = doublets) <- 'Doublet'
    # object@meta.data$hash.ID <- Idents(object)
    hash_ids <- classification.metadata$hto_classification
    hash_ids[classification.metadata$hto_classification.global == "Doublet"] <- "Doublet"
    hash_ids[classification.metadata$hto_classification.global == "Negative"] <- "Negative"
    object$hash.ID <- hash_ids
    return(object)
  }
  
  morris_dataset_1 <- HTODemux_scater(morris_dataset_1, count_transformation = "counts", kfunc = "clara")
  
  cells_keep_demultiplex <- morris_dataset_1$hash.ID != "Negative" & morris_dataset_1$hash.ID != "Doublet"
  
  # gdo thresholds
  quantiles_gdo <- get_quantiles(morris_dataset_1$altexps_cre_pert_sum, 0.01, 0.99)
  qplot(x = "1", y = morris_dataset_1$altexps_cre_pert_sum) + geom_violin() + geom_hline(yintercept = quantiles_gdo, col = "red") + scale_y_log10()
  
  cells_keep_gdo <- (morris_dataset_1$altexps_cre_pert_sum > quantiles_gdo[[1]] & morris_dataset_1$altexps_cre_pert_sum < quantiles_gdo[[2]])
  
  # filter cells for final dataset
  morris_dataset_1_filtered <- morris_dataset_1[ , cells_keep_quality & cells_keep_demultiplex & cells_keep_gdo ]
  
  # gdo annotations
  umi_thresholds <- read.csv("../data/GSE171452_RAW/GSM7108121_STINGseq-v2_GDO-A.umi_thresholds.txt") %>%
    dplyr::filter(grepl("^SNP-", Protospacer)) %>% column_to_rownames("Protospacer")
  
  grna_hits <- do.call("rbind", lapply(rownames(umi_thresholds), function(grna){
    threshold = umi_thresholds[grna, ]
    counts = counts(altExps(morris_dataset_1_filtered)[["cre_pert"]])[grna, ]
    return(as.numeric(counts > threshold & counts > 3))
  }))
  colnames(grna_hits) <- colnames(morris_dataset_1_filtered)
  rownames(grna_hits) <- rownames(umi_thresholds)
  
  # add target genes for each cre
  readxl::read_excel("~/Desktop/science.adh7699_tables_s1_to_s4(1)/science.adh7699_table_s3.xlsx", sheet = 6, skip = 2) %>% colnames()
  readxl::read_excel("~/Desktop/science.adh7699_tables_s1_to_s4(1)/science.adh7699_table_s3.xlsx", sheet = 6, skip = 2) %>% dim()
  
  readxl::read_excel("~/Desktop/science.adh7699_tables_s1_to_s4(1)/science.adh7699_table_s3.xlsx", sheet = 6, skip = 2) %>%
    # dplyr::select(c("SNP", "SNP Coordinates (hg19)", "Gene", "Ensembl ID", "gRNAs")) %>% 
    dplyr::select(c("Gene", "Ensembl ID", "gRNAs")) %>% 
    separate_longer_delim(gRNAs, delim = "__") %>% 
    group_by(gRNAs) %>% 
    summarize(target_genes = list(`Ensembl ID`), gene_symbol = list(Gene)) %>% 
    column_to_rownames("gRNAs") %>% 
    mutate(gRNAs = rownames(.)) %>%
    mutate(cre_target = gsub("-[0-9]*$", "", gRNAs)) %>% 
    distinct() -> annotations
  
  grna_hits <- grna_hits[rownames(grna_hits) %in% rownames(annotations), ]
  
  grna_hits_se <- SingleCellExperiment(assays = list("counts" = grna_hits))
  rowData(grna_hits_se) <- annotations[rownames(grna_hits), ]
  rownames(grna_hits_se) <- gsub("_", "-", rownames(grna_hits_se))
  
  altExps(morris_dataset_1_filtered)[["cre_pert"]] <- grna_hits_se
  
  # remove unexpressed genes
  # morris_dataset_1_filtered <- morris_dataset_1_filtered[rowSums(counts(morris_dataset_1_filtered)) > 20, ]
  
  # estimate means and dispersions
  # morris_dataset_1_filtered <- fit_negbinom_deseq2(morris_dataset_1_filtered, size_factors = "ratio", fit_type = "parametric")
  
  morris_dataset_1_filtered_empty <- morris_dataset_1_filtered
  assays(morris_dataset_1_filtered_empty) <- list()
  
  saveRDS(morris_dataset_1_filtered, paste0("../data/morris_largescreen_processed_", i, ".rds"))
  saveRDS(morris_dataset_1_filtered, paste0("../data/morris_largescreen_processed_empty_", i, ".rds"))
}

## aggregate data

data_list <- lapply(1:4, function(i){readRDS(paste0("../data/morris_largescreen_processed_empty_", i, ".rds"))})
data <- do.call("cbind", data)



saveRDS(data, "../data_morris_largescreen_processed_full_empty.rds")
