library(scater)
library(scran)
library(tidyverse)
library(DropletUtils)

setwd(paste0(dirname(rstudioapi::getSourceEditorContext()$path)))

get_quantiles <- function(values, qlow = .1, qhigh = .9){
  return(c(quantile(values, qlow), quantile(values, qhigh)))
}

morris_dataset_1 <- readRDS("../data/morris_raw_small_screen.rds")

# calculate qc metrics and set filters
mito_genes <- rowData(morris_dataset_1) %>% data.frame() %>% dplyr::filter(grepl("^MT-", Symbol)) %>% pull(ID)
colData(morris_dataset_1) <- cbind(colData(morris_dataset_1), perCellQCMetrics(morris_dataset_1, subset = list("mito" = mito_genes)))

quantiles_sum <- get_quantiles(morris_dataset_1$sum, 0.20, 0.99)
qplot(x = "1", y = morris_dataset_1$sum) + geom_violin() + geom_hline(yintercept = quantiles_sum, col = "red")

quantiles_detected <- get_quantiles(morris_dataset_1$detected, 0.15, 0.99)
qplot(x = "1", y = morris_dataset_1$detected) + geom_violin() + geom_hline(yintercept = quantiles_detected, col = "red")

quantiles_mito <- get_quantiles(morris_dataset_1$subsets_mito_percent, 0.05, 0.90)
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

# reimplementation of HTODemux for scater objects
HTODemux_scater <- function(
    object,
    assay = "hto",
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
  counts <- counts(data)

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
        x = t(x = counts(data)),
        centers = ncenters,
        nstart = nstarts
      )
      #identify positive and negative signals for all HTO
      Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
    },
    'clara' = {
      #use fast k-medoid clustering
      init.clusters <- clara(
        x = t(x = counts(data)),
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
  average.expression <- counts(aggregateAcrossCells(object, object$clustering))
  
  #checking for any cluster with all zero counts for any barcode
  if (any(colSums(average.expression == 0) > 0)) {
    stop("Cells with zero counts exist as a cluster.")
  }
  #create a matrix to store classification result
  discrete <- counts(object)
  discrete[discrete > 0] <- 0
  # for each HTO, we will use the minimum cluster for fitting
  for (iter in rownames(x = data)) {
    values <- counts[iter, colnames(object)]
    #commented out if we take all but the top cluster as background
    #values_negative=values[setdiff(object@cell.names,WhichCells(object,which.max(average.expression[iter,])))]
    values.use <- values[WhichCells(
      object = object,
      idents = levels(x = Idents(object = object))[[which.min(x = average.expression[iter, ])]]
    )]
    fit <- suppressWarnings(expr = fitdist(data = values.use, distr = "nbinom"))
    cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
    discrete[iter, names(x = which(x = values > cutoff))] <- 1
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
  hash.max <- apply(X = data, MARGIN = 2, FUN = max)
  hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
  hash.second <- apply(X = data, MARGIN = 2, FUN = MaxN, N = 2)
  hash.maxID <- as.character(x = donor.id[sapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(which(x = data[, x] == hash.max[x])[1])
    }
  )])
  hash.secondID <- as.character(x = donor.id[sapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(which(x = data[, x] == hash.second[x])[1])
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
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste0(assay, '_classification')
  # Idents(object, cells = rownames(object@meta.data[object@meta.data$classification.global == "Doublet", ])) <- "Doublet"
  doublets <- rownames(x = object[[]])[which(object[[paste0(assay, "_classification.global")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- 'Doublet'
  # object@meta.data$hash.ID <- Idents(object)
  object$hash.ID <- Idents(object = object)
  return(object)
}

counts(altExps(morris_dataset_1)[["hto"]]) %>%
  t() %>%
  as.matrix() %>%
  data.frame() %>%
  pivot_longer(-c(HTO21)) %>%
  ggplot(aes(x = HTO21, y = value)) + geom_point(size = .1) + facet_wrap(~name) + scale_x_log10() + scale_y_log10()

logcounts(altExps(morris_dataset_1)[["hto"]]) %>%
  t() %>%
  as.matrix() %>%
  data.frame() %>%
  pivot_longer(-c()) %>%
  ggplot(aes(x = value)) + facet_wrap(~name) + geom_density() + scale_x_log10()


