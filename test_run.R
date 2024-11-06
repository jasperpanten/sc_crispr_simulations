### Power analysis of perturbation effects (Sams project)
## this is just some tests and playing with our data

library(tidyr)

setwd(paste0(dirname(rstudioapi::getSourceEditorContext()$path)))
# setwd("/cfs/klemming/projects/supr/lappalainen_lab1/users/panten/projects/sam_simulations/sc_crispr_simulations/")

source("./differential_expression_fun.R")
source("./power_simulations_fun.R")
# source("../MKPdiff/scripts/general_purpose/general_purpose_functions.R")

### Outline: 
## Load libraries
## Get data from Gasperini paper
# Get cell number and expression level distributions
# Get estimates of real effect sizes
# Fit dispersion parameters

# parameters
# n_umis_per_gene <- 1
# n_cells_per_perturbation <- 10
# fdr <- .1
# 
# # load data
# data <- readRDS("~/Desktop/dosage_paper/domingo_sce.rds")
# 
# # subset on CRISPRi
# data <- data[ , data$cell_line == "CRISPRi"]
# 
# # remove the Cas9 transcripts
# data <- data[ !grepl("CRISPR", rownames(data)), ]
# 
# # clean a bit... 
# data_cleaned <- SingleCellExperiment(
#   assays = list("counts" = counts(data)),
#   colData = data.frame(
#     gene = ifelse(data$gene == "ntc", "NTC", data$gene),
#     guide = data$guide_call
#   )
# )
# altExps(data_cleaned) <- altExps(data)
# 
# # make perturbation matrix
# pert_matrix <- t(do.call("cbind", lapply(rownames(altExps(data_cleaned)[["GDO"]]), function(x){as.numeric(x == data_cleaned$guide)})))
# colnames(pert_matrix) <- colnames(data_cleaned)
# rownames(pert_matrix) <- rownames(altExps(data_cleaned)[["guide"]])
# gdo_exp <- altExps(data_cleaned)[["GDO"]]
# assays(gdo_exp, withDimnames = F)[["perts"]] <- pert_matrix
# gdo_exp$gene <- data$gene
# altExps(data_cleaned)[["guide"]] <- gdo_exp
# altExps(data_cleaned)[["grna_perts"]] <- gdo_exp
# 
# # add gene coordinates
# library(EnsDb.Hsapiens.v86)
# all_genes <- genes(EnsDb.Hsapiens.v86)
# all_genes <- c(all_genes, GRanges("seqnames" = "None", ranges = IRanges(start = 1, end = 2), gene_name = "NTC"))
# 
# genes_remove <- rownames(data_cleaned)[is.na(match(rownames(data_cleaned), all_genes$gene_name))] # some genes in dataset are not in EnsDb??
# data_cleaned <- data_cleaned[ !rownames(data_cleaned) %in% genes_remove, ]
# genes_in_dataset <- all_genes[ match(rownames(data_cleaned), all_genes$gene_name), ]
# 
# rowData(data_cleaned) <- genes_in_dataset
# rowRanges(data_cleaned) <- genes_in_dataset
# 
# # add gene coordinates to perturbation data
# gene_per_guide <- str_extract(rownames(altExps(data_cleaned)[["guide"]]), "GFI1B|MYB|NFE2|TET2|NTC")
# all_gene_matching <- match(gene_per_guide, all_genes$gene_name)
# rowdata_add <- all_genes[ ifelse(is.na(all_gene_matching), length(all_genes), all_gene_matching), ]
# names(rowdata_add) <- rownames(altExps(data_cleaned)[["guide"]])
# rowdata_add$name <- names(rowdata_add)
# rowRanges(altExps(data_cleaned)[["guide"]]) <- rowdata_add
# rowData(altExps(data_cleaned)[["guide"]]) <- data.frame(rowdata_add) %>% rename("seqnames" = "chr") %>% 
#   mutate(target_name = ifelse(is.na(symbol), NA, paste0(symbol, "_promoter")))
# rowRanges(altExps(data_cleaned)[["grna_perts"]]) <- rowdata_add
# rowData(altExps(data_cleaned)[["grna_perts"]]) <- data.frame(rowdata_add) %>% rename("seqnames" = "chr") %>% 
#   mutate(target_name = ifelse(is.na(symbol), NA, paste0(symbol, "_promoter")))
# 
# data <- data_cleaned
# 
# ## Get realistic statistics for parameters
# 
# # number of guides
# colData(data) %>%
#   data.frame() %>%
#   group_by(guide) %>%
#   summarize(n = n()) %>% 
#   {
#     ggplot(data = ., aes(x = n)) + geom_histogram(fill = "white", col = "black") + theme_paper() + 
#       xlab("Number of cells") + ylab("Number of guides") + 
#       geom_vline(xintercept = median(.$n), col = "red") + 
#       annotate(x = 200, y = 15, geom = "text", label = paste0("Median number of guides: \n", median(.$n)), col = "red", size = 8)
#   }
#   
# # expression level
# counts(data) %>%
#   rowMeans() %>%
#   data.frame() %>%
#   rownames_to_column("Gene") %>% 
#   {
#     ggplot(data = ., aes(x = .)) + geom_histogram(fill = "white", col = "black") + theme_paper() + 
#       xlab("Number of reads") + ylab("Number of genes") + 
#       geom_vline(xintercept = median(.$.), col = "red") + 
#       annotate(x = 0.1, y = 10, geom = "text", label = paste0("Median number of reads \n per gene: ", round(median(.$.), digits = 2)), col = "red", size = 8) + 
#       scale_x_log10()
#   }

# test for DE for each perturbation

# first filter cells for minimum and maximum number total UMIs per cell
# sce <- filter_umis_per_cell(data, min_umis = n_umis_per_gene, max_umis = Inf, recompute = T)

# filter for minimum number of cells per perturbation
# sce <- filter_cells_per_pert(data, min_cells = n_cells_per_perturbation, pert_level = "guide", altExp = F)

# normalize data
# data <- normalize_cens_mean(data)
# assay(data, "logcounts") <- log1p(assay(data, "normcounts"))

# output <- test_differential_expression(data,
#                                        pert_level = "guide",
#                                        max_dist = NULL,
#                                        de_function = de_MAST,
#                                        formula = ~pert,
#                                        n_ctrl = F,
#                                        p_adj_method = "bonferroni")

# analyze perturbation effects... 
# output %>%
#   data.frame() %>%
#   ggplot(aes(x = logFC, y = -log10(pvalue))) + geom_point(size = .1) + theme_paper()


# fit dispersion parameters
# library(DESeq2)
# data <- fit_negbinom_deseq2(data, size_factors = "ratio", fit_type = "parametric")

# run simulations
# effect_size <- 2.0
# 
# guides_check <- c("GFI1B.chr9:132985047-132985066.+.WT", "NTC-1", "NTC-2", "NTC-3", "NTC-4", "NTC-5")
# data_subset <- data[ , data$guide %in% guides_check]
# altExps(data_subset)[["guide"]] <- altExps(data_subset)[["guide"]][guides_check, ]
# altExps(data_subset)[["grna_perts"]] <- altExps(data_subset)[["grna_perts"]][guides_check, ]
# 
# output <- simulate_diff_expr(data_subset,
#                              effect_size = effect_size,
#                              pert_level = "grna_perts",
#                              max_dist = 1e6,
#                              genes_iter = FALSE,
#                              guide_sd = 0,
#                              center = FALSE,
#                              rep = 10,
#                              norm = "real",
#                              de_function = de_MAST,
#                              formula = ~pert,
#                              n_ctrl = F,
#                              cell_batches = NULL)

# saveRDS(output, "~/Desktop/output_test.rds")

# pval_cutoff <- output %>% dplyr::filter(iteration == 1) %>% mutate(p_adj = p.adjust(pvalue, method = "BH")) %>%
#   dplyr::filter(p_adj < .1) %>% arrange(-p_adj) %>% pull(pvalue)
# pval_cutoff <- pval_cutoff[[1]]
# output %>% 
#   group_by(perturbation, gene) %>%
#   summarize(n_nonsig = sum(pvalue >= .01), n_sig = sum(pvalue < pval_cutoff)) %>%
#   mutate(power = n_sig / (n_sig + n_nonsig)) %>%
#   ggplot(aes(x = power)) + geom_histogram()
# 
# output %>% 
#   ggplot(aes(x = logFC)) + geom_histogram()
# 
# hist(log2(exp(output$logFC)), breaks = 100)

### for sams data: procedure: 
# implement test with sceptre
# question: what is the power for all tested CRE-gene links to discover an effect at different effect sizes?
# - need actual single-cell data
# - for each tested CRE-gene pair: 
# --- what is the power to see an effect for b %in% c(0.01, 0.1, 0.5, 1, 5, ..., whatever realistic effect sizes are ), given real n_cells, mu, disp
# --- for each CRE, do we see an effect stratified by power, how many effects are we likely missing
# --- show power landscape as heatmap against n_cells / expression and for different b
# --- relationship between power and eqtl being identified (talk to sam)

## morris
# sting_seq_data <- read.csv("~/Desktop/PostDoc_TL_Lab/Projects/Sam/All_STING_seq_CREs.csv", check.names = F)
# 
# pval_cutoff <- sting_seq_data %>%
#   dplyr::filter(`Q-value (1 Mb)` < .1) %>%
#   arrange(- `Skew fit p-value`) %>%
#   pull(`Skew fit p-value`)
# 
# # volcano plot of differnetial effects
# sting_seq_data %>% 
#   ggplot(aes(x = `Log2 fold-change`, y = -log10(`Skew fit p-value`))) + geom_point() + 
#     theme_paper() + geom_hline(yintercept = -log10(pval_cutoff[[1]]), linetype = 'dashed')
# 
# median_effect_size <- 
#   sting_seq_data %>%
#   dplyr::filter(`Q-value (1 Mb)` < .1) %>% 
#   pull(`Log2 fold-change`) %>% median()
# sting_seq_data %>%
#   dplyr::filter(`Q-value (1 Mb)` < .1) %>% 
#   mutate(`Log2 fold-change` = ifelse(`Log2 fold-change` < -2, -2, `Log2 fold-change`)) %>% 
#   {
#   ggplot(data = ., aes(x = `Log2 fold-change`)) + geom_histogram(fill = "white", col = "black") + theme_paper() + 
#     ylab("Number of effects") + geom_vline(xintercept = median_effect_size, col = "red", linetype = 'dashed') 
#   }

## gasperini
# gasperini_data <- fst::read.fst("~/Desktop/PostDoc_TL_Lab/Projects/Sam/resampling_results.fst")
# 
# pval_cutoff <- gasperini_data %>%
#   dplyr::filter(rejected) %>%
#   arrange(- p_value) %>%
#   pull(p_value)
# 
# # volcano plot of differnetial effects
# gasperini_data %>% 
#   ggplot(aes(x = xi, y = -log10(p_value))) + geom_point(size = .1) + 
#   theme_paper() + geom_hline(yintercept = -log10(pval_cutoff[[1]]), linetype = 'dashed')
# 
# median_effect_size <- gasperini_data %>%
#   dplyr::filter(rejected) %>%
#   pull(xi) %>% median()
# gasperini_data %>%
#   dplyr::filter(rejected) %>%
#   mutate(xi = ifelse(xi < -4, -4, xi)) %>% {
#     ggplot(data = ., aes(x = xi)) + geom_histogram(fill = "white", col = "black") + theme_paper() + 
#       ylab("Number of effects") + geom_vline(xintercept = median_effect_size, col = "red", linetype = 'dashed') 
#   }


# test sceptre for DE instead

# BiocManager::install("sceptre")

# tutorial data
# library(sceptre)
# library(sceptredata)
# data(highmoi_example_data)
# data(grna_target_data_frame_highmoi)
# # import data
# sceptre_object <- import_data(
#   response_matrix = highmoi_example_data$response_matrix,
#   grna_matrix = highmoi_example_data$grna_matrix,
#   grna_target_data_frame = grna_target_data_frame_highmoi,
#   moi = "high",
#   extra_covariates = highmoi_example_data$extra_covariates,
#   response_names = highmoi_example_data$gene_names
# )
# # set analysis parameters, assign grnas, run qc
# discovery_pairs <- construct_cis_pairs(sceptre_object)
# sceptre_object <- sceptre_object |>
#   set_analysis_parameters(
#     side = "left",
#     resampling_mechanism = "permutations",
#     discovery_pairs = discovery_pairs
#   ) |>
#   assign_grnas(method = "thresholding") |>
#   run_qc() |>
#   run_discovery_analysis(
#     parallel = F,
#     n_processors = 2
#   )

# results from sam: 
# perturbations_sam <- read_csv("../data/gasperini_target_sites_Jasper_long.txt")[[1]]
# gasperini_data <- fst::read.fst("../data/resampling_results.fst") %>% 
#   dplyr::filter(target_site %in% perturbations_sam)
# 
# gasperini_data %>%
#   group_by(target_site) %>% 
#   summarize(n = n())
# 
# gasperini_data %>%
#   ggplot(aes(x = rejected, y = xi)) + geom_violin() + ylim(c(-2.5, 2.5))

# perform power analysis with sceptre instead of mast
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(sceptre)

#data_here <- readRDS("../data/sce_gasperini_sam.rds")
# data_here <- readRDS("../data/sce_gasperini_subset.rds")

#perturbation_test <- c("chr1.7428_top_two", "chr1.9538_top_two")
#perturbation_test_collapsed <- paste0(perturbation_test, collapse = "|")
#cells_here <- grepl(perturbation_test_collapsed, data_here$gene)
#data_here_test <- data_here[ , cells_here]
#altExps(data_here_test)[["cre_pert"]] <- altExps(data_here)[["cre_pert"]][perturbation_test , cells_here]

#data_here_test <- data_here
#data_here_test <- logNormCounts(data_here_test)
#data_here_test <- fit_negbinom_deseq2(data_here_test)

#saveRDS(data_here_test, "../data/sce_gasperini_sam_dispersions_test.rds")
data_here_test <- readRDS("../data/sce_gasperini_sam_dispersions.rds")

genes(EnsDb.Hsapiens.v86) %>% data.frame() %>% dplyr::select("seqnames", "start", "end", "gene_id") %>% rename("gene_id" = "id") -> gene_coordinates
data_here_test <- data_here_test[rownames(data_here_test) %in% gene_coordinates$id, ]

row_data_temp <- rowData(data_here_test) %>% data.frame() %>% left_join(gene_coordinates) %>% column_to_rownames("id")
rowRanges(data_here_test) <- makeGRangesFromDataFrame(data.frame(seqnames = row_data_temp$seqnames, start = row_data_temp$start, end = row_data_temp$end, row.names = rownames(row_data_temp)))
rowData(data_here_test) <- row_data_temp

saveRDS(data_here_test, "../data/sce_gasperini_sam_finished.rds")

data_here_test <- readRDS("../data/sce_gasperini_sam_finished.rds")

output <- simulate_diff_expr(sce = data_here_test,
                             effect_size = .5,
                             pert_level = "cre_pert",
                             pert_test = "chr1.7428_top_two",
                             max_dist = NULL,
                             genes_iter = F,
                             guide_sd = 0,
                             center = FALSE,
                             rep = 10,
                             norm = "real",
                             de_function = de_SCEPTRE,
                             formula = ~pert,
                             n_ctrl = F,
                             cell_batches = NULL)

# saveRDS(output, "../results/simulation_output.rds")

# output %>%
#   group_by(gene, perturbation) %>%
#   mutate(padj = pvalue * 20000 * 20) %>%
#   summarize(
#     sig = sum(padj < .1),
#     nonsig = sum(padj >= .1)
#   ) %>%
#   mutate(fraction_sig = sig / (sig + nonsig)) %>%
#   ggplot(aes(x = sig)) + geom_histogram()
