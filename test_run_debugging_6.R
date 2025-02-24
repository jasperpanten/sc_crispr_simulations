#!/bin/ Rscript
library(tidyr)
library(sceptre)

# run as
# sbatch 

setwd("/cfs/klemming/projects/supr/lappalainen_lab1/users/panten/projects/sam_simulations/sc_crispr_simulations/")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("./differential_expression_fun.R")
source("./power_simulations_fun.R")

args = commandArgs(trailingOnly=TRUE)

effect_size <- as.numeric(args[[1]])
reps <- as.numeric(args[[2]])
pert_here <- as.character(args[[3]])
repetition <- as.numeric(args[[4]])

# effect_size <- 0.5
# reps <- 10
# pert_here <- "chr1.7428_top_two"
# repetition <- 1

print(effect_size)
print(reps)
print(pert_here)

#data_here_test <- readRDS("../data/sce_gasperini_sam_finished.rds")
#data_here_test <- readRDS("../data/morris_smallscreen_processed_empty.rds")
data_here_test_morris <- readRDS("../data/morris_small_screen_processed_empty.rds")
data_here_test <- readRDS("../data/sce_gasperini_sam_finished_test.rds")
data_here_test <- readRDS("../data/sce_gasperini_sam_finished_empty.rds")

## collapse guides for morris dataset to see if that fixes it: 

alt_exp_here <- altExps(data_here_test_morris)[["cre_pert"]]
row_data_old <- rowData(alt_exp_here)
row_data_old %>% data.frame() %>%
  mutate(gRNAs = cre_target) %>%
  distinct() -> row_data_new
rownames(row_data_new) <- row_data_new$gRNAs

counts_aggregated <- aggregate(counts(alt_exp_here), list(row_data_old$cre_target), max)
rownames(counts_aggregated) <- counts_aggregated[, 1]
counts_aggregated <- counts_aggregated[, -1]

data_here_test_morris_collapsed <- data_here_test_morris
altExps(data_here_test_morris_collapsed)[["cre_pert"]] <- SingleCellExperiment(assays = list("counts" = as.matrix(counts_aggregated[rownames(row_data_new), ])), 
                                                                               colData = colData(alt_exp_here), rowData = row_data_new)

output <- simulate_diff_expr_pooled(sce = data_here_test_morris,
                                    effect_size = effect_size,
                                    pert_level = "cre_pert",
                                    pert_test = NULL,
                                    max_dist = NULL,
                                    genes_iter = F,
                                    guide_sd = 0,
                                    center = FALSE,
                                    rep = reps,
                                    norm = "real",
                                    de_function = de_SCEPTRE_pooled,
                                    formula = ~pert,
                                    n_ctrl = F,
                                    cell_batches = NULL)

output_2 <- simulate_diff_expr_pooled(sce = data_here_test_morris,
                                    effect_size = effect_size,
                                    pert_level = "cre_pert",
                                    pert_test = NULL,
                                    max_dist = NULL,
                                    genes_iter = F,
                                    guide_sd = 0,
                                    center = FALSE,
                                    rep = reps,
                                    norm = "real",
                                    de_function = de_SCEPTRE_pooled,
                                    formula = ~pert,
                                    n_ctrl = F,
                                    cell_batches = NULL)

output_c <- simulate_diff_expr_pooled(sce = data_here_test_morris_collapsed,
                                      effect_size = effect_size,
                                      pert_level = "cre_pert",
                                      pert_test = NULL,
                                      max_dist = NULL,
                                      genes_iter = F,
                                      guide_sd = 0,
                                      center = FALSE,
                                      rep = reps,
                                      norm = "real",
                                      de_function = de_SCEPTRE_pooled,
                                      formula = ~pert,
                                      n_ctrl = F,
                                      cell_batches = NULL)

# saveRDS(output, paste0("../results/morris_small_results/sim_sceptre_res_all_", effect_size, "_", repetition, ".rds"))


### 

# output %>% 
#   distinct() %>%
#   arrange(pvalue)
# 
# output_c %>% 
#   distinct() %>%
#   arrange(pvalue)
# 
# output_2 %>% 
#   distinct() %>%
#   arrange(pvalue)
# 
# output %>% 
#   ungroup() %>%
#   distinct() %>%
#   dplyr::select(c("gene", "cre_pert", "pvalue")) %>%
#   right_join((output_c[ , c("gene", "cre_pert", "pvalue")] %>% ungroup()), by = c("gene", "cre_pert")) -> merged_data
# 
# output %>% 
#   ungroup() %>%
#   distinct() %>%
#   dplyr::select(c("gene", "cre_pert", "pvalue")) %>%
#   right_join((output_2[ , c("gene", "cre_pert", "pvalue")] %>% ungroup()), by = c("gene", "cre_pert")) -> merged_data
# 
# merged_data %>%
#   ggplot(aes(x = -log10(pvalue.x), y = -log10(pvalue.y))) + geom_point() + xlim(c(0, 50)) + ylim(c(0, 50))
# 
# # morris data small: 
# morris_small_data <- readRDS("../data/morris_small_screen_processed.rds")
# 
# per_gene_stats_morris_small <- rowData(morris_small_data)[ , c("mean", "dispersion") ]
# ncells_per_perturbation_morris_small <- rowSums(counts(altExps(morris_small_data)[["cre_pert"]])) %>% 
#   data.frame() %>% rownames_to_column("SNP") %>% mutate(SNP = gsub("-[0-9]$", "", SNP)) %>%
#   group_by(SNP) %>% summarize(ncells = sum(.)) %>%
#   pull(ncells, name = SNP)
# 
# nguides_per_perturbation_morris_small <- rowSums(counts(altExps(morris_small_data)[["cre_pert"]])) %>% 
#   data.frame() %>% rownames_to_column("SNP") %>% mutate(SNP = gsub("-[0-9]$", "", SNP)) %>%
#   group_by(SNP) %>% summarize(nguides = n()) %>%
#   pull(nguides, name = SNP)
# 
# ### 
# 
# cor(-log10(merged_data$pvalue.x), -log10(merged_data$pvalue.y))
# 
# merged_data %>%
#   mutate(ncells = ncells_per_perturbation_morris_small[cre_pert]) %>%
#   mutate(exp_level = per_gene_stats_morris_small[gene, ]$mean) %>%
#   mutate(ncells = cut_number(ncells, n = 10)) %>%
#   mutate(exp_level = cut_number(exp_level, n = 10)) %>%
#   ggplot(aes(x = exp_level, fill = pvalue.x < .1)) + geom_bar(position = "fill")

