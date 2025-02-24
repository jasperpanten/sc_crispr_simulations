#!/bin/ Rscript
library(tidyr)
library(sceptre)

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
# reps <- 1
# pert_here <- "chr1.7428_top_two"
# repetition <- 1

print(effect_size)
print(reps)
print(pert_here)

print("gets here 1")

# data_here_test <- readRDS("../data/sce_gasperini_sam_finished.rds")
data_here_test <- readRDS("../data/sce_gasperini_sam_finished_test.rds")

print("gets here 2")

# split across perturbations: 
output <- simulate_diff_expr(data_here_test,
                             effect_size = effect_size,
                             pert_level = "cre_pert",
                             pert_test = NULL,
                             max_dist = NULL,
                             genes_iter = F,
                             guide_sd = 0,
                             center = FALSE,
                             rep = reps,
                             norm = "real",
                             de_function = de_DESeq,
                             formula = ~pert,
                             n_ctrl = F,
                             cell_batches = NULL)

#saveRDS(output, paste0("../results/gasperini_results/sim_res_", pert_here, "_", effect_size, "_", repetition, ".rds"))
