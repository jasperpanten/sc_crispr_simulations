library(tidyr)
library(sceptre)

setwd("~/Desktop/PostDoc_TL_Lab/Projects/Sam/sc_crispr_simulations/")

source("./differential_expression_fun.R")
source("./power_simulations_fun.R")

effect_size <- .5
reps <- 10

data_here_test <- readRDS("../data/sce_gasperini_sam_finished.rds")

# split across perturbation../data/GSE120861_50k_reference_cells.rds# split across perturbations: 

output <- simulate_diff_expr(data_here_test,
                             effect_size = effect_size,
                             pert_level = "cre_pert",
                             max_dist = NULL,
                             genes_iter = F,
                             guide_sd = 0,
                             center = FALSE,
                             rep = reps,
                             norm = "real",
                             de_function = de_SCEPTRE,
                             formula = ~pert,
                             n_ctrl = F,
                             cell_batches = NULL)

#saveRDS(output, "../results/simulation_output_test.rds")
