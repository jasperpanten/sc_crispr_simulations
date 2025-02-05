# get statistics on gene expression and number of cells per perturbation in each dataset: 

setwd(paste0(dirname(rstudioapi::getSourceEditorContext()$path)))

# gasperini data:
gasperini_data <- readRDS("../data/sce_gasperini_sam_finished_empty.rds")

per_gene_stats_gasperini <- rowData(gasperini_data)[ , c("mean", "dispersion") ]
ncells_per_perturbation_gasperini <- rowSums(counts(altExps(gasperini_data)[["cre_pert"]]))


# morris data small: 
morris_small_data <- readRDS("../data/morris_small_screen_processed.rds")

per_gene_stats_morris_small <- rowData(morris_small_data)[ , c("mean", "dispersion") ]
ncells_per_perturbation_morris_small <- rowSums(counts(altExps(morris_small_data)[["cre_pert"]])) %>% 
  data.frame() %>% rownames_to_column("SNP") %>% mutate(SNP = gsub("-[0-9]$", "", SNP)) %>%
  group_by(SNP) %>% summarize(ncells = sum(.)) %>%
  pull(ncells, name = SNP)


# morris data large: 
morris_large_data <- readRDS("../data/morris_largescreen_processed_full_empty.rds")

per_gene_stats_morris_large <- rowData(morris_large_data)[ , c("mean", "dispersion") ]
ncells_per_perturbation_morris_large <- rowSums(counts(altExps(morris_large_data)[["cre_pert"]])) %>% 
  data.frame() %>% rownames_to_column("SNP") %>% mutate(SNP = gsub("-[0-9]$", "", SNP)) %>%
  group_by(SNP) %>% summarize(ncells = sum(.)) %>%
  pull(ncells, name = SNP)


# process and plot results of power analysis 
output_list <- lapply(list.files("~/Desktop/gasperini_results/", full.name = T), function(x){
  data = readRDS(x)
  data$effect_size <- as.numeric(gsub("sim_sceptre_res_all_|_[0-9]*.rds", "", basename(x)))
  data
}) %>%
  do.call("rbind", .)

output_list %>%
  group_by(gene, cre_pert) %>%
  # mutate(padj = pvalue * 20000 * 20) %>%
  mutate(padj = pvalue) %>%
  summarize(
    sig = sum(padj < .1),
    nonsig = sum(padj >= .1)
  ) %>%
  mutate(fraction_sig = sig / (sig + nonsig)) %>%
  ggplot(aes(x = fraction_sig)) + geom_histogram()

# read output 

output_new <- output_list

n_tests <- output_new %>%
  group_by(gene, cre_pert) %>% 
  summarize(n = n()) %>% nrow()

testy <- output_new %>%
  group_by(iteration, effect_size) %>%
  mutate(padj = p.adjust(pvalue, method = 'BH'))

output_new %>%
  group_by(effect_size) %>%
  dplyr::filter(!is.infinite(logFC)) %>%
  summarize(logFC = mean(logFC, na.rm = T)) %>%
  mutate(effect_size_measured = 2 ** logFC)

# first plot without pvalue adjustment: 
output_new %>%
  dplyr::filter(!is.na(pvalue)) %>%
  dplyr::filter(!is.na(effect_size)) %>%
  #mutate(padj = pvalue) %>%
  ## mutate(padj = p.adjust(pvalue)) %>%
  group_by(iteration, effect_size) %>%
  mutate(padj = pvalue) %>%
  # mutate(padj = pvalue) %>%
  group_by(gene, cre_pert, effect_size) %>%
  summarize(
    sig = sum(padj < .1),
    nonsig = sum(padj >= .1)
  ) %>% 
  mutate(fraction_sig = sig / (sig + nonsig)) %>% 
  ungroup() -> output_processed

output_processed %>%
  group_by(effect_size, fraction_sig) %>%
  summarize(n_hits = n()) %>% 
  ungroup() %>%
  complete(fraction_sig, effect_size, fill = list(n_hits = 0)) %>%
  group_by(effect_size) %>%
  arrange(-fraction_sig) %>%
  mutate(n_total = cumsum(n_hits)) %>%
  ungroup() %>%
  mutate(fraction_total = n_total / max(n_total)) %>%
  add_row(effect_size = unique(.$effect_size), fraction_sig = 1, fraction_total = 0) %>%
  ggplot(aes(x = fraction_total, y = fraction_sig, col = as.factor(effect_size))) + geom_step(linewidth = 2) + xlim(c(0, 1)) + ylim(c(0, 1)) + 
  geom_hline(yintercept = 0.8, linetype = 'dashed') + theme_classic(base_size = 15) + 
  xlab("Proportion of CRE-gene pairs") + ylab("Power") + labs("color" = "Effect size") + 
  ggtitle("CRISPRi effect detection power for tested \n CRE−gene links (SCEPTRE, Gasperini data)")
ggsave("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/plots/241216_gasperini_5it_no_adjust.pdf")

output_new %>%
  dplyr::filter(!is.na(pvalue)) %>%
  dplyr::filter(!is.na(effect_size)) %>%
  #mutate(padj = pvalue) %>%
  ## mutate(padj = p.adjust(pvalue)) %>%
  group_by(iteration, effect_size) %>%
  mutate(padj = p.adjust(pvalue, method = 'BH')) %>%
  # mutate(padj = pvalue) %>%
  group_by(gene, cre_pert, effect_size) %>%
  summarize(
    sig = sum(padj < .1),
    nonsig = sum(padj >= .1)
  ) %>% 
  mutate(fraction_sig = sig / (sig + nonsig)) %>% 
  ungroup() -> output_processed

output_processed %>%
  group_by(effect_size, fraction_sig) %>%
  summarize(n_hits = n()) %>% 
  ungroup() %>%
  complete(fraction_sig, effect_size, fill = list(n_hits = 0)) %>%
  group_by(effect_size) %>%
  arrange(-fraction_sig) %>%
  mutate(n_total = cumsum(n_hits)) %>%
  ungroup() %>%
  mutate(fraction_total = n_total / max(n_total)) %>%
  add_row(effect_size = unique(.$effect_size), fraction_sig = ifelse(unique(.$effect_size) == 0.9, 0, 1), fraction_total = 0) %>%
  ggplot(aes(x = fraction_total, y = fraction_sig, col = as.factor(effect_size))) + geom_step(linewidth = 2) + xlim(c(0, 1)) + ylim(c(0, 1)) + 
  geom_hline(yintercept = 0.8, linetype = 'dashed') + theme_classic(base_size = 15) + 
  xlab("Proportion of CRE-gene pairs") + ylab("Power") + labs("color" = "Effect size") + 
  ggtitle("CRISPRi effect detection power for tested \n CRE−gene links (SCEPTRE, Gasperini data)")
ggsave("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/plots/241216_gasperini_5it_benjaminihochberg.pdf")

output_processed %>% write_csv("../results/gasperini_power_per_pair.csv")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### morris small analysis
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

output_list <- lapply(list.files("~/Desktop/morris_small_results/", full.name = T), function(x){
  data = readRDS(x)
  data$effect_size <- as.numeric(gsub("sim_sceptre_res_all_|_[0-9]*.rds", "", basename(x)))
  data
}) %>%
  do.call("rbind", .) %>%
  distinct()

output_list %>%
  group_by(gene, cre_pert) %>%
  # mutate(padj = pvalue * 20000 * 20) %>%
  mutate(padj = pvalue) %>%
  summarize(
    sig = sum(padj < .1),
    nonsig = sum(padj >= .1)
  ) %>%
  mutate(fraction_sig = sig / (sig + nonsig)) %>%
  ggplot(aes(x = sig)) + geom_histogram()

# read output 

output_new <- output_list

n_tests <- output_new %>%
  group_by(gene, cre_pert) %>% 
  summarize(n = n()) %>% nrow()

testy <- output_new %>%
  group_by(iteration, effect_size) %>%
  mutate(padj = p.adjust(pvalue, method = 'BH'))

output_new %>%
  group_by(effect_size) %>%
  dplyr::filter(!is.infinite(logFC)) %>%
  summarize(logFC = mean(logFC, na.rm = T)) %>%
  mutate(effect_size_measured = 2 ** logFC)

# first plot without pvalue adjustment: 
output_new %>%
  dplyr::filter(!is.na(pvalue)) %>%
  #mutate(padj = pvalue) %>%
  ## mutate(padj = p.adjust(pvalue)) %>%
  group_by(iteration, effect_size) %>%
  mutate(padj = pvalue) %>%
  # mutate(padj = pvalue) %>%
  group_by(gene, cre_pert, effect_size) %>%
  summarize(
    sig = sum(padj < .1),
    nonsig = sum(padj >= .1)
  ) %>% 
  mutate(fraction_sig = sig / (sig + nonsig)) %>% 
  ungroup() -> output_processed

output_processed %>%
  group_by(effect_size, fraction_sig) %>%
  summarize(n_hits = n()) %>% 
  ungroup() %>%
  complete(fraction_sig, effect_size, fill = list(n_hits = 0)) %>%
  group_by(effect_size) %>%
  arrange(-fraction_sig) %>%
  mutate(n_total = cumsum(n_hits)) %>%
  ungroup() %>%
  mutate(fraction_total = n_total / max(n_total)) %>%
  add_row(effect_size = unique(.$effect_size), fraction_sig = ifelse(unique(.$effect_size) == 0.9, 1, 1), fraction_total = 0) %>%
  ggplot(aes(x = fraction_total, y = fraction_sig, col = as.factor(effect_size))) + geom_step(linewidth = 2) + xlim(c(0, 1)) + ylim(c(0, 1)) + 
  geom_hline(yintercept = 0.8, linetype = 'dashed') + theme_classic(base_size = 15) + 
  xlab("Proportion of CRE-gene pairs") + ylab("Power") + labs("color" = "Effect size") + 
  ggtitle("CRISPRi effect detection power for tested \n CRE−gene links (SCEPTRE, Morris small data)")
ggsave("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/plots/241222_morris_small_5it_no_adjust.pdf")

output_new %>%
  dplyr::filter(!is.na(pvalue)) %>%
  #mutate(padj = pvalue) %>%
  ## mutate(padj = p.adjust(pvalue)) %>%
  group_by(iteration, effect_size) %>%
  mutate(padj = p.adjust(pvalue, method = 'BH')) %>%
  # mutate(padj = pvalue) %>%
  group_by(gene, cre_pert, effect_size) %>%
  summarize(
    sig = sum(padj < .1),
    nonsig = sum(padj >= .1)
  ) %>% 
  mutate(fraction_sig = sig / (sig + nonsig)) %>% 
  ungroup() -> output_processed

output_processed %>%
  group_by(effect_size, fraction_sig) %>%
  summarize(n_hits = n()) %>% 
  ungroup() %>%
  complete(fraction_sig, effect_size, fill = list(n_hits = 0)) %>%
  group_by(effect_size) %>%
  arrange(-fraction_sig) %>%
  mutate(n_total = cumsum(n_hits)) %>%
  ungroup() %>%
  mutate(fraction_total = n_total / max(n_total)) %>%
  add_row(effect_size = unique(.$effect_size), fraction_sig = ifelse(unique(.$effect_size) == 0.9, 1, 1), fraction_total = 0) %>%
  ggplot(aes(x = fraction_total, y = fraction_sig, col = as.factor(effect_size))) + geom_step(linewidth = 2) + xlim(c(0, 1)) + ylim(c(0, 1)) + 
  geom_hline(yintercept = 0.8, linetype = 'dashed') + theme_classic(base_size = 15) + 
  xlab("Proportion of CRE-gene pairs") + ylab("Power") + labs("color" = "Effect size") + 
  ggtitle("CRISPRi effect detection power for tested \n CRE−gene links (SCEPTRE, Morris small data)")
ggsave("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/plots/241222_morris_small_5it_benjaminihochberg.pdf")

output_processed %>% write_csv("../results/morris_small_power_per_pair.csv")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### morris large analysis
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

output_list <- lapply(list.files("~/Desktop/morris_large_results//", full.name = T), function(x){
  data = readRDS(x)
  data$effect_size <- as.numeric(gsub("sim_sceptre_res_all_|_[0-9]*.rds", "", basename(x)))
  data
}) %>%
  do.call("rbind", .) %>%
  distinct()

output_list %>%
  group_by(gene, cre_pert) %>%
  # mutate(padj = pvalue * 20000 * 20) %>%
  mutate(padj = pvalue) %>%
  summarize(
    sig = sum(padj < .1),
    nonsig = sum(padj >= .1)
  ) %>%
  mutate(fraction_sig = sig / (sig + nonsig)) %>%
  ggplot(aes(x = sig)) + geom_histogram()

# read output 

output_new <- output_list

n_tests <- output_new %>%
  group_by(gene, cre_pert) %>% 
  summarize(n = n()) %>% nrow()

output_new %>%
  group_by(effect_size) %>%
  dplyr::filter(!is.infinite(logFC)) %>%
  summarize(logFC = mean(logFC, na.rm = T)) %>%
  mutate(effect_size_measured = 2 ** logFC)

# first plot without pvalue adjustment: 
output_new %>%
  dplyr::filter(!is.na(effect_size)) %>%
  dplyr::filter(!is.na(pvalue)) %>%
  #mutate(padj = pvalue) %>%
  ## mutate(padj = p.adjust(pvalue)) %>%
  group_by(iteration, effect_size) %>%
  mutate(padj = pvalue) %>%
  # mutate(padj = pvalue) %>%
  group_by(gene, cre_pert, effect_size) %>%
  summarize(
    sig = sum(padj < .1),
    nonsig = sum(padj >= .1)
  ) %>% 
  mutate(fraction_sig = sig / (sig + nonsig)) %>% 
  ungroup() -> output_processed

output_processed %>%
  group_by(effect_size, fraction_sig) %>%
  summarize(n_hits = n()) %>% 
  ungroup() %>%
  complete(fraction_sig, effect_size, fill = list(n_hits = 0)) %>%
  group_by(effect_size) %>%
  arrange(-fraction_sig) %>%
  mutate(n_total = cumsum(n_hits)) %>%
  ungroup() %>%
  mutate(fraction_total = n_total / max(n_total)) %>%
  add_row(effect_size = unique(.$effect_size), fraction_sig = 1, fraction_total = 0) %>%
  ggplot(aes(x = fraction_total, y = fraction_sig, col = as.factor(effect_size))) + geom_step(linewidth = 2) + xlim(c(0, 1)) + ylim(c(0, 1)) + 
  geom_hline(yintercept = 0.8, linetype = 'dashed') + theme_classic(base_size = 15) + 
  xlab("Proportion of CRE-gene pairs") + ylab("Power") + labs("color" = "Effect size") + 
  ggtitle("CRISPRi effect detection power for tested \n CRE−gene links (SCEPTRE, Morris data)")
ggsave("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/plots/241222_morris_large_5it_no_adjust.pdf")

output_new %>%
  dplyr::filter(!is.na(effect_size)) %>%
  dplyr::filter(!is.na(pvalue)) %>%
  #mutate(padj = pvalue) %>%
  ## mutate(padj = p.adjust(pvalue)) %>%
  group_by(iteration, effect_size) %>%
  mutate(padj = p.adjust(pvalue, method = 'BH')) %>%
  # mutate(padj = pvalue) %>%
  group_by(gene, cre_pert, effect_size) %>%
  summarize(
    sig = sum(padj < .1),
    nonsig = sum(padj >= .1)
  ) %>% 
  mutate(fraction_sig = sig / (sig + nonsig)) %>% 
  ungroup() -> output_processed

output_processed %>%
  group_by(effect_size, fraction_sig) %>%
  summarize(n_hits = n()) %>% 
  ungroup() %>%
  complete(fraction_sig, effect_size, fill = list(n_hits = 0)) %>%
  group_by(effect_size) %>%
  arrange(-fraction_sig) %>%
  mutate(n_total = cumsum(n_hits)) %>%
  ungroup() %>%
  mutate(fraction_total = n_total / max(n_total)) %>%
  add_row(effect_size = unique(.$effect_size), fraction_sig = ifelse(unique(.$effect_size) %in% c(0.8, 0.9), 0, 1), fraction_total = 0) %>%
  ggplot(aes(x = fraction_total, y = fraction_sig, col = as.factor(effect_size))) + geom_step(linewidth = 2) + xlim(c(0, 1)) + ylim(c(0, 1)) + 
  geom_hline(yintercept = 0.8, linetype = 'dashed') + theme_classic(base_size = 15) + 
  xlab("Proportion of CRE-gene pairs") + ylab("Power") + labs("color" = "Effect size") + 
  ggtitle("CRISPRi effect detection power for tested \n CRE−gene links (SCEPTRE, Morris data)")
ggsave("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/plots/241222_morris_large_5it_benjaminihochberg.pdf")

output_processed %>% write_csv("../results/morris_large_power_per_pair.csv")

### read per perturbation results
output_processed_gasperini <- read.csv("../results/gasperini_power_per_pair.csv") %>% mutate(dataset = "Gasperini")
output_processed_morris_small <- read.csv("../results/morris_small_power_per_pair.csv") %>% mutate(dataset = "Morris_small")
output_processed_morris_large <- read.csv("../results/morris_large_power_per_pair.csv") %>% mutate(dataset = "Morris_large")

# add expression level and number of perturbed cells
output_processed_gasperini <- output_processed_gasperini %>% 
  mutate(expression_level = per_gene_stats_gasperini[gene, ]$mean) %>%
  mutate(expression_dispersion = per_gene_stats_gasperini[gene, ]$dispersion) %>%
  mutate(number_of_cells = ncells_per_perturbation_gasperini[cre_pert])

output_processed_morris_small <- output_processed_morris_small %>% 
  mutate(expression_level = per_gene_stats_morris_small[gene, ]$mean) %>%
  mutate(expression_dispersion = per_gene_stats_morris_small[gene, ]$dispersion) %>%
  mutate(number_of_cells = ncells_per_perturbation_morris_small[cre_pert])

output_processed_morris_large <- output_processed_morris_large %>% 
  mutate(expression_level = per_gene_stats_morris_large[gene, ]$mean) %>%
  mutate(expression_dispersion = per_gene_stats_morris_large[gene, ]$dispersion) %>%
  mutate(number_of_cells = ncells_per_perturbation_morris_large[cre_pert])

output_processed <- rbind(output_processed_gasperini, output_processed_morris_small, output_processed_morris_large)

output_processed %>% 
  dplyr::filter(!is.na(number_of_cells)) %>%
  dplyr::filter(effect_size == 0.5) %>%
  sample_frac(1) %>%
  ggplot(aes(x = expression_level, y = number_of_cells, col = fraction_sig >=.8)) + 
    geom_point(size = 1) + facet_wrap(~dataset) + scale_x_log10() + xlab("Expression Level") + ylab("Number of cells")
ggsave("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/plots/250110_power_expressionlevel_cellnumber.pdf")

output_processed %>% 
  dplyr::filter(!is.na(number_of_cells)) %>%
  dplyr::filter(effect_size == 0.5) %>%
  ggplot(aes(x = factor(round(fraction_sig, digits = 1)), y = expression_level)) + geom_boxplot() + facet_wrap(~dataset) + scale_y_log10() + 
  xlab("Power") + ylab("Expression level")
ggsave("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/plots/250110_power_box_expression_level.pdf")

output_processed %>% 
  dplyr::filter(!is.na(number_of_cells)) %>%
  dplyr::filter(effect_size == 0.5) %>%
  ggplot(aes(x = factor(round(fraction_sig, digits = 1)), y = number_of_cells)) + geom_boxplot() + facet_wrap(~dataset) + 
  xlab("Power") + ylab("Number of cells")
ggsave("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/plots/250110_power_box_cellnumber.pdf")

