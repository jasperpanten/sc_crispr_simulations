# get statistics on gene expression and number of cells per perturbation in each dataset: 

library(tidyverse)

theme_paper <- function(textsize = 20){theme_bw(base_size = textsize)}

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

nguides_per_perturbation_morris_small <- rowSums(counts(altExps(morris_small_data)[["cre_pert"]])) %>% 
  data.frame() %>% rownames_to_column("SNP") %>% mutate(SNP = gsub("-[0-9]$", "", SNP)) %>%
  group_by(SNP) %>% summarize(nguides = n()) %>%
  pull(nguides, name = SNP)


# morris data large: 
morris_large_data <- readRDS("../data/morris_largescreen_processed_full_empty_new.rds")

per_gene_stats_morris_large <- rowData(morris_large_data)[ , c("mean", "dispersion") ]
ncells_per_perturbation_morris_large <- rowSums(counts(altExps(morris_large_data)[["cre_pert"]])) %>% 
  data.frame() %>% rownames_to_column("SNP") %>% mutate(SNP = gsub("-[0-9]$", "", SNP)) %>%
  group_by(SNP) %>% summarize(ncells = sum(.)) %>%
  pull(ncells, name = SNP)

nguides_per_perturbation_morris_large <- rowSums(counts(altExps(morris_large_data)[["cre_pert"]])) %>% 
  data.frame() %>% rownames_to_column("SNP") %>% mutate(SNP = gsub("-[0-9]$", "", SNP)) %>%
  group_by(SNP) %>% summarize(nguides = n()) %>%
  pull(nguides, name = SNP)

guide_distr_test <- rowSums(counts(altExps(morris_large_data)[["cre_pert"]])) %>% 
  data.frame() %>% rownames_to_column("SNP") %>% mutate(SNP_2 = gsub("-[0-9]$", "", SNP))

# process and plot results of power analysis 
output_list <- lapply(list.files("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/simulation_results/gasperini_results_new/", full.name = T), function(x){
  data = readRDS(x)
  data$effect_size <- as.numeric(gsub("sim_sceptre_res_all_|_[0-9]*.rds", "", basename(x)))
  data
}) %>%
  do.call("rbind", .)

output_list %>%
  group_by(gene, cre_pert) %>%
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

# barplots
ntested <- paste0(output_processed$gene, "_", output_processed$cre_pert) %>% unique() %>% length()
output_processed %>%
  dplyr::filter(fraction_sig >= .8) %>%
  mutate(effect_size = factor(effect_size, levels = sort(unique(effect_size), decreasing = T))) %>%
  ggplot(aes(x = effect_size, fill = effect_size)) + 
    ylim(c(0, ntested)) + 
    geom_hline(yintercept = ntested) + 
    theme_classic(base_size = 15) + 
    xlab("Effect Size") + ylab("Number of powered CRE-gene pairs") + labs("fill" = "Effect size") + 
    ggtitle("CRISPRi effect detection power for tested \n CRE−gene links (SCEPTRE, Gasperini data)")
ggsave("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/plots/241216_gasperini_5it_benjaminihochberg_barplot.pdf")

output_processed %>%
  dplyr::mutate(fraction_sig = fraction_sig >= .8) %>%
  group_by(fraction_sig, effect_size) %>%
  summarize(n = n()) %>%
  pivot_wider(values_from = n, names_from = fraction_sig) %>%
  mutate(frac = `TRUE` / (`FALSE` + `TRUE`))

output_processed %>% write_csv("../results/gasperini_power_per_pair.csv")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### morris small analysis
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

output_list <- lapply(list.files("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/simulation_results/morris_small_results_new/", full.name = T), function(x){
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

output_list <- lapply(list.files("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/simulation_results/morris_large_results_new/", full.name = T), function(x){
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
  ggplot(aes(x = fraction_sig)) + geom_histogram()

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
  add_row(effect_size = unique(.$effect_size), fraction_sig = 1, fraction_total = 0) %>%
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
  group_by(gene, dataset) %>%
  summarize(expression_level = unique(expression_level)) %>%
  ggplot(aes(x = expression_level, col = dataset)) + geom_density() + scale_x_log10() + theme_paper()

output_processed %>%
  group_by(gene, dataset) %>%
  summarize(expression_level = unique(expression_level), expression_dispersion = unique(expression_dispersion)) %>%
  ggplot(aes(x = expression_level, y = expression_dispersion, col = dataset)) + geom_point(alpha = .05) + scale_x_log10() + theme_paper() + 
    scale_x_log10() + scale_y_log10() + geom_smooth()

output_processed %>%
  group_by(cre_pert, dataset) %>%
  summarize(number_of_cells = unique(number_of_cells)) %>%
  ggplot(aes(x = number_of_cells, col = dataset)) + geom_density() + theme_paper()

output_processed %>%
  group_by(gene, dataset) %>%
  summarize(expression_level = unique(expression_level)) %>%
  ggplot(aes(x = expression_level, fill = dataset)) + geom_histogram() + scale_x_log10()

### 

output_processed %>%
  dplyr::filter(effect_size == 0.5) %>%
  mutate(number_of_cells = cut(number_of_cells, breaks = 0:20 * 100)) %>%
  dplyr::filter(!is.na(number_of_cells)) %>%
  ggplot(aes(x = number_of_cells, y = fraction_sig)) + stat_summary() + facet_wrap(~dataset) + 
    theme_paper(textsize = 15) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

output_processed %>%
  dplyr::filter(!is.na(number_of_cells)) %>%
  dplyr::filter(effect_size == 0.5) %>%
  ggplot(aes(x = number_of_cells, y = expression_level)) + 
    stat_density_2d(geom = "polygon", aes(alpha = ..level..)) + 
    # geom_point(size = .1, alpha = .1) + 
    facet_wrap(~dataset) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    theme_paper(textsize = 10) + scale_color_manual(values = c("white", "black")) + 
    scale_y_log10()

output_processed %>%
  dplyr::filter(!is.na(number_of_cells)) %>%
  dplyr::filter(effect_size == .5) %>%
  mutate(number_of_cells = cut(number_of_cells, breaks = 1:20 * 100)) %>%
  mutate(expression_level = cut(expression_level, c(0, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000))) %>%
  group_by(number_of_cells, expression_level, dataset) %>%
  summarize(fraction_sig = mean(fraction_sig)) %>%
  dplyr::filter(!is.na(number_of_cells)) %>%
  ggplot(aes(x = number_of_cells, y = expression_level)) + geom_tile(aes(fill = fraction_sig, col = fraction_sig >= .8), size = 1) + facet_wrap(~dataset) + 
    scale_fill_gradient(low = "white", high = "darkred") + 
    theme_paper(textsize = 10) + scale_color_manual(values = c("white", "black")) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle("Effect size = 0.5") + 
    xlab("Number of cells") + ylab("Expression level")
ggsave("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/plots/summarized_power_05.pdf")

output_processed %>%
  dplyr::filter(!is.na(number_of_cells)) %>%
  dplyr::filter(effect_size == .9) %>%
  mutate(number_of_cells = cut(number_of_cells, breaks = 1:20 * 100)) %>%
  mutate(expression_level = cut(expression_level, c(0, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000))) %>%
  group_by(number_of_cells, expression_level, dataset) %>%
  summarize(fraction_sig = mean(fraction_sig)) %>%
  dplyr::filter(!is.na(number_of_cells)) %>%
  ggplot(aes(x = number_of_cells, y = expression_level)) + geom_tile(aes(fill = fraction_sig, col = fraction_sig >= .8), size = 1) + facet_wrap(~dataset) + 
  scale_fill_gradient(low = "white", high = "darkred") + 
  theme_paper(textsize = 10) + scale_color_manual(values = c("white", "black")) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle("Effect size = 0.9") + 
  xlab("Number of cells") + ylab("Expression level")
ggsave("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/plots/summarized_power_09.pdf")

output_processed %>%
  dplyr::filter(!is.na(number_of_cells)) %>%
  dplyr::filter(effect_size == .5) %>%
  mutate(number_of_cells = cut(number_of_cells, breaks = 1:20 * 100)) %>%
  mutate(expression_level = cut(expression_level, c(0, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000))) %>%
  group_by(number_of_cells, expression_level, dataset) %>%
  summarize(fraction_sig = mean(fraction_sig)) %>%
  dplyr::filter(!is.na(number_of_cells)) %>%
  group_by(number_of_cells, dataset, expression_level) %>%
  summarize(fraction_sig = mean(fraction_sig)) %>%
  ggplot(aes(x = number_of_cells, y = fraction_sig, col = expression_level, group = expression_level)) + 
    geom_point() + geom_line() + 
    facet_wrap(~dataset)  + 
    theme_paper(textsize = 10) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle("Effect size = 0.5") + 
    xlab("") + ylab("Average Power")
ggsave("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/plots/summarized_power_lines_05.pdf")

output_processed %>%
  dplyr::filter(!is.na(number_of_cells)) %>%
  dplyr::filter(effect_size == .9) %>%
  mutate(number_of_cells = cut(number_of_cells, breaks = 1:20 * 100)) %>%
  mutate(expression_level = cut(expression_level, c(0, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000))) %>%
  group_by(number_of_cells, expression_level, dataset) %>%
  summarize(fraction_sig = mean(fraction_sig)) %>%
  dplyr::filter(!is.na(number_of_cells)) %>%
  group_by(number_of_cells, dataset, expression_level) %>%
  summarize(fraction_sig = mean(fraction_sig)) %>%
  ggplot(aes(x = number_of_cells, y = fraction_sig, col = expression_level, group = expression_level)) + 
  geom_point() + geom_line() + 
  facet_wrap(~dataset)  + 
  theme_paper(textsize = 10) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle("Effect size = 0.9") + 
  xlab("") + ylab("Average Power")
ggsave("~/Desktop/PostDoc_TL_Lab/Projects/Sam/results/plots/summarized_power_lines_09.pdf")
