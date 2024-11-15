### 

gasperini <- readRDS("~/Desktop/gasperini_combined_0101124.rds")

# how many cre-gene pairs were tested?
n_tests <- gasperini %>% 
  group_by(perturbation, gene) %>%
  summarize(n_tests = n()) %>% nrow()

p_adj_cutoff <- 0.05

gasperini %>% 
  mutate(pvalue_bonf = pvalue * n_tests) %>%
  group_by(perturbation, gene) %>%
  summarize(n_sig = sum(pvalue < p_adj_cutoff), n_non_sig = sum(pvalue >= p_adj_cutoff)) %>%
  mutate(frac_sig = n_sig / (n_sig + n_non_sig)) %>%
  ggplot(aes(x = frac_sig)) + geom_histogram(bins = 20)
