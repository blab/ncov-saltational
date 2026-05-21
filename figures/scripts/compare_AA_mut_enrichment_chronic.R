library(tidyverse)
library(ggrepel)
library(ggpubr)
source('config_plot.R')

# Load OR by AA mutation
df_OR_AA_mut <- read_tsv('../../results/mutation_enrichment/df_AA_mut_enrichment_saltational_opportunity_CV_1.8_BF_20.tsv') %>% 
  filter(threshold_for_inclusion == 1) %>% 
  filter(p_value_BH < 0.05)

df_OR_AA_mut %>% filter(gene == 'S') %>%
  dplyr::select(mutation, OR, OR_lower, OR_upper, n_mut_of_interest_adapt_branches) %>% 
  group_by(mutation) %>% mutate(position = as.integer(substr(mutation, start = 1, stop = nchar(mutation) - 1))) %>% 
  arrange(position) %>% 
  View()
