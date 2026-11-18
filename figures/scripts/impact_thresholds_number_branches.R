### This script explores the impact of different threshold definitions on the number of saltational branches

library(tidyverse)
library(ggrepel)
library(glue)
library(ggtext)
library(patchwork)
source('config_plot.R')

# Load child parent relationships
df_child_parent <- read_tsv('../../results/bte/df_parent_child_relationship.tsv')

# Load number of AA muts per branch
df_aa_mut_counts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  rename(strain = node_id) %>% group_by(strain) %>% 
  summarise(n_muts = n(), n_syn = sum(aa_initial == aa_new), n_nonsyn = sum(aa_initial != aa_new))

df_child_parent <- df_child_parent %>% 
  left_join(df_aa_mut_counts, by = 'strain') %>% 
  mutate(n_muts = replace_na(n_muts, 0), n_syn = replace_na(n_syn, 0), n_nonsyn = replace_na(n_nonsyn, 0))


input_data_genome <- read_tsv('../../results/input_model/data_viridian_Genome.tsv')
input_data_genome_group <- read_tsv('../../results/input_model/data_viridian_Genome_grouped.tsv')
df_param_estim_genome <- read_tsv('../../results/output_model/df_BF_Genome_CV_1.8.tsv') %>% filter(param == 'omega')

input_data_genome <- input_data_genome %>% 
  left_join(input_data_genome_group) %>% 
  dplyr::select(strain, i_row) %>% 
  left_join(df_param_estim_genome, by = c('i_row'))

vec_BF <- c('5', '10', '20', '40', '80')

# Load number of saltational branches per combination
df_n_saltational_branches_per_threshold <- bind_rows(
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    n_saltational_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_', curr_BF, '.tsv'), 
                                         header = F) %>% unlist() %>% as.character() %>% length()
    tibble(BF_threshold = curr_BF,
           nonsyn_count_percentile = 0.995,
           n_saltational_branches)
  })),
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    n_saltational_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_perc_99.9_CV_1.8_BF_', curr_BF, '.tsv'), 
                                         header = F) %>% unlist() %>% as.character() %>% length()
    tibble(BF_threshold = curr_BF,
           nonsyn_count_percentile = 0.999,
           n_saltational_branches)
  })),
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    n_saltational_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_perc_99_CV_1.8_BF_', curr_BF, '.tsv'), 
                                         header = F) %>% unlist() %>% as.character() %>% length()
    tibble(BF_threshold = curr_BF,
           nonsyn_count_percentile = 0.99,
           n_saltational_branches)
  }))
)

# Load dN/dS at the genome level per combination
df_dNdS_by_branch_type_by_threshold <- bind_rows(
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    name_salt_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_', curr_BF, '.tsv'), 
                                     header = F) %>% unlist() %>% as.character()
    name_adapt_branch_genome_salt_branches <- read.table(paste0('../../results/adaptive_branches/name_adaptive_branches_Genome_CV_1.8_BF_', curr_BF, '.tsv'), 
                                                         header = F) %>% unlist() %>% as.character()
    
    input_data_genome %>% 
      mutate(is_saltational = (strain %in% name_salt_branches) & (strain %in% name_adapt_branch_genome_salt_branches)) %>% 
      group_by(is_saltational) %>% 
      rename(y = posterior_mean) %>% 
      summarise(y0_posterior_mean = quantile(y, 0.025),
                y25_posterior_mean = quantile(y, 0.25),
                y50_posterior_mean= median(y),
                y75_posterior_mean = quantile(y, 0.75),
                y100_posterior_mean = quantile(y, 0.975)) %>% 
      mutate(BF_threshold = curr_BF,
             nonsyn_count_percentile = 0.995)
  })),
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    name_salt_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_perc_99.9_CV_1.8_BF_', curr_BF, '.tsv'), 
                                     header = F) %>% unlist() %>% as.character()
    name_adapt_branch_genome_salt_branches <- read.table(paste0('../../results/adaptive_branches/name_adaptive_branches_Genome_CV_1.8_BF_', curr_BF, '.tsv'), 
                                                         header = F) %>% unlist() %>% as.character()
    
    input_data_genome %>% 
      mutate(is_saltational = (strain %in% name_salt_branches) & (strain %in% name_adapt_branch_genome_salt_branches)) %>% 
      group_by(is_saltational) %>% 
      rename(y = posterior_mean) %>% 
      summarise(y0_posterior_mean = quantile(y, 0.025),
                y25_posterior_mean = quantile(y, 0.25),
                y50_posterior_mean= median(y),
                y75_posterior_mean = quantile(y, 0.75),
                y100_posterior_mean = quantile(y, 0.975)) %>% 
      mutate(BF_threshold = curr_BF,
             nonsyn_count_percentile = 0.999)
  })),
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    name_salt_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_perc_99_CV_1.8_BF_', curr_BF, '.tsv'), 
                                     header = F) %>% unlist() %>% as.character()
    name_adapt_branch_genome_salt_branches <- read.table(paste0('../../results/adaptive_branches/name_adaptive_branches_Genome_CV_1.8_BF_', curr_BF, '.tsv'), 
                                                         header = F) %>% unlist() %>% as.character()
    input_data_genome %>% 
      mutate(is_saltational = (strain %in% name_salt_branches) & (strain %in% name_adapt_branch_genome_salt_branches)) %>% 
      group_by(is_saltational) %>% 
      rename(y = posterior_mean) %>% 
      summarise(y0_posterior_mean = quantile(y, 0.025),
                y25_posterior_mean = quantile(y, 0.25),
                y50_posterior_mean= median(y),
                y75_posterior_mean = quantile(y, 0.75),
                y100_posterior_mean = quantile(y, 0.975)) %>% 
      mutate(BF_threshold = curr_BF,
             nonsyn_count_percentile = 0.99)
  }))
)

# Load mutation count per combination
df_mean_mut_count_by_branch_type_by_threshold <- bind_rows(
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    name_salt_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_', curr_BF, '.tsv'), 
                                     header = F) %>% unlist() %>% as.character()
    
    df_child_parent %>% mutate(is_salt = strain %in% name_salt_branches) %>% 
      summarise(mean_nonsyn_salt = mean(n_nonsyn[is_salt]),
                mean_nonsyn_non_salt = mean(n_nonsyn[!is_salt]),
                median_nonsyn_salt = median(n_nonsyn[is_salt]),
                median_nonsyn_non_salt = median(n_nonsyn[!is_salt]),
                lower_IQR_salt = quantile(n_nonsyn[is_salt], 0.25),
                lower_IQR_non_salt = quantile(n_nonsyn[!is_salt], 0.25),
                upper_IQR_salt = quantile(n_nonsyn[is_salt], 0.75),
                upper_IQR_non_salt = quantile(n_nonsyn[!is_salt], 0.75)) %>% 
      mutate(BF_threshold = curr_BF,
             nonsyn_count_percentile = 0.995)
  })),
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    name_salt_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_perc_99.9_CV_1.8_BF_', curr_BF, '.tsv'), 
                                     header = F) %>% unlist() %>% as.character()
    df_child_parent %>% mutate(is_salt = strain %in% name_salt_branches) %>% 
      summarise(mean_nonsyn_salt = mean(n_nonsyn[is_salt]),
                mean_nonsyn_non_salt = mean(n_nonsyn[!is_salt]),
                median_nonsyn_salt = median(n_nonsyn[is_salt]),
                median_nonsyn_non_salt = median(n_nonsyn[!is_salt]),
                lower_IQR_salt = quantile(n_nonsyn[is_salt], 0.25),
                lower_IQR_non_salt = quantile(n_nonsyn[!is_salt], 0.25),
                upper_IQR_salt = quantile(n_nonsyn[is_salt], 0.75),
                upper_IQR_non_salt = quantile(n_nonsyn[!is_salt], 0.75)) %>% 
      mutate(BF_threshold = curr_BF,
             nonsyn_count_percentile = 0.999)
  })),
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    name_salt_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_perc_99_CV_1.8_BF_', curr_BF, '.tsv'), 
                                     header = F) %>% unlist() %>% as.character()
    df_child_parent %>% mutate(is_salt = strain %in% name_salt_branches) %>% 
      summarise(mean_nonsyn_salt = mean(n_nonsyn[is_salt]),
                mean_nonsyn_non_salt = mean(n_nonsyn[!is_salt]),
                median_nonsyn_salt = median(n_nonsyn[is_salt]),
                median_nonsyn_non_salt = median(n_nonsyn[!is_salt]),
                lower_IQR_salt = quantile(n_nonsyn[is_salt], 0.25),
                lower_IQR_non_salt = quantile(n_nonsyn[!is_salt], 0.25),
                upper_IQR_salt = quantile(n_nonsyn[is_salt], 0.75),
                upper_IQR_non_salt = quantile(n_nonsyn[!is_salt], 0.75)
      ) %>% 
      mutate(BF_threshold = curr_BF,
             nonsyn_count_percentile = 0.99)
  }))
)
    
## Visualize impact of criterion definition on branch characteristics
highlight_labels <- function(x, pat, color = "black", family = "") {
  ifelse(grepl(pat, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}



plt_dNdS_func_threshold <- df_dNdS_by_branch_type_by_threshold %>% 
  mutate(BF_threshold_label = case_when(
    (nonsyn_count_percentile == 0.995 & BF_threshold == "20") ~ "**20**",
    TRUE ~ as.character(BF_threshold)),
    BF_threshold_label = factor(BF_threshold_label, c('5', '10', '20', '**20**', '40', '80'))) %>%
  ggplot(aes(x = BF_threshold_label, group = interaction(is_saltational, BF_threshold),
             fill = is_saltational, colour = is_saltational)) +
  geom_boxplot(
    aes(ymin = y0_posterior_mean, 
        lower = y25_posterior_mean, 
        middle = y50_posterior_mean, 
        upper = y75_posterior_mean, 
        ymax = y100_posterior_mean),
    stat = "identity", width = 0.5, alpha = 0.3
  ) +
  scale_x_discrete(name = 'Bayes Factor threshold') +
  scale_y_continuous(expression(atop("Estimated genome level", d[N]/d[S]~"(posterior mean)"))) +
  scale_color_manual(name = 'Branch type', values = c(col_adaptive, col_non_adaptive),
                     breaks = c(T, F),
                     labels = c('Saltational\n(genome level)', 'Other')) +
  scale_fill_manual(name = 'Branch type', values = c(col_adaptive, col_non_adaptive), 
                    breaks = c(T, F), 
                    labels = c('Saltational\n(genome level)', 'Other')) +
  facet_wrap(.~ nonsyn_count_percentile, 
             scales = 'free_x',
             labeller = as_labeller(function(x) highlight_labels(paste0('Percentile nonsyn count: ', scales::percent(as.numeric(x))),
                                                                 'Percentile nonsyn count: 99.5%'))) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        strip.text.x = element_markdown(),
        axis.text.x = element_markdown())


plt_nonsyn_mut_count <- df_mean_mut_count_by_branch_type_by_threshold %>% 
  mutate(BF_threshold_label = case_when(
    (nonsyn_count_percentile == 0.995 & BF_threshold == "20") ~ "**20**",
    TRUE ~ as.character(BF_threshold)),
    BF_threshold_label = factor(BF_threshold_label, c('5', '10', '20', '**20**', '40', '80'))) %>%
  ggplot(aes(x = BF_threshold_label)) +
  geom_point(aes(y = median_nonsyn_salt, colour = 'Saltational')) +
  geom_linerange(aes(ymin = lower_IQR_salt, ymax = upper_IQR_salt, colour = 'Saltational')) +
  geom_point(aes(y = median_nonsyn_non_salt, colour = 'Non-saltational')) +
  geom_linerange(aes(ymin = lower_IQR_non_salt, ymax = upper_IQR_non_salt, colour = 'Non-saltational')) +
  scale_x_discrete(name = 'Bayes Factor threshold') +
  scale_y_continuous(name = 'Nonsynonymous\nmutation count') +
  scale_color_manual(name = 'Branch type', values = c(col_adaptive, col_non_adaptive), breaks = c('Saltational', 'Non-saltational')) +
  facet_wrap(.~ nonsyn_count_percentile, scales = 'free_x',
             labeller = as_labeller(function(x) highlight_labels(paste0('Percentile nonsyn count: ', scales::percent(as.numeric(x))),
                                                                 'Percentile nonsyn count: 99.5%'))) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        strip.text.x = element_markdown(),
        axis.text.x = element_markdown())

plt_nb_salt_branch_func_threshold <- df_n_saltational_branches_per_threshold %>% 
  mutate(BF_threshold_label = case_when(
    (nonsyn_count_percentile == 0.995 & BF_threshold == "20") ~ "**20**",
    TRUE ~ as.character(BF_threshold)),
    BF_threshold_label = factor(BF_threshold_label, c('5', '10', '20', '**20**', '40', '80'))) %>%
  ggplot(aes(x = BF_threshold_label, y = n_saltational_branches)) +
  geom_point(stat = 'identity') +
  geom_text_repel(aes(label = n_saltational_branches)) +
  facet_wrap(.~ nonsyn_count_percentile, 
             scales = 'free_x',
             labeller = as_labeller(function(x) highlight_labels(paste0('Percentile nonsyn count: ', scales::percent(as.numeric(x))),
                                                                 'Percentile nonsyn count: 99.5%'))) +
  scale_x_discrete(name = 'Bayes Factor threshold') +
  scale_y_continuous(name = 'Number of\nsaltational branches', trans = 'log',
                     breaks = c(10, 100, 1000, 10000)) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        strip.text.x = element_markdown(),
        axis.text.x = element_markdown())

# Make panel and save results  
panel_impact_threshold <- (plt_nb_salt_branch_func_threshold / plt_nonsyn_mut_count / plt_dNdS_func_threshold) &
  theme(plot.tag = element_text(face = "bold", size = 16)) & 
  plot_annotation(tag_levels = 'A') &
  theme(plot.margin = unit(c(.01,.01,.01,.01), "cm"))

plot(panel_impact_threshold)

png('../figure_model_illustration/impact_threshold_char_salt_branches.png',
    height = 7., width = 8.5, res = 350, units = 'in')
plot(panel_impact_threshold)
dev.off()
  


