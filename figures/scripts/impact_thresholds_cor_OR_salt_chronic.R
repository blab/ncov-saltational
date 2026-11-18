library(tidyverse)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(binom)
library(glue)
library(ggtext)
library(patchwork)
source('utils_plot_chronic.R')
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

vec_BF <- c('5', '10', '20', '40', '80')

## Impact of thresholds used to define saltational branches
## on the cor between the log OR of a mut occurring in 
## saltational vs non-saltational branches and in persistent
## infections vs background sequences

df_pearson_cor_log_OR <- bind_rows(
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    name_salt_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_', curr_BF, '.tsv'), 
                                     header = F) %>% unlist() %>% as.character()
    
    df_OR_chronic_saltational <- compute_OR_muts_chronic_infections_saltational('../../data/frequent_mutations/ryan_epci_hqcs_AA_ratio.csv',
                                                                                '../../results/mutations/nodes_aa_muts_summary_coded.tsv',
                                                                                name_salt_branches)
    
    pearson_cor_OR <- cor(log(df_OR_chronic_saltational$OR_saltational_modif),
                          log(df_OR_chronic_saltational$OR_chronic_modif),
                          method = 'pearson')
    tibble(BF = curr_BF, 
           pearson_cor_log_OR = pearson_cor_OR,
           nonsyn_count_percentile = 0.995)
  })),
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    name_salt_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_perc_99.9_CV_1.8_BF_', curr_BF, '.tsv'), 
                                     header = F) %>% unlist() %>% as.character()
    
    df_OR_chronic_saltational <- compute_OR_muts_chronic_infections_saltational('../../data/frequent_mutations/ryan_epci_hqcs_AA_ratio.csv',
                                                                                '../../results/mutations/nodes_aa_muts_summary_coded.tsv',
                                                                                name_salt_branches)
    
    pearson_cor_OR <- cor(log(df_OR_chronic_saltational$OR_saltational_modif),
                          log(df_OR_chronic_saltational$OR_chronic_modif),
                          method = 'pearson')
    tibble(BF = curr_BF, 
           pearson_cor_log_OR = pearson_cor_OR,
           nonsyn_count_percentile = 0.999)
  })),
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    name_salt_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_perc_99_CV_1.8_BF_', curr_BF, '.tsv'), 
                                     header = F) %>% unlist() %>% as.character()
    
    df_OR_chronic_saltational <- compute_OR_muts_chronic_infections_saltational('../../data/frequent_mutations/ryan_epci_hqcs_AA_ratio.csv',
                                                                                '../../results/mutations/nodes_aa_muts_summary_coded.tsv',
                                                                                name_salt_branches)
    
    pearson_cor_OR <- cor(log(df_OR_chronic_saltational$OR_saltational_modif),
                          log(df_OR_chronic_saltational$OR_chronic_modif),
                          method = 'pearson')
    tibble(BF = curr_BF, 
           pearson_cor_log_OR = pearson_cor_OR,
           nonsyn_count_percentile = 0.99)
  }))
)

# Display results
highlight_labels <- function(x, pat, color = "black", family = "") {
  ifelse(grepl(pat, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}

panel_cor_log_OR_impact_threshold <- df_pearson_cor_log_OR %>% 
  mutate(BF_threshold_label = case_when(
    (nonsyn_count_percentile == 0.995 & BF == "20") ~ "**20**",
    TRUE ~ as.character(BF)),
    BF_threshold_label = factor(BF_threshold_label, c('5', '10', '20', '**20**', '40', '80'))) %>%
  ggplot(aes(x = BF_threshold_label, y = pearson_cor_log_OR)) +
  geom_bar(stat = 'identity', fill = 'darkcyan') +
  geom_text(aes(label = round(pearson_cor_log_OR, 2)), vjust = -.8) +
  facet_wrap(.~ nonsyn_count_percentile, 
             scales = 'free_x',
             labeller = as_labeller(function(x) highlight_labels(paste0('Percentile nonsyn count: ', scales::percent(as.numeric(x))),
            'Percentile nonsyn count: 99.5%'))) +
  scale_x_discrete(name = 'Bayes Factor threshold') +
  scale_y_continuous(name = 'Pearson r',
                     breaks = seq(0., 1., 0.2),
                     limits = c(0., 1.)) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        strip.text.x = element_markdown(),
        axis.text.x = element_markdown())
# 
png('../figure_model_illustration/impact_threshold_char_cor_log_OR.png',
    height = 2.8, width = 6.5, res = 350, units = 'in')
plot(panel_cor_log_OR_impact_threshold)
dev.off()

