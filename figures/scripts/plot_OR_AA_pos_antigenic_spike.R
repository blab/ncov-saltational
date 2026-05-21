### This script reproduces figures depicting the 
### AA positions in Spike that are more likely to 
### have nonsynonymous mutation in saltational 
### than non saltattional branches (Figure 3E)

library(tidyverse)
library(ggrepel)
source('config_plot.R')

#############################
### LOADING MODEL RESULTS ###
#############################

# Load name of saltational branches
name_saltational_branches <- 
  read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
             header = F) %>% unlist() %>% as.character()

# Load AA mutations information and get proportion of adaptive branches with each mutation
df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  filter(is_coded) %>% # Coded AA mut
  filter(aa_new != aa_initial) %>% # Non-synonymous
  rename(strain = node_id) %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(has_mut = T) %>% 
  select(strain, gene, aa_muts_position, has_mut) %>% 
  group_by(gene, aa_muts_position) %>% 
  summarise(n_adapt_branches_with_mut = n()) %>% 
  mutate(prop_adapt_branches_with_mut = n_adapt_branches_with_mut / length(name_saltational_branches))

# Load sites associated with antibody escape
df_antigenic_sites <- read_csv('../../data/antigenic_sites_escape_antibodies.csv')

# Load OR of NS mutation in saltational branches by AA position
df_OR_AA_pos_spike <- read_tsv('../../results/mutation_enrichment/df_AA_pos_enrichment_saltational_CV_1.8_BF_20.tsv') %>% 
  filter(threshold_for_inclusion == 2) %>% 
  filter(gene == 'S') %>% 
  mutate(is_antigenic = position %in% df_antigenic_sites$position)

###################
### MAKE FIGURE ###
###################
get_figure_3e <- function(){
  plt_OR_func_n_mut_salt_branches <- df_OR_AA_pos_spike %>% 
    left_join(df_aa_muts, by = c('gene', 'position' = 'aa_muts_position')) %>% 
    mutate(region = case_when(position >= 13 & position <= 302 ~ 'NTD',
                              position >= 319 & position <= 541 ~ 'RBD',
                              position >= 13 & position <= 685 ~ 'Other S1',
                              position >= 686 & position <= 1273 ~ 'S2',
                              T ~ 'Other S')) %>% 
    mutate(region = factor(region, levels = names(palette_spike))) %>% 
    ggplot(aes(x = n_adapt_branches_with_mut, 
               y = OR, 
               alpha = p_value_BH < 0.05)) +
    geom_point(shape = 21,
               aes(stroke = is_antigenic,
                   fill = region),
               size = 2) +
    geom_text_repel(data = df_OR_AA_pos_spike %>%
                      left_join(df_aa_muts, by = c('gene', 'position' = 'aa_muts_position')) %>%
                      filter(p_value_BH < 0.05, is_antigenic),
                    aes(label = position)) +
    scale_x_continuous(name = 'Number of saltational branches with\nnonsynonymous mutation at Spike position') +
    scale_y_continuous(name = 'OR', trans = 'log',
                       breaks = c(0.1, 0.2, 0.5, 1.0, 2., 5., 10., 20, 50, 100)) +
    scale_alpha_manual(values = c(0.2, 1.0),
                       breaks = c(F, T), guide = 'none') +
    scale_discrete_manual(
      aesthetics = "stroke",
      values = c('FALSE' = 0, 'TRUE' = .8),
      name = '',
      labels = c('', 'Antibody escape site')
    ) +
    scale_fill_manual(values = palette_spike, name = 'Spike region / domain') +
    coord_cartesian(ylim = c(1., NA)) +
    theme_classic() +
    guides(fill = guide_legend(override.aes = list(label = '', stroke = 0))) +
    theme(legend.background = element_blank(),
          axis.text = element_text(size = 11))
  
  return(plt_OR_func_n_mut_salt_branches)
}

plt_OR_func_n_mut_salt_branches <- get_figure_3e()

pdf('../figure_mutation_enrichment/scatter_plot_OR_func_pos_count_saltational_within_S.pdf',
    height = 3., width = 5.5)
plot(plt_OR_func_n_mut_salt_branches)
dev.off()
png('../figure_mutation_enrichment/scatter_plot_OR_func_pos_count_saltational_within_S.png',
    height = 3., width = 5.5, res = 350, units = 'in')
plot(plt_OR_func_n_mut_salt_branches)
dev.off()


# df_OR_AA_pos_spike %>% 
#   filter(p_value_BH < 0.05, OR > 1.) %>% 
#   mutate(region = case_when(position >= 13 & position <= 302 ~ 'NTD',
#                             position >= 319 & position <= 541 ~ 'RBD',
#                             position >= 13 & position <= 685 ~ 'Other S1',
#                             position >= 686 & position <= 1273 ~ 'S2',
#                             T ~ 'Other S')) %>% 
#   filter(region == 'Other S1') %>% 
#   arrange(position) %>% select(position) %>% unlist() %>% as.character() %>% paste(collapse = '+')
