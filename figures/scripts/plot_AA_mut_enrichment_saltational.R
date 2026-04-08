library(tidyverse)
library(ggrepel)
library(ggpubr)
source('config_plot.R')

# Load OR by AA mutation
df_OR_AA_mut <- read_tsv('../../results/mutation_enrichment/df_AA_mut_enrichment_saltational_CV_1.7_BF_20.tsv') %>% 
  filter(threshold_for_inclusion == 1)

# Load name of adaptive paths
vec_name_node_adapt <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.7_BF_20.tsv') %>% unlist() %>% as.character()

# Load AA mutations information and get proportion of adaptive branches with each mutation
df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  filter(is_coded) %>% # Coded AA mut
  filter(aa_new != aa_initial) %>% # Non-synonymous
  mutate(full_mut = paste0(aa_muts_position, aa_new)) %>% 
  rename(strain = node_id) %>% 
  filter(strain %in% vec_name_node_adapt) %>% 
  mutate(has_mut = T) %>% 
  select(strain, gene, full_mut, has_mut) %>% 
  group_by(gene, full_mut) %>% 
  summarise(n_adapt_branches_with_mut = n()) %>% 
  mutate(prop_adapt_branches_with_mut = n_adapt_branches_with_mut / length(vec_name_node_adapt))

# Scatter plots of OR as a function of the proportion of branches with that mutation
get_plot_OR_func_prop_adapt_branches_labels_zoom_S_regions <- function(df_OR_AA_mut, df_aa_muts,
                                                                       prop_threshold = 0., OR_threshold = 1.,
                                                                       xmax = NA, ymin = NA, ymax = NA){
  OR_value_inf <- max(df_OR_AA_mut$OR[is.finite(df_OR_AA_mut$OR)])*3
  
  label_125 <- function(x) {
    m <- signif(x / 10^floor(log10(x)), 12)   # mantissa
    labs <- scales::label_percent(accuracy = 0.1)(x)
    ifelse(m %in% c(1, 2, 5), labs, "")
  }
  
  tmp_df <- df_OR_AA_mut %>% 
    filter(gene == 'S') %>% 
    group_by(gene, mutation) %>%
    mutate(position_within_gene = as.integer(substr(mutation, start = 1, stop = nchar(mutation) - 1))) %>% 
    mutate(region = case_when(position_within_gene >= 13 & position_within_gene <= 302 ~ 'NTD',
                              position_within_gene >= 319 & position_within_gene <= 541 ~ 'RBD',
                              position_within_gene >= 13 & position_within_gene <= 685 ~ 'Other S1',
                              position_within_gene >= 686 & position_within_gene <= 1273 ~ 'S2',
                              T ~ 'Other S')) %>% 
    left_join(df_aa_muts, by = c('mutation' = 'full_mut', 'gene')) 
  

  tmp_df <- tmp_df %>% 
    mutate(OR = ifelse(is.infinite(OR), OR_value_inf, OR))
  
  plot_OR_func_prop_branches_with_mut <- tmp_df %>% 
    ggplot(aes(x = prop_adapt_branches_with_mut, 
               y = OR, colour = region)) +
    geom_point(aes(alpha = p_value_BH < 0.05)) +
    geom_text_repel(data = tmp_df %>% 
                      filter(prop_adapt_branches_with_mut > prop_threshold, 
                             OR > OR_threshold, 
                             p_value_BH < 0.05),
                    aes(label = mutation), 
                    size = 3.5,
                    max.time = 20, min.segment.length = 0., max.overlaps = Inf) +
    scale_x_continuous(labels = label_125,
                       trans = 'log',
                       expand = expansion(mult = c(0.065, 0.06)),
                       breaks = c(seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 9e-1, 1e-1)),
                       name = 'Proportion of adaptive branches with mutation')  +
    scale_y_continuous(name = 'OR',
                       trans = 'log', 
                       expand = expansion(mult = c(0.02, 0.02)),
                       breaks = c(0.1, 0.2, 0.5, 1.0, 2.0, 5, 10, 20, 50, 100, 200, 500, 1000, OR_value_inf),
                       labels = c(0.1, 0.2, 0.5, 1.0, 2.0, 5, 10, 20, 50, 100, 200, 500, 1000, 'Inf')
    ) +
    scale_colour_manual(breaks = names(palette_spike), values = palette_spike, name = 'Region') +
    scale_alpha_manual(breaks = c(T, F), values = c(0.9, .05), guide = 'none') +
    theme_classic() +
    coord_cartesian(xlim = c(NA, xmax), ylim = c(ymin, ymax)) +
    theme(axis.title = element_text(size = 11), axis.text = element_text(size = 11),
          strip.background = element_blank(), strip.text = element_blank()) +
    guides(colour = guide_legend(override.aes = list(label = '')))
  
  return(plot_OR_func_prop_branches_with_mut)
}

get_plot_OR_func_prop_adapt_branches_labels_zoom_genes <- function(df_OR_AA_mut, df_aa_muts, 
                                                                   vec_genes, 
                                                                   prop_threshold = 0., OR_threshold = 1., 
                                                                   xmax = NA, ymin = NA, ymax = NA){
  
  OR_value_inf <- Inf
  
  label_125 <- function(x) {
    m <- signif(x / 10^floor(log10(x)), 12) 
    labs <- scales::label_percent(accuracy = 0.1)(x)
    ifelse(m %in% c(1, 2, 5), labs, "")
  }
  label_125_no_percent <- function(x) {
    m <- signif(x / 10^floor(log10(x)), 12) 
    labs <- x
    ifelse(m %in% c(1, 2, 5), labs, "")
  }
  
  plot_OR_func_prop_branches_with_mut <- df_OR_AA_mut %>% 
    filter(gene %in% vec_genes) %>% 
    left_join(df_aa_muts, by = c('mutation' = 'full_mut', 'gene')) %>% 
    mutate(OR = ifelse(is.infinite(OR), OR_value_inf, OR)) %>% 
    ggplot(aes(x = prop_adapt_branches_with_mut, 
               y = OR, colour = gene)) +
    geom_point(aes(alpha = p_value_BH < 0.05)) +
    geom_text_repel(data = df_OR_AA_mut %>% 
                      left_join(df_aa_muts, by = c('mutation' = 'full_mut', 'gene')) %>% 
                      mutate(OR = ifelse(is.infinite(OR), OR_value_inf, OR)) %>% 
                      filter(gene %in% vec_genes) %>% 
                      filter(prop_adapt_branches_with_mut > prop_threshold, 
                             OR > OR_threshold, 
                             p_value_BH < 0.05),
                    aes(label = mutation),
                    size = 3.5,
                    max.time = 20, min.segment.length = 0., max.overlaps = Inf) +
    scale_x_continuous(labels = label_125,
                       trans = 'log',
                       expand = expansion(mult = c(0.065, 0.06)),
                       breaks = c(seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 9e-1, 1e-1)),
                       name = 'Proportion of adaptive branches with mutation')  +
    scale_y_continuous(name = 'OR',
                       labels = label_125_no_percent,
                       trans = 'log', 
                       expand = expansion(mult = c(0.02, 0.02)),
                       breaks = c(seq(1e-1, 9e-1, 1e-1), seq(1, 9, 1), seq(10, 90, 10),
                                  seq(100, 900, 100), seq(1000, 9000, 1000))
    ) +
    scale_colour_manual(breaks = names(palette_genes), values = palette_genes, name = 'Gene') +
    scale_alpha_manual(breaks = c(T, F), values = c(0.9, .05), guide = 'none') +
    theme_classic() +
    coord_cartesian(ylim = c(ymin, ymax)) +
    theme(axis.title = element_text(size = 11), axis.text = element_text(size = 11),
          strip.background = element_blank(), strip.text = element_blank()) +
    guides(colour = guide_legend(override.aes = list(label = '')))
  
  return(plot_OR_func_prop_branches_with_mut)
}


prop_threshold_labels <- 0.015
OR_threshold_labels <- 10
x_max_plot <- NA
y_min_plot <- 1.
y_max_plot <- NA

scatter_plot_outside_S <- get_plot_OR_func_prop_adapt_branches_labels_zoom_genes(df_OR_AA_mut, df_aa_muts,
                                                                                 c('ORF1a', 'ORF1b', 'ORF7a', 'ORF8', 'N', 'E', 'M', 'ORF3a', 'ORF7b', 'ORF9b', 'ORF6'),
                                                                                 prop_threshold_labels, OR_threshold_labels, 
                                                                                 NA, y_min_plot)

scatter_plot_within_S <- get_plot_OR_func_prop_adapt_branches_labels_zoom_S_regions(df_OR_AA_mut, df_aa_muts,
                                                                                    prop_threshold_labels, OR_threshold_labels, 
                                                                                    x_max_plot, y_min_plot)

scatter_plot_all <- get_plot_OR_func_prop_adapt_branches_labels_zoom_genes(df_OR_AA_mut, df_aa_muts,
                                                                            c('S', 'ORF1a', 'ORF1b', 'ORF7a', 'ORF8', 'N', 'E', 'M', 'ORF3a', 'ORF7b', 'ORF9b', 'ORF6'),
                                                                            prop_threshold_labels, OR_threshold_labels, 
                                                                            NA, y_min_plot)

set.seed(15345)
panel_across_genes <- ggarrange(scatter_plot_within_S + 
                                  ggtitle('Within S') + 
                                  guides(colour = guide_legend(nrow = 3, override.aes = list(label = ''))) +
                                  theme(legend.position = 'bottom'),
                                scatter_plot_outside_S + 
                                  ggtitle('Outside S') + 
                                  guides(colour = guide_legend(nrow = 3, override.aes = list(label = ''))) +
                                  theme(legend.position = 'bottom'),
                                ncol = 2, widths = c(1., 1.2))


pdf('../figure_mutation_enrichment/scatter_plot_OR_func_mut_prevalence_saltational.pdf',
    height = 4.5, width = 6.)
plot(scatter_plot_all)
dev.off()
png('../figure_mutation_enrichment/scatter_plot_OR_func_mut_prevalence_saltational.png',
    height = 4.5, width = 6., res = 350, units = 'in')
plot(scatter_plot_all)
dev.off()
