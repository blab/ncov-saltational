library(tidyverse)
library(ggrepel)
library(ggpubr)
source('config_plot.R')

# Load OR by AA mutation
df_OR_AA_mut <- read_tsv('../../results/mutation_enrichment/df_AA_mut_enrichment_saltational_opportunity_CV_1.8_BF_20.tsv') %>% 
  filter(threshold_for_inclusion == 2)

# Load name of adaptive paths
vec_name_node_adapt <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv') %>% unlist() %>% as.character()

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

get_plot_OR_func_n_mut_within_salt_branches <- function(df_OR_AA_mut, df_aa_muts, 
                                                        vec_genes, 
                                                        mut_threshold = 0., OR_threshold = 1., 
                                                        xmax = NA, ymin = NA, ymax = NA){
  
  OR_value_inf <- Inf
 
  label_125 <- function(x) {
    m <- signif(x / 10^floor(log10(x)), 12) 
    labs <- x
    ifelse(m %in% c(1, 2, 5), labs, "")
  }
  
  plot_OR_func_prop_branches_with_mut <- df_OR_AA_mut %>% 
    filter(gene %in% vec_genes) %>% 
    left_join(df_aa_muts, by = c('mutation' = 'full_mut', 'gene')) %>% 
    mutate(OR = ifelse(is.infinite(OR), OR_value_inf, OR)) %>% 
    ggplot(aes(x = n_mut_of_interest_adapt_branches, 
               y = OR, colour = gene)) +
    geom_point(aes(alpha = p_value_BH < 0.05)) +
    geom_text_repel(data = df_OR_AA_mut %>% 
                      left_join(df_aa_muts, by = c('mutation' = 'full_mut', 'gene')) %>% 
                      mutate(OR = ifelse(is.infinite(OR), OR_value_inf, OR)) %>% 
                      filter(gene %in% vec_genes) %>% 
                      filter(n_mut_of_interest_adapt_branches > mut_threshold, 
                             OR > OR_threshold, 
                             p_value_BH < 0.05),
                    aes(label = mutation),
                    size = 3.5,
                    max.time = 20, min.segment.length = 0., max.overlaps = Inf) +
    scale_x_continuous(labels = label_125,
                       trans = 'log',
                       expand = expansion(mult = c(0.065, 0.06)),
                       breaks = c(seq(1, 9, 1), seq(1e1, 9e1, 1e1), seq(1e2, 9e2, 1e2)),
                       name = 'Number of saltational branches with mutation')  +
    scale_y_continuous(name = 'OR',
                       labels = label_125,
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

get_plot_OR_func_n_mut_within_salt_branches_zoom_S_regions <- function(df_OR_AA_mut, df_aa_muts,
                                                                       mut_threshold = 0., OR_threshold = 1.,
                                                                       xmax = NA, ymin = NA, ymax = NA){
  OR_value_inf <- Inf
  
  label_125 <- function(x) {
    m <- signif(x / 10^floor(log10(x)), 12) 
    labs <- x
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
    ggplot(aes(x = n_mut_of_interest_adapt_branches, 
               y = OR, colour = region)) +
    geom_point(aes(alpha = p_value_BH < 0.05)) +
    geom_text_repel(data = tmp_df %>% 
                      filter(n_mut_of_interest_adapt_branches > mut_threshold, 
                             OR > OR_threshold, 
                             p_value_BH < 0.05),
                    aes(label = mutation), 
                    size = 3.5,
                    max.time = 20, min.segment.length = 0., max.overlaps = Inf) +
    scale_x_continuous(labels = label_125,
                       trans = 'log',
                       expand = expansion(mult = c(0.065, 0.06)),
                       breaks = c(seq(1, 9, 1), seq(1e1, 9e1, 1e1), seq(1e2, 9e2, 1e2)),
                       name = 'Number of saltational branches with mutation')  +
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
prop_threshold_labels <- 0.005
OR_threshold_labels <- 5
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

scatter_plot_all_func_n_mut <- get_plot_OR_func_n_mut_within_salt_branches(df_OR_AA_mut, df_aa_muts,
                                                                           c('S', 'ORF1a', 'ORF1b', 'ORF7a', 'ORF8', 'N', 'E', 'M', 'ORF3a', 'ORF7b', 'ORF9b', 'ORF6'),
                                                                           5, OR_threshold_labels, 
                                                                           NA, y_min_plot)
scatter_plot_outside_S_func_n_mut <- get_plot_OR_func_n_mut_within_salt_branches(df_OR_AA_mut, df_aa_muts,
                                                                                 c('ORF1a', 'ORF1b', 'ORF7a', 'ORF8', 'N', 'E', 'M', 'ORF3a', 'ORF7b', 'ORF9b', 'ORF6'),
                                                                                 5, OR_threshold_labels, 
                                                                                 NA, y_min_plot)

scatter_plot_S_func_n_mut <- get_plot_OR_func_n_mut_within_salt_branches_zoom_S_regions(df_OR_AA_mut, df_aa_muts,
                                                                                        5, OR_threshold_labels, 
                                                                                        NA, y_min_plot)

# Save results
pdf('../figure_mutation_enrichment/scatter_plot_OR_func_mut_count_saltational.pdf',
    height = 4.5, width = 6.)
plot(scatter_plot_all_func_n_mut)
dev.off()
png('../figure_mutation_enrichment/scatter_plot_OR_func_mut_count_saltational.png',
    height = 4.5, width = 6., res = 350, units = 'in')
plot(scatter_plot_all_func_n_mut)
dev.off()

pdf('../figure_mutation_enrichment/scatter_plot_OR_func_mut_count_saltational_outside_S.pdf',
    height = 4.5, width = 6.)
plot(scatter_plot_outside_S_func_n_mut)
dev.off()
set.seed(1)
png('../figure_mutation_enrichment/scatter_plot_OR_func_mut_count_saltational_outside_S.png',
    height = 4.5, width = 6., res = 350, units = 'in')
plot(scatter_plot_outside_S_func_n_mut)
dev.off()

pdf('../figure_mutation_enrichment/scatter_plot_OR_func_mut_count_saltational_within_S.pdf',
    height = 4.5, width = 7.5)
plot(scatter_plot_S_func_n_mut)
dev.off()
set.seed(1)
png('../figure_mutation_enrichment/scatter_plot_OR_func_mut_count_saltational_within_S.png',
    height = 4.5, width = 7.5, res = 350, units = 'in')
plot(scatter_plot_S_func_n_mut)
dev.off()


## Barplot with counts of enriched mutations per positions
get_barplot_nb_enriched_mutations <- function(df_OR_AA_mut, threshold_inclusion){
  barplot_nb_enriched_positions_by_gene <- tibble(gene = names(palette_genes)[names(palette_genes) != 'Genome']) %>% 
    left_join(
      df_OR_AA_mut %>% 
        filter(threshold_for_inclusion == threshold_inclusion) %>% filter(p_value_BH < 0.05) %>% 
        group_by(gene) %>% summarise(n_pos = n()),
      by = 'gene'
    ) %>% 
    mutate(gene = factor(gene, levels = names(palette_genes)[names(palette_genes) != 'Genome']),
           n_pos = replace_na(n_pos, 0)) %>% 
    ggplot(aes(x = gene, y = n_pos, fill = gene)) +
    geom_bar(stat = 'identity') +
    geom_text(aes(label = n_pos, colour = gene), vjust = -0.5) +
    scale_x_discrete(name = 'Gene') +
    scale_y_continuous(name = 'Number of enriched mutations', 
                       expand = expansion(mult = c(0., 0.1))) +
    scale_fill_manual(breaks = names(palette_genes), values = palette_genes, name = 'Gene', 
                      guide = 'none') +
    scale_colour_manual(breaks = names(palette_genes), values = palette_genes, name = 'Gene', 
                        guide = 'none') +
    theme_classic() +
    theme(axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1., vjust = 1.))
  
  return(barplot_nb_enriched_positions_by_gene)
}
get_barplot_nb_enriched_mutations_within_S <- function(df_OR_AA_mut, threshold_inclusion){
  
  df_nb_positions_per_S_region <- df_OR_AA_mut %>% 
    filter(gene == 'S', threshold_for_inclusion == threshold_inclusion) %>% filter(p_value_BH < 0.05) %>% 
    group_by(mutation) %>% 
    mutate(position = as.integer(substr(mutation, start = 1, stop = nchar(mutation) - 1))) %>% 
    mutate(region = case_when(position >= 13 & position <= 302 ~ 'NTD',
                              position >= 319 & position <= 541 ~ 'RBD',
                              position >= 13 & position <= 685 ~ 'Other S1',
                              position >= 686 & position <= 1273 ~ 'S2',
                              T ~ 'Other S')) %>% 
    group_by(region) %>% summarise(n_pos = n()) %>% 
    mutate(region = factor(region, levels = c( 'Other S1', 'RBD', 'NTD', 'S2', 'Other S')),
           region_plot = case_when(region %in% c('NTD', 'RBD', 'Other S1') ~ 'S1',
                                   region == 'S2' ~ 'S2',
                                   T ~ 'Other'),
           region_plot = factor(region_plot, levels = c('S1', 'S2', 'Other')))
  
  barplot_nb_enriched_positions_by_gene <- df_nb_positions_per_S_region %>% 
    ggplot(aes(x = region_plot, y = n_pos, fill = region, group = region)) +
    geom_bar(stat = 'identity') +
    scale_x_discrete(name = 'Subunit') +
    scale_y_continuous(name = 'Number of enriched mutations', 
                       expand = expansion(mult = c(0., 0.1))) +
    scale_fill_manual(breaks = c('NTD', 'RBD'), values = c(palette_spike['NTD'], palette_spike['RBD']), name = 'Domain', 
                      na.value = palette_genes['S']) +
    theme_classic() +
    theme(axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1., vjust = 1.))
  
  return(barplot_nb_enriched_positions_by_gene)
}
get_barplot_nb_enriched_mutations_within_ORF1ab <- function(df_OR_AA_mut, threshold_inclusion){
  
  df_nb_positions_per_NSPs <- df_OR_AA_mut %>% 
    filter(gene %in% c('ORF1a', 'ORF1b'), threshold_for_inclusion == threshold_inclusion) %>% filter(p_value_BH < 0.05) %>% 
    group_by(mutation) %>% 
    mutate(position = as.integer(substr(mutation, start = 1, stop = nchar(mutation) - 1))) %>% 
    mutate(NSP = case_when((gene == 'ORF1a') & (position < 181) ~ 'NSP1',
                           (gene == 'ORF1a') & (position < 819) ~ 'NSP2',
                           (gene == 'ORF1a') & (position < 2764) ~ 'NSP3',
                           (gene == 'ORF1a') & (position < 3264) ~ 'NSP4',
                           (gene == 'ORF1a') & (position < 3570) ~ 'NSP5',
                           (gene == 'ORF1a') & (position < 3860) ~ 'NSP6',
                           (gene == 'ORF1a') & (position < 3943) ~ 'NSP7',
                           (gene == 'ORF1a') & (position < 4141) ~ 'NSP8',
                           (gene == 'ORF1a') & (position < 4254) ~ 'NSP9',
                           (gene == 'ORF1a') & (position < 4393) ~ 'NSP10',
                           (gene == 'ORF1a') ~ 'NSP12',
                           (gene == 'ORF1b') & (position < 924) ~ 'NSP12',
                           (gene == 'ORF1b') & (position < 1525) ~ 'NSP13',
                           (gene == 'ORF1b') & (position < 2052) ~ 'NSP14',
                           (gene == 'ORF1b') & (position < 2398) ~ 'NSP15',
                           (gene == 'ORF1b') & (position <= 2696) ~ 'NSP16',
                           T ~ 'Other'
                           ))  %>% 
    group_by(NSP, gene) %>% summarise(n_pos = n()) %>% 
    mutate(NSP = factor(NSP, levels = paste0('NSP', 16:1)))
  
    barplot_nb_enriched_positions_by_gene <- df_nb_positions_per_NSPs %>% 
    ggplot(aes(x = gene, y = n_pos, fill = NSP, group = NSP)) +
    geom_bar(stat = 'identity') +
    scale_x_discrete(name = 'Gene') +
    scale_y_continuous(name = 'Number of enriched mutations', 
                       expand = expansion(mult = c(0., 0.1)),
                       breaks = seq(0, 50, 5)) +
    scale_fill_manual(values = palette_NSPs, limits = rev) +
    theme_classic() +
    theme(axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1., vjust = 1.))
  
  return(barplot_nb_enriched_positions_by_gene)
}

barplot_nb_enriched_mutations_by_gene <- get_barplot_nb_enriched_mutations(df_OR_AA_mut, 2)
barplot_nb_enriched_mutations_within_S <- get_barplot_nb_enriched_mutations_within_S(df_OR_AA_mut, 2)
barplot_nb_enriched_mutations_within_ORF1ab <- get_barplot_nb_enriched_mutations_within_ORF1ab(df_OR_AA_mut, 2)

pdf('../figure_mutation_enrichment/barplot_nb_enriched_mutations_by_gene_saltational.pdf',
    height = 2.8, width = 4.)
plot(barplot_nb_enriched_mutations_by_gene)
dev.off()
png('../figure_mutation_enrichment/barplot_nb_enriched_mutations_by_gene_saltational.png',
    height = 2.8, width = 4., res = 350, units = 'in')
plot(barplot_nb_enriched_mutations_by_gene)
dev.off()

pdf('../figure_mutation_enrichment/barplot_nb_enriched_mutations_S_by_subunit_saltational.pdf',
    height = 2.8, width = 2.8)
plot(barplot_nb_enriched_mutations_within_S)
dev.off()
png('../figure_mutation_enrichment/barplot_nb_enriched_mutations_S_by_subunit_saltational.png',
    height = 2.8, width = 2.8, res = 350, units = 'in')
plot(barplot_nb_enriched_mutations_within_S)
dev.off()

pdf('../figure_mutation_enrichment/barplot_nb_enriched_mutations_ORF1ab_by_NSP_saltational.pdf',
    height = 3.8, width = 2.8)
plot(barplot_nb_enriched_mutations_within_ORF1ab)
dev.off()
png('../figure_mutation_enrichment/barplot_nb_enriched_mutations_ORF1ab_by_NSP_saltational.png',
    height = 3.8, width = 2.8, res = 350, units = 'in')
plot(barplot_nb_enriched_mutations_within_ORF1ab)
dev.off()


######### Get dataframe with enriched AA muts
df_OR_AA_mut %>% 
  filter(threshold_for_inclusion == 2) %>% filter(p_value_BH < 0.05, OR > 1.) %>% 
  arrange(- n_mut_of_interest_adapt_branches) %>% 
  mutate(OR_char = paste0(round(OR, 2), ' (', round(OR_lower, 2), '-', round(OR_upper, 2), ')')) %>% 
  select(gene, mutation, OR_char, n_mut_of_interest_adapt_branches) %>% 
  write_csv('../figure_mutation_enrichment/summary_enriched_mutations.csv')
  