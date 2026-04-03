library(tidyverse)
library(ggsignif)
source('config_plot.R')

vec_genes <- names(palette_genes)

# Load adaptive branches names
name_adaptive_branches <- read.table('../../results/adaptive_branches/concat_name_adaptive_branches_CV_1.7_BF_31.6227766.tsv', 
                                     header = F) %>% unlist() %>% as.character()

# Load number of mutations by branch
mut_count_genome <- read_tsv('../../results/input_model/data_viridian_Genome_aggregated.tsv')
mut_count_S <- read_tsv('../../results/input_model/data_viridian_S.tsv')

# Load AA mutation information
df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  rename(strain = node_id) %>% filter(is_coded) %>% 
  mutate(is_adaptive = strain %in% name_adaptive_branches) %>% 
  group_by(gene, is_adaptive) %>% 
  summarise(n_tot_mut = n(), n_non_syn = sum(aa_initial != aa_new))

# Compare number of mutations by branch type
get_boxplot_mut_count_per_branch_type <- function(mut_count_genome, name_adaptive_branches){
  
  df_for_boxplot_tot_muts <- mut_count_genome %>% 
    mutate(is_adaptive = strain %in% name_adaptive_branches) %>% 
    mutate(tot_muts = synonymous_muts + non_synonymous_muts) %>% 
    group_by(is_adaptive = strain %in% name_adaptive_branches) %>% 
    summarise(quantile_2.5 = quantile(tot_muts, 0.025),
              quantile_25 = quantile(tot_muts, 0.25),
              quantile_50 = quantile(tot_muts, 0.5),
              quantile_75 = quantile(tot_muts, 0.75),
              quantile_97.5 = quantile(tot_muts, 0.975)) %>% 
    mutate(type = 'All')
  
  df_for_boxplot_syn_muts <- mut_count_genome %>% 
    mutate(is_adaptive = strain %in% name_adaptive_branches) %>% 
    mutate(tot_muts = synonymous_muts) %>% 
    group_by(is_adaptive = strain %in% name_adaptive_branches) %>% 
    summarise(quantile_2.5 = quantile(tot_muts, 0.025),
              quantile_25 = quantile(tot_muts, 0.25),
              quantile_50 = quantile(tot_muts, 0.5),
              quantile_75 = quantile(tot_muts, 0.75),
              quantile_97.5 = quantile(tot_muts, 0.975)) %>% 
    mutate(type = 'Synonymous')
  
  df_for_boxplot_non_syn_muts <- mut_count_genome %>% 
    mutate(is_adaptive = strain %in% name_adaptive_branches) %>% 
    mutate(tot_muts = non_synonymous_muts) %>% 
    group_by(is_adaptive = strain %in% name_adaptive_branches) %>% 
    summarise(quantile_2.5 = quantile(tot_muts, 0.025),
              quantile_25 = quantile(tot_muts, 0.25),
              quantile_50 = quantile(tot_muts, 0.5),
              quantile_75 = quantile(tot_muts, 0.75),
              quantile_97.5 = quantile(tot_muts, 0.975)) %>% 
    mutate(type = 'Non-synonymous')
  
  # Wilcoxon test to compare results
  tmp_df_tot_muts <- mut_count_genome %>% 
    mutate(is_adaptive = strain %in% name_adaptive_branches,
           mut_tot = synonymous_muts + non_synonymous_muts)
  w_test_tot_muts <- wilcox.test(x = tmp_df_tot_muts$mut_tot[tmp_df_tot_muts$is_adaptive],
                                 y = tmp_df_tot_muts$mut_tot[! tmp_df_tot_muts$is_adaptive])
  w_test_syn_muts <- wilcox.test(x = tmp_df_tot_muts$synonymous_muts[tmp_df_tot_muts$is_adaptive],
                                 y = tmp_df_tot_muts$synonymous_muts[! tmp_df_tot_muts$is_adaptive])
  w_test_non_syn_muts <- wilcox.test(x = tmp_df_tot_muts$non_synonymous_muts[tmp_df_tot_muts$is_adaptive],
                                     y = tmp_df_tot_muts$non_synonymous_muts[! tmp_df_tot_muts$is_adaptive])
  
  p_to_stars <- function(p) {
    case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    )
  }
  
  signif_df <- tibble(
    type = c("All", "Synonymous", "Non-synonymous"),
    xmin = c("All", "Synonymous", "Non-synonymous"),
    xmax = c("All", "Synonymous", "Non-synonymous"),
    y_position = c(
      max(df_for_boxplot_tot_muts$quantile_97.5) + 2,
      max(df_for_boxplot_syn_muts$quantile_97.5) + 2,
      max(df_for_boxplot_non_syn_muts$quantile_97.5) + 2
    ),
    annotations = c(
      p_to_stars(w_test_tot_muts$p.value),
      p_to_stars(w_test_syn_muts$p.value),
      p_to_stars(w_test_non_syn_muts$p.value)
    )
  ) %>% 
    mutate(x_num = c(3, 1, 2))
  
  
  boxplot_mut_count_func_branch_type <- bind_rows(df_for_boxplot_tot_muts, df_for_boxplot_syn_muts, df_for_boxplot_non_syn_muts) %>% 
    ggplot(aes(x = type, group = interaction(type, is_adaptive), 
               fill = is_adaptive, colour = is_adaptive)) +
    geom_boxplot(aes(ymin = quantile_2.5, lower = quantile_25, middle = quantile_50, upper = quantile_75, ymax = quantile_97.5),
                 stat = 'identity', 
                 position = position_dodge(0.7),
                 width = 0.65, alpha = 0.2) +
    geom_signif(data = signif_df,
      aes(xmin = xmin, xmax = xmax, annotations = annotations, y_position = y_position),
      manual = TRUE,
      tip_length = 0.,
      textsize = 4,
      vjust = 0.2,
      inherit.aes = FALSE
    ) +
    geom_segment(data = signif_df, aes(x = x_num - 0.15, xend = x_num + 0.15, y = y_position, yend = y_position),
                 inherit.aes = FALSE) + 
    scale_x_discrete(name = '', limits = rev) +
    scale_y_continuous(name = 'Number of mutations', expand = expansion(mult = c(0.01, 0.05)),
                       breaks = seq(0., 25, 5)) + 
    scale_fill_manual(name = '', breaks = c(T, F), values = c(col_adaptive, col_non_adaptive),
                      labels = c('Adaptive', 'Non adaptive')) +
    scale_colour_manual(name = '', breaks = c(T, F), values = c(col_adaptive, col_non_adaptive),
                        labels = c('Adaptive', 'Non adaptive')) +
    theme_classic() +
    theme(legend.position = 'inside',
          legend.background = element_blank(),
          legend.position.inside = c(0.8, 0.3),
          axis.text = element_text(size = 11)) +
    coord_flip()
 
  return(boxplot_mut_count_func_branch_type) 
}

boxplot_mut_count_by_branch_type <- get_boxplot_mut_count_per_branch_type(mut_count_genome, name_adaptive_branches)

plot(boxplot_mut_count_by_branch_type)

pdf('../figure_adaptive_branches/boxplot_mut_counts_per_branch_type.pdf',
    height = 2.2, width = 5.)
plot(boxplot_mut_count_by_branch_type)
dev.off()
png('../figure_adaptive_branches/boxplot_mut_counts_per_branch_type.png',
    height = 2.2, width = 5., res = 350, units = 'in')
plot(boxplot_mut_count_by_branch_type)
dev.off()


# Distribution of mutations oer gene
get_mut_frequency_by_gene_by_branch_type <- function(df_aa_muts){
  hist_mut_freq_by_branch_type <- df_aa_muts %>% 
    group_by(is_adaptive) %>% 
    mutate(prop_tot_mut = n_tot_mut / sum(n_tot_mut),
           prop_non_syn = n_non_syn / sum(n_non_syn)) %>% 
    mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
    ggplot(aes(x = gene, y = prop_tot_mut, fill = is_adaptive)) +
    geom_bar(stat = 'identity', width = 0.8,
             position = position_dodge(0.8)) + 
    scale_x_discrete(name = 'Gene') +
    scale_y_continuous(name = 'Proportion of mutations',
                       expand = expansion(mult = c(0.01, 0.05))) +
    scale_fill_manual(name = '', breaks = c(T, F), values = c(col_adaptive, col_non_adaptive),
                      labels = c('Adaptive', 'Non adaptive')) +
    theme_classic() +
    theme(axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, vjust = 1., hjust = 1.),
          legend.position = 'inside',
          legend.position.inside = c(0.8, 0.8),
          legend.background = element_blank())
  
  return(hist_mut_freq_by_branch_type)
}

hist_mut_freq_by_branch_type <- get_mut_frequency_by_gene_by_branch_type(df_aa_muts)

pdf('../figure_adaptive_branches/mut_frequency_by_gene_by_branch_type.pdf',
    height = 2.8, width = 4.)
plot(hist_mut_freq_by_branch_type)
dev.off()
png('../figure_adaptive_branches/mut_frequency_by_gene_by_branch_type.png',
    height = 2.8, width = 4., res = 350, units = 'in')
plot(hist_mut_freq_by_branch_type)
dev.off()
