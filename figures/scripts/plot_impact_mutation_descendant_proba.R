library(tidyverse)
library(scales)
library(binom)
source('utils_descendants_analysis.R')
source('config_plot.R')

# Load saltational branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', header = F) %>%
  unlist() %>% as.character()

name_saltational_branches_less_100_desc <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches) %>% 
  filter(n_terminal_descendants <= 100) %>% 
  select(strain) %>% unlist() %>% as.character()

name_saltational_branches_less_1000_desc <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches) %>% 
  filter(n_terminal_descendants <= 1000) %>% 
  select(strain) %>% unlist() %>% as.character()


# Load descendants information
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(has_descendants = n_terminal_descendants > 1) %>% 
  select(strain, has_descendants)

inferred_country <- read_tsv('../../results/bte/inferred_countries.tsv') %>% select(- location) %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(spread_to_multiple_countries = (inferred_location == 'Multiple descendant locations')) %>% 
  select(strain, spread_to_multiple_countries)
  
inferred_subregion <- read_tsv('../../results/bte/inferred_subregions.tsv') %>% select(- location) %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(spread_to_multiple_subregions = (inferred_location == 'Multiple descendant locations')) %>% 
  select(strain, spread_to_multiple_subregions)

inferred_continent <- read_tsv('../../results/bte/inferred_continents.tsv') %>% select(- location) %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(spread_to_multiple_continents = (inferred_location == 'Multiple descendant locations')) %>% 
  select(strain, spread_to_multiple_continents)

df_bool_descendants_across_scale <- df_descendants %>% 
  left_join(inferred_country) %>% 
  left_join(inferred_subregion) %>% 
  left_join(inferred_continent)

# Load mutation enrichment information
enriched_AA_mut <- read_tsv('../../results/mutation_enrichment/df_AA_mut_enrichment_saltational_opportunity_CV_1.8_BF_20.tsv') %>% 
  filter(threshold_for_inclusion == 2) %>% filter(p_value_BH < 0.05) %>% select(gene, mutation, OR) %>% 
  mutate(mut_full = paste0(gene, ':', mutation))

AA_mut_observed_salt_branches <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  filter(aa_initial != aa_new, node_id %in% name_saltational_branches) %>% 
  mutate(mut_full = paste0(gene, ':', aa_muts_position, aa_new)) %>% 
  mutate(mutation = paste0(aa_muts_position, aa_new)) %>% 
  select(gene, mutation, mut_full) %>% distinct()

# Load AA mutation information
df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  rename(strain = node_id) %>% 
  filter(is_coded, strain %in% name_saltational_branches) %>% 
  mutate(is_synonymous = (aa_initial == aa_new),
         aa_muts_direction = paste0(aa_muts_position, aa_new)) %>% 
  filter(! is_synonymous) %>% 
  left_join(df_bool_descendants_across_scale)

##
get_prop_branches_having_muts_with_desc <- function(df_aa_muts, 
                                                    name_saltational_branches,
                                                    df_AA_mut_for_prop_comput){
  
  df_aa_muts_crop <- df_aa_muts %>% filter(strain %in% name_saltational_branches)
  df_prop_with_descendants_by_mut <- Reduce('bind_rows', lapply(1:nrow(df_AA_mut_for_prop_comput), FUN = function(i_row){
    tmp_df_has_mut <- df_aa_muts_crop %>% 
      filter(gene == df_AA_mut_for_prop_comput$gene[i_row],
             aa_muts_direction == df_AA_mut_for_prop_comput$mutation[i_row]) %>% 
      mutate(has_mut = T) %>% 
      select(strain, has_mut, has_descendants)
  
    n_with_mut_with_descendants <- sum((tmp_df_has_mut$has_descendants) & (tmp_df_has_mut$has_mut))
    n_with_mut_without_descendants <- sum((!tmp_df_has_mut$has_descendants) & (tmp_df_has_mut$has_mut))
    n_without_mut_with_descendants <- sum((tmp_df_has_mut$has_descendants) & (!tmp_df_has_mut$has_mut))
    n_without_mut_without_descendants <- sum((!tmp_df_has_mut$has_descendants) & (!tmp_df_has_mut$has_mut))
  

    tibble(gene = df_AA_mut_for_prop_comput$gene[i_row],
           aa_muts_direction = df_AA_mut_for_prop_comput$mutation[i_row],
           n_with_mut_with_descendants = n_with_mut_with_descendants,
           n_with_mut_without_descendants = n_with_mut_without_descendants,
           n_without_mut_with_descendants = n_without_mut_with_descendants,
           n_without_mut_without_descendants = n_without_mut_without_descendants)
  })) %>% 
    mutate(scale = 'Multiple descendants')
  
  df_prop_with_descendants_multiple_countries_by_mut <- Reduce('bind_rows', lapply(1:nrow(df_AA_mut_for_prop_comput), FUN = function(i_row){
    tmp_df_has_mut <- df_aa_muts_crop %>% 
      filter(gene == df_AA_mut_for_prop_comput$gene[i_row],
             aa_muts_direction == df_AA_mut_for_prop_comput$mutation[i_row]) %>% 
      mutate(has_mut = T) %>% 
      select(strain, has_mut, spread_to_multiple_countries)
    
    n_with_mut_with_descendants <- sum((tmp_df_has_mut$spread_to_multiple_countries) & (tmp_df_has_mut$has_mut))
    n_with_mut_without_descendants <- sum((!tmp_df_has_mut$spread_to_multiple_countries) & (tmp_df_has_mut$has_mut))
    n_without_mut_with_descendants <- sum((tmp_df_has_mut$spread_to_multiple_countries) & (!tmp_df_has_mut$has_mut))
    n_without_mut_without_descendants <- sum((!tmp_df_has_mut$spread_to_multiple_countries) & (!tmp_df_has_mut$has_mut))
    
    
    tibble(gene = df_AA_mut_for_prop_comput$gene[i_row],
           aa_muts_direction = df_AA_mut_for_prop_comput$mutation[i_row],
           n_with_mut_with_descendants = n_with_mut_with_descendants,
           n_with_mut_without_descendants = n_with_mut_without_descendants,
           n_without_mut_with_descendants = n_without_mut_with_descendants,
           n_without_mut_without_descendants = n_without_mut_without_descendants)
  })) %>% 
    mutate(scale = 'Multiple countries')
  
  df_prop_with_descendants_multiple_subregions_by_mut <- Reduce('bind_rows', lapply(1:nrow(df_AA_mut_for_prop_comput), FUN = function(i_row){
    tmp_df_has_mut <- df_aa_muts_crop %>% 
      filter(gene == df_AA_mut_for_prop_comput$gene[i_row],
             aa_muts_direction == df_AA_mut_for_prop_comput$mutation[i_row]) %>% 
      mutate(has_mut = T) %>% 
      select(strain, has_mut, spread_to_multiple_subregions)
    
    n_with_mut_with_descendants <- sum((tmp_df_has_mut$spread_to_multiple_subregions) & (tmp_df_has_mut$has_mut))
    n_with_mut_without_descendants <- sum((!tmp_df_has_mut$spread_to_multiple_subregions) & (tmp_df_has_mut$has_mut))
    n_without_mut_with_descendants <- sum((tmp_df_has_mut$spread_to_multiple_subregions) & (!tmp_df_has_mut$has_mut))
    n_without_mut_without_descendants <- sum((!tmp_df_has_mut$spread_to_multiple_subregions) & (!tmp_df_has_mut$has_mut))
    
    
    tibble(gene = df_AA_mut_for_prop_comput$gene[i_row],
           aa_muts_direction = df_AA_mut_for_prop_comput$mutation[i_row],
           n_with_mut_with_descendants = n_with_mut_with_descendants,
           n_with_mut_without_descendants = n_with_mut_without_descendants,
           n_without_mut_with_descendants = n_without_mut_with_descendants,
           n_without_mut_without_descendants = n_without_mut_without_descendants)
  })) %>% 
    mutate(scale = 'Multiple subregions')
  
  df_prop_with_descendants_multiple_continents_by_mut <- Reduce('bind_rows', lapply(1:nrow(df_AA_mut_for_prop_comput), FUN = function(i_row){
    tmp_df_has_mut <- df_aa_muts %>% 
      filter(gene == df_AA_mut_for_prop_comput$gene[i_row],
             aa_muts_direction == df_AA_mut_for_prop_comput$mutation[i_row]) %>% 
      mutate(has_mut = T) %>% 
      select(strain, has_mut, spread_to_multiple_continents)
    
    n_with_mut_with_descendants <- sum((tmp_df_has_mut$spread_to_multiple_continents) & (tmp_df_has_mut$has_mut))
    n_with_mut_without_descendants <- sum((!tmp_df_has_mut$spread_to_multiple_continents) & (tmp_df_has_mut$has_mut))
    n_without_mut_with_descendants <- sum((tmp_df_has_mut$spread_to_multiple_continents) & (!tmp_df_has_mut$has_mut))
    n_without_mut_without_descendants <- sum((!tmp_df_has_mut$spread_to_multiple_continents) & (!tmp_df_has_mut$has_mut))
    
    
    tibble(gene = df_AA_mut_for_prop_comput$gene[i_row],
           aa_muts_direction = df_AA_mut_for_prop_comput$mutation[i_row],
           n_with_mut_with_descendants = n_with_mut_with_descendants,
           n_with_mut_without_descendants = n_with_mut_without_descendants,
           n_without_mut_with_descendants = n_without_mut_with_descendants,
           n_without_mut_without_descendants = n_without_mut_without_descendants)
  })) %>% 
    mutate(scale = 'Multiple continents')
  
  return(
    bind_rows(
      df_prop_with_descendants_by_mut,
      df_prop_with_descendants_multiple_countries_by_mut, 
      df_prop_with_descendants_multiple_subregions_by_mut,
      df_prop_with_descendants_multiple_continents_by_mut
    )
  )
}

df_prop_branches_with_mut_among_enriched_muts <- get_prop_branches_having_muts_with_desc(df_aa_muts, 
                                                                                         name_saltational_branches,
                                                                                         enriched_AA_mut) %>% 
  mutate(full_mut = paste0(gene, ':', aa_muts_direction))

df_prop_branches_with_mut_among_enriched_muts_branches_less_100_muts <- get_prop_branches_having_muts_with_desc(df_aa_muts, 
                                                                                                                name_saltational_branches_less_100_desc,
                                                                                                                enriched_AA_mut) %>% 
  mutate(full_mut = paste0(gene, ':', aa_muts_direction))

df_prop_branches_with_mut_among_enriched_muts_branches_less_1000_muts <- get_prop_branches_having_muts_with_desc(df_aa_muts, 
                                                                                                                name_saltational_branches_less_1000_desc,
                                                                                                                enriched_AA_mut) %>% 
  mutate(full_mut = paste0(gene, ':', aa_muts_direction))



## Probably should only depict what happens for the most enriched mutations
plt_prop_branches_with_desc_func_nb_branches_with_mut <- df_prop_branches_with_mut_among_enriched_muts %>% 
  filter(scale == 'Multiple descendants') %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  ggplot(aes(x = n_with_mut_with_descendants + n_with_mut_without_descendants, 
             y = n_with_mut_with_descendants/(n_with_mut_with_descendants + n_with_mut_without_descendants), 
             colour = gene)) +
  geom_point() +
  ggrepel::geom_text_repel(data = df_prop_branches_with_mut_among_enriched_muts %>% 
                             filter(scale == 'Multiple descendants') %>% 
                             left_join(enriched_AA_mut, by = c('gene', 'aa_muts_direction' = 'mutation')) %>% 
                             filter(OR > 5., n_with_mut_with_descendants + n_with_mut_without_descendants >= 5),
                           aes(label = aa_muts_direction), 
                           min.segment.length = 0.,
                           max.overlaps = Inf) +
  scale_x_continuous(trans = 'log',
                     breaks = c(1:10, 15, 20, 30, 40),
                     name = 'Number of saltational branches with mutation') +
  scale_y_continuous(name = 'Proportion of branches with descendants', 
                     breaks = seq(0., 1., 0.2)) +
  scale_colour_manual(values = palette_genes,
                      name = 'Gene') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        axis.text = element_text(size = 11)) +
  guides(colour = guide_legend(override.aes = list(label = '')))

plot(plt_prop_branches_with_desc_func_nb_branches_with_mut)

pdf('../figure_descendant_analysis/prop_branches_with_desc_func_nb_branches_with_mut.pdf',
    height = 5.5, width = 8.5)
plot(plt_prop_branches_with_desc_func_nb_branches_with_mut)
dev.off()  
png('../figure_descendant_analysis/prop_branches_with_desc_func_nb_branches_with_mut.png',
    height = 5.5, width = 8.5, res = 350, units = 'in')
plot(plt_prop_branches_with_desc_func_nb_branches_with_mut)
dev.off()  

#plt_prop_branches_with_desc_func_nb_branches_with_mut <- 
value_zero_desc <- 1e-2

test_plot <- df_prop_branches_with_mut_among_enriched_muts %>% 
  filter(scale == 'Multiple descendants') %>% 
  mutate(y_axis_to_plot = ifelse(n_with_mut_with_descendants == 0, value_zero_desc, n_with_mut_with_descendants)) %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  ggplot(aes(x = n_with_mut_with_descendants + n_with_mut_without_descendants, 
             y = y_axis_to_plot)) +
  geom_point(aes(colour = gene)) +
  ggrepel::geom_text_repel(data = df_prop_branches_with_mut_among_enriched_muts %>% 
                             mutate(y_axis_to_plot = ifelse(n_with_mut_with_descendants == 0, value_zero_desc, n_with_mut_with_descendants)) %>% 
                             filter(scale == 'Multiple descendants') %>% 
                             left_join(enriched_AA_mut, by = c('gene', 'aa_muts_direction' = 'mutation')) %>% 
                             filter(OR > 5., n_with_mut_with_descendants + n_with_mut_without_descendants >= 5),
                           aes(label = aa_muts_direction, colour = gene), 
                           min.segment.length = 0.,
                           max.overlaps = Inf) +
  facet_grid((n_with_mut_with_descendants <= 0) ~ ., scales = 'free', space = 'free_y') +
  scale_x_continuous(trans = 'log',
                     minor_breaks = minor_breaks_log(smallest = 1, detail = 1),
                     breaks = c(1:10, 12, 15, 20, 30, 40),
                     name = 'Number of saltational branches with mutation') +
  scale_y_continuous(name = 'Number of saltational branches with mutation with descendants',
                     trans = 'log',
                     minor_breaks = minor_breaks_log(smallest = 1, detail = 1),
                     expand = expansion(mult = c(0.05, 0.05), 
                                        add = c(0.2, 0.2)),
                     breaks = c(value_zero_desc, 1:10, 12, 15, 20, 30, 40),
                     labels = c('0', 1:10, 12, 15, 20, 30, 40)
                     ) +
  scale_colour_manual(values = palette_genes,
                      name = 'Gene') +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid = element_blank()
        ) +
  guides(colour = guide_legend(override.aes = list(label = '')))

pdf('../figure_descendant_analysis/nb_branches_with_desc_func_nb_branches_with_mut.pdf',
    height = 5., width = 7.)
plot(test_plot)
dev.off()






##### Compare with global trends
recombinants_clades_to_drop <- c('24F (XEC)', '25G (XFC)', '22F (XBB)', '25C (XFG)')
name_clades_after_mid_2024 <- c('24H (LF.7)', '25A (LP.8.1)', '25B (NB.1.8.1)', '25E (PY.1)', '25F (NW.1.2)') # Some are missing
df_mut_counts_by_nextstrain_clade <- read_tsv('../../results/mutation_enrichment/df_aa_freq_combined_CV_1.8_BF_20.tsv') %>% 
  mutate(prop_AA = count / clade_size, full_mut = paste0(gene, ':', pos, AA)) %>% 
  filter(AA != '-', ! clade %in% recombinants_clades_to_drop, is_lineage_defining | (is_within_clade_mut & (prop_AA > 0.8))) %>% 
  group_by(full_mut) %>% summarise(n_clades = n())

df_mut_counts_by_nextstrain_clade_post_mid_2024 <- read_tsv('../../results/mutation_enrichment/df_aa_freq_combined_CV_1.8_BF_20.tsv') %>% 
  mutate(prop_AA = count / clade_size, full_mut = paste0(gene, ':', pos, AA)) %>% 
  filter(AA != '-', clade %in% name_clades_after_mid_2024, is_lineage_defining | (is_within_clade_mut & (prop_AA > 0.8))) %>% 
  group_by(full_mut) %>% summarise(n_clades = n())


bootstrap_mean_function <-  function(x, B = 2000) {
  boot_means <- replicate(B, mean(sample(x, replace = TRUE), na.rm = TRUE))
  
  tibble(
    mean = mean(x, na.rm = TRUE),
    lwr = quantile(boot_means, 0.025, na.rm = TRUE),
    upr = quantile(boot_means, 0.975, na.rm = TRUE)
  )
}

df_for_plot <- df_prop_branches_with_mut_among_enriched_muts %>% 
  filter(scale == 'Multiple descendants') %>% 
  left_join(df_mut_counts_by_nextstrain_clade %>% 
              filter(full_mut %in% df_prop_branches_with_mut_among_enriched_muts$full_mut)) %>% 
  mutate(n_clades = replace_na(n_clades, 0)) %>% 
  mutate(quantile_prop = ntile(n_with_mut_with_descendants/(n_with_mut_with_descendants + n_with_mut_without_descendants), 5)) %>% 
  group_by(quantile_prop) %>% 
  reframe(
    mean_number_clades = mean(n_clades),
    mean_prop_within_quantile = mean(n_with_mut_with_descendants/(n_with_mut_with_descendants + n_with_mut_without_descendants), 
                                     na.rm = TRUE),
    bootstrap_mean_function(n_clades)
  )

df_for_plot_crop_100 <- df_prop_branches_with_mut_among_enriched_muts_branches_less_100_muts %>% 
  filter(scale == 'Multiple descendants') %>% 
  left_join(df_mut_counts_by_nextstrain_clade %>% 
              filter(full_mut %in% df_prop_branches_with_mut_among_enriched_muts_branches_less_100_muts$full_mut)) %>% 
  mutate(n_clades = replace_na(n_clades, 0)) %>% 
  mutate(quantile_prop = ntile(n_with_mut_with_descendants/(n_with_mut_with_descendants + n_with_mut_without_descendants), 5)) %>% 
  group_by(quantile_prop) %>% 
  reframe(
    mean_number_clades = mean(n_clades),
    mean_prop_within_quantile = mean(n_with_mut_with_descendants/(n_with_mut_with_descendants + n_with_mut_without_descendants), 
                                     na.rm = TRUE),
    bootstrap_mean_function(n_clades)
  )

df_for_plot_crop_1000 <- df_prop_branches_with_mut_among_enriched_muts_branches_less_1000_muts %>% 
  filter(scale == 'Multiple descendants') %>% 
  left_join(df_mut_counts_by_nextstrain_clade %>% 
              filter(full_mut %in% df_prop_branches_with_mut_among_enriched_muts_branches_less_1000_muts$full_mut)) %>% 
  mutate(n_clades = replace_na(n_clades, 0)) %>% 
  mutate(quantile_prop = ntile(n_with_mut_with_descendants/(n_with_mut_with_descendants + n_with_mut_without_descendants), 5)) %>% 
  group_by(quantile_prop) %>% 
  reframe(
    mean_number_clades = mean(n_clades),
    mean_prop_within_quantile = mean(n_with_mut_with_descendants/(n_with_mut_with_descendants + n_with_mut_without_descendants), 
                                     na.rm = TRUE),
    bootstrap_mean_function(n_clades)
  )


plt_mean_nb_clades_in_which_mutations_fixed <- 
  bind_rows(df_for_plot %>% mutate(type = 'All branches'), 
            df_for_plot_crop_100 %>% mutate(type = 'Less than 100 desc.'), 
            df_for_plot_crop_1000 %>% mutate(type = 'Less than 1000 desc.')) %>% 
  mutate(type = factor(type, levels = c('All branches', 'Less than 1000 desc.', 'Less than 100 desc.'))) %>% 
  ggplot(aes(x = quantile_prop, y = mean_number_clades, colour = type)) +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = lwr, ymax = upr), position = position_dodge(0.4)) +
  scale_x_continuous(name = 'Quintile of proportion of saltational\nbranches with descendants') +
  scale_y_continuous(name = 'Mean number of clades\nin which mutation fixed',
                     breaks = seq(0, 2., 0.2)) + 
  scale_colour_manual(breaks = c('All branches', 'Less than 1000 desc.', 'Less than 100 desc.'),
                      values = c('darkred', 'firebrick2', 'goldenrod1'),
                      name = '') +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        panel.grid.minor = element_blank(),
        legend.position = 'right')

png('../figure_descendant_analysis/nb_clades_where_mut_fixed_func_prop_of_salt_branches_per_desc.png',
    height = 2.8, width = 5., res = 350, units = 'in')
plot(plt_mean_nb_clades_in_which_mutations_fixed)
dev.off()
pdf('../figure_descendant_analysis/nb_clades_where_mut_fixed_func_prop_of_salt_branches_per_desc.pdf',
    height = 2.8, width = 5.)
plot(plt_mean_nb_clades_in_which_mutations_fixed)
dev.off()
