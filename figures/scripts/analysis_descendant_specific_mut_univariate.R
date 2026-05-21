library(tidyverse)
library(scales)
library(binom)
source('utils_descendants_analysis.R')
source('config_plot.R')

# Load saltational branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', header = F) %>%
  unlist() %>% as.character()

# Load descendants information
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(has_descendants = n_terminal_descendants > 1) %>% 
  select(strain, has_descendants)

df_descendant_counts <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(has_descendants = n_terminal_descendants > 1)

# Get proportion of saltational branches with descendants
n_salt_with_descendants <- sum(df_descendant_counts$has_descendants)
n_salt <- length(df_descendant_counts$has_descendants)
prop_salt_with_descendants <- n_salt_with_descendants/n_salt
CI_prop_salt_with_descendants <- prop.test(x = n_salt_with_descendants, n = n_salt)$conf.int

  
# Load mutation enrichment information
df_OR_AA_mut <- read_tsv('../../results/mutation_enrichment/df_AA_mut_enrichment_saltational_opportunity_CV_1.8_BF_20.tsv') %>% 
  filter(threshold_for_inclusion == 2)

enriched_AA_mut <- df_OR_AA_mut %>% filter(p_value_BH < 0.05) %>% select(gene, mutation) %>% 
  mutate(mut_full = paste0(gene, ':', mutation))

# Load AA mutation information
df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  rename(strain = node_id) %>% 
  filter(is_coded, strain %in% name_saltational_branches) %>% 
  mutate(is_synonymous = (aa_initial == aa_new),
         aa_muts_direction = paste0(aa_muts_position, aa_new)) %>% 
  filter(! is_synonymous)

get_OR_having_descendants_having_muts <- function(df_aa_muts, 
                                                  df_descendants,
                                                  threshold_mutation_occurence_inclusion = 4){
  
  df_muts_to_include <- df_aa_muts %>% 
    group_by(gene, aa_muts_direction, is_synonymous) %>% 
    summarise(n_occ = n()) %>% 
    filter(n_occ > threshold_mutation_occurence_inclusion)
  
  df_OR_having_descendants_by_mut <- Reduce('bind_rows', lapply(1:nrow(df_muts_to_include), FUN = function(i_row){
    tmp_df_has_mut <- df_aa_muts %>% 
      filter(gene == df_muts_to_include$gene[i_row],
             aa_muts_direction == df_muts_to_include$aa_muts_direction[i_row],
             is_synonymous == df_muts_to_include$is_synonymous[i_row]) %>% 
      mutate(has_mut = T) %>% 
      select(strain, has_mut)
    
    tmp_df_descendants_muts <- df_descendants %>% 
      left_join(tmp_df_has_mut, by = 'strain') %>% 
      mutate(has_mut = replace_na(has_mut, F))
    
    n_with_mut_with_descendants <- sum((tmp_df_descendants_muts$has_descendants) & (tmp_df_descendants_muts$has_mut))
    n_with_mut_without_descendants <- sum((!tmp_df_descendants_muts$has_descendants) & (tmp_df_descendants_muts$has_mut))
    n_without_mut_with_descendants <- sum((tmp_df_descendants_muts$has_descendants) & (!tmp_df_descendants_muts$has_mut))
    n_without_mut_without_descendants <- sum((!tmp_df_descendants_muts$has_descendants) & (!tmp_df_descendants_muts$has_mut))
    
    mat_fisher <- matrix(c(n_with_mut_with_descendants, n_with_mut_without_descendants, n_without_mut_with_descendants, n_without_mut_without_descendants),
                         nrow = 2, ncol = 2)
    f_test <- fisher.test(x = mat_fisher)
    
    tibble(gene = df_muts_to_include$gene[i_row],
           aa_muts_direction = df_muts_to_include$aa_muts_direction[i_row],
           is_synonymous = df_muts_to_include$is_synonymous[i_row],
           OR = as.numeric(f_test$estimate),
           OR_lower = as.numeric(f_test$conf.int[1]),
           OR_upper = as.numeric(f_test$conf.int[2]),
           p_value = as.numeric(f_test$p.value),
           n_with_mut_with_descendants = n_with_mut_with_descendants,
           n_with_mut_without_descendants = n_with_mut_without_descendants,
           n_without_mut_with_descendants = n_without_mut_with_descendants,
           n_without_mut_without_descendants = n_without_mut_without_descendants)
  }))
  
  return(df_OR_having_descendants_by_mut)
}


get_OR_having_descendants_having_muts_among_enriched_muts_non_syn <- function(df_aa_muts, 
                                                                              df_descendants,
                                                                              enriched_AA_mut){
  
  df_OR_having_descendants_by_mut <- Reduce('bind_rows', lapply(1:nrow(enriched_AA_mut), FUN = function(i_row){
    tmp_df_has_mut <- df_aa_muts %>% 
      filter(gene == enriched_AA_mut$gene[i_row],
             aa_muts_direction == enriched_AA_mut$mutation[i_row]) %>% 
      mutate(has_mut = T) %>% 
      select(strain, has_mut)
    
    tmp_df_descendants_muts <- df_descendants %>% 
      left_join(tmp_df_has_mut, by = 'strain') %>% 
      mutate(has_mut = replace_na(has_mut, F))
    
    n_with_mut_with_descendants <- sum((tmp_df_descendants_muts$has_descendants) & (tmp_df_descendants_muts$has_mut))
    n_with_mut_without_descendants <- sum((!tmp_df_descendants_muts$has_descendants) & (tmp_df_descendants_muts$has_mut))
    n_without_mut_with_descendants <- sum((tmp_df_descendants_muts$has_descendants) & (!tmp_df_descendants_muts$has_mut))
    n_without_mut_without_descendants <- sum((!tmp_df_descendants_muts$has_descendants) & (!tmp_df_descendants_muts$has_mut))
    
    mat_fisher <- matrix(c(n_with_mut_with_descendants, n_with_mut_without_descendants, n_without_mut_with_descendants, n_without_mut_without_descendants),
                         nrow = 2, ncol = 2)
    f_test <- fisher.test(x = mat_fisher)
    
    tibble(gene = enriched_AA_mut$gene[i_row],
           aa_muts_direction = enriched_AA_mut$mutation[i_row],
           OR = as.numeric(f_test$estimate),
           OR_lower = as.numeric(f_test$conf.int[1]),
           OR_upper = as.numeric(f_test$conf.int[2]),
           p_value = as.numeric(f_test$p.value),
           n_with_mut_with_descendants = n_with_mut_with_descendants,
           n_with_mut_without_descendants = n_with_mut_without_descendants,
           n_without_mut_with_descendants = n_without_mut_with_descendants,
           n_without_mut_without_descendants = n_without_mut_without_descendants)
  }))
  
  return(df_OR_having_descendants_by_mut)
}
get_OR_having_descendants_having_muts_at_pos_among_enriched_pos_non_syn <- function(df_aa_muts, 
                                                                                    df_descendants,
                                                                                    enriched_AA_pos){
  
  df_OR_having_descendants_by_mut <- Reduce('bind_rows', lapply(1:nrow(enriched_AA_pos), FUN = function(i_row){
    tmp_df_has_mut <- df_aa_muts %>% 
      filter(gene == enriched_AA_pos$gene[i_row],
             aa_muts_position == enriched_AA_pos$position[i_row]) %>% 
      mutate(has_mut = T) %>% 
      select(strain, has_mut)
    
    tmp_df_descendants_muts <- df_descendants %>% 
      left_join(tmp_df_has_mut, by = 'strain') %>% 
      mutate(has_mut = replace_na(has_mut, F))
    
    n_with_mut_with_descendants <- sum((tmp_df_descendants_muts$has_descendants) & (tmp_df_descendants_muts$has_mut))
    n_with_mut_without_descendants <- sum((!tmp_df_descendants_muts$has_descendants) & (tmp_df_descendants_muts$has_mut))
    n_without_mut_with_descendants <- sum((tmp_df_descendants_muts$has_descendants) & (!tmp_df_descendants_muts$has_mut))
    n_without_mut_without_descendants <- sum((!tmp_df_descendants_muts$has_descendants) & (!tmp_df_descendants_muts$has_mut))
    
    mat_fisher <- matrix(c(n_with_mut_with_descendants, n_with_mut_without_descendants, n_without_mut_with_descendants, n_without_mut_without_descendants),
                         nrow = 2, ncol = 2)
    f_test <- fisher.test(x = mat_fisher)
    
    tibble(gene = enriched_AA_pos$gene[i_row],
           position = enriched_AA_pos$position[i_row],
           OR = as.numeric(f_test$estimate),
           OR_lower = as.numeric(f_test$conf.int[1]),
           OR_upper = as.numeric(f_test$conf.int[2]),
           p_value = as.numeric(f_test$p.value),
           n_with_mut_with_descendants = n_with_mut_with_descendants,
           n_with_mut_without_descendants = n_with_mut_without_descendants,
           n_without_mut_with_descendants = n_without_mut_with_descendants,
           n_without_mut_without_descendants = n_without_mut_without_descendants)
  }))
  
  return(df_OR_having_descendants_by_mut)
}

df_prop_with_descendants_by_mut <- 
  get_OR_having_descendants_having_muts_among_enriched_muts_non_syn(df_aa_muts, 
                                                                    df_descendants,
                                                                    enriched_AA_mut) %>% 
  mutate(n_with_mut = n_with_mut_with_descendants + n_with_mut_without_descendants,
         prop_with_descendants = n_with_mut_with_descendants/n_with_mut) %>% 
  mutate(full_mut = paste0(gene, ":", aa_muts_direction))

df_prop_with_descendants_by_mut$p_value_BH_descendants <- p.adjust(df_prop_with_descendants_by_mut$p_value, method = 'BH')

scatter_plot_prop_with_desc_func_nb_branches_with_muts <- df_prop_with_descendants_by_mut %>% rename(OR_desc = OR) %>% 
  left_join(df_OR_AA_mut %>% select(gene, mutation, OR) %>% rename(OR_mut = OR),
            by = c('gene', 'aa_muts_direction' = 'mutation')) %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  ggplot(aes(x = n_with_mut, 
             y = prop_with_descendants, colour = gene)) +
  geom_hline(yintercept = prop_salt_with_descendants, color = 'darkgrey') +
  annotate("rect", 
           xmin = 0, xmax = Inf, 
           ymin = CI_prop_salt_with_descendants[1],
           ymax = CI_prop_salt_with_descendants[2], 
           alpha = 0.2, fill = 'darkgrey') +
  geom_point() +
  ggrepel::geom_text_repel(data = df_prop_with_descendants_by_mut %>%
                             rename(OR_desc = OR) %>% 
                             left_join(df_OR_AA_mut %>% select(gene, mutation, OR) %>% rename(OR_mut = OR),
                                       by = c('gene', 'aa_muts_direction' = 'mutation')) %>% 
                             filter(OR_mut > 5, n_with_mut > 4),
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
  theme(panel.grid.minor.x = element_blank()) +
  guides(colour = guide_legend(override.aes = list(label = '')))

set.seed(2)
pdf('../figure_descendant_analysis/prop_branches_with_desc_func_nb_branches_with_mut.pdf',
    height = 5.5, width = 8.5)
plot(scatter_plot_prop_with_desc_func_nb_branches_with_muts)
dev.off()  
png('../figure_descendant_analysis/prop_branches_with_desc_func_nb_branches_with_mut.png',
    height = 5.5, width = 8.5, res = 350, units = 'in')
plot(scatter_plot_prop_with_desc_func_nb_branches_with_muts)
dev.off()  


##
df_aa_reference <- read_tsv('../../results/mutation_enrichment/df_aa_reference_sequence.tsv') %>% 
  mutate(full_aa = paste0(gene, ':', AA_pos, AA))
recombinants_clades_to_drop <- c('24F (XEC)', '25G (XFC)', '22F (XBB)', '25C (XFG)')
df_mut_counts_by_nextstrain_clade <- read_tsv('../../results/mutation_enrichment/df_aa_freq_combined_CV_1.8_BF_20.tsv') %>% 
  mutate(prop_AA = count / clade_size, full_mut = paste0(gene, ':', pos, AA)) %>% 
  filter(AA != '-', ! clade %in% recombinants_clades_to_drop, is_lineage_defining | (is_within_clade_mut & (prop_AA > 0.8))) %>% 
  group_by(full_mut) %>% summarise(n_clades = n())


bootstrap_mean_function <-  function(x, B = 2000) {
  boot_means <- replicate(B, mean(sample(x, replace = TRUE), na.rm = TRUE))
  
  tibble(
    mean = mean(x, na.rm = TRUE),
    lwr = quantile(boot_means, 0.025, na.rm = TRUE),
    upr = quantile(boot_means, 0.975, na.rm = TRUE)
  )
}

  
df_for_plot <- df_prop_with_descendants_by_mut %>% 
  left_join(df_mut_counts_by_nextstrain_clade %>% 
              filter(full_mut %in% df_prop_with_descendants_by_mut$full_mut)) %>% 
  mutate(n_clades = replace_na(n_clades, 0)) %>% 
  mutate(quantile_prop = ntile(prop_with_descendants, 5)) 
  

df_for_plot %>% 
  group_by(quantile_prop) %>% 
  reframe(
    mean_number_clades = mean(n_clades),
    mean_prop_within_quantile = mean(prop_with_descendants, na.rm = TRUE),
    bootstrap_mean_function(n_clades)
  ) %>% 
  ggplot(aes(x = quantile_prop, y = mean_number_clades)) +
  geom_point() +
  geom_linerange(aes(ymin = lwr, ymax = upr)) +
  scale_x_continuous(name = 'Quintile of proportion of saltational\nbranches with descendants') +
  scale_y_continuous(name = 'Mean number of clades\nin which mutation fixed',
                     breaks = seq(0, 2., 0.2)) + 
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        panel.grid.minor = element_blank())
    


## Reversions
df_prop_with_descendants_by_mut %>% 
  mutate(is_reversion = full_mut %in% df_aa_reference$full_aa) %>% 
  group_by(is_reversion) %>% summarise(n_with_mut_with_descendants = sum(n_with_mut_with_descendants), n_with_mut = sum(n_with_mut))
  ggplot(aes(x = is_reversion, 
             y = n_with_mut_with_descendants / n_with_mut)) +
  geom_point()
  geom_hline(yintercept = prop_salt_with_descendants, color = 'darkgrey') +
  annotate("rect", 
           xmin = 0, xmax = Inf, 
           ymin = CI_prop_salt_with_descendants[1],
           ymax = CI_prop_salt_with_descendants[2], 
           alpha = 0.2, fill = 'darkgrey') +
  geom_point() 
  

### Volcano plot
df_prop_with_descendants_by_mut %>% rename(OR_desc = OR) %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  ggplot(aes(x = log10(OR_desc), y = - log10(p_value), colour = gene)) +
  geom_point() +
  scale_colour_manual(values = palette_genes,
                      name = 'Gene') +
  theme_bw() 


df_prop_with_descendants_by_mut %>% 
  filter(n_with_mut_with_descendants + n_with_mut_without_descendants >= 4) %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  ggplot(aes(x = reorder(paste0(gene , ':', aa_muts_direction), OR), y = OR, colour = gene)) +
  geom_point() +
  geom_linerange(aes(ymin = OR_lower, ymax = OR_upper)) +
  facet_wrap(. ~ gene == 'S', scales = 'free') +
  scale_y_continuous(trans = 'log') +
  scale_colour_manual(values = palette_genes, 
                      name = 'Gene') +
  coord_flip()

