library(tidyverse)
library(scales)
library(glmnet)
source('utils_descendants_analysis.R')
source('config_plot.R')

# Load saltational branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.7_BF_20.tsv', header = F) %>%
  unlist() %>% as.character()

# Load descendants information
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(has_descendants = n_terminal_descendants > 1) %>% 
  select(strain, has_descendants)

df_descendant_counts <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(has_descendants = n_terminal_descendants > 1)


# Load mutation enrichment information
df_OR_AA_pos <- read_tsv('../../results/mutation_enrichment/df_AA_pos_enrichment_saltational_CV_1.7_BF_20.tsv') %>% 
  filter(threshold_for_inclusion == 1)
df_OR_AA_mut <- read_tsv('../../results/mutation_enrichment/df_AA_mut_enrichment_saltational_CV_1.7_BF_20.tsv') %>% 
  filter(threshold_for_inclusion == 1)

enriched_AA_mut <- df_OR_AA_mut %>% filter(p_value_BH < 0.05) %>% select(gene, mutation) %>% 
  mutate(mut_full = paste0(gene, ':', mutation))
enriched_AA_pos <- df_OR_AA_pos %>% filter(p_value_BH < 0.05) %>% select(gene, position) %>% 
  mutate(mut_full = paste0(gene, ':', position))

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

df_OR_having_descendants_by_mut <- 
  get_OR_having_descendants_having_muts_among_enriched_muts_non_syn(df_aa_muts, 
                                                                    df_descendants,
                                                                    enriched_AA_mut)
df_OR_having_descendants_by_mut$p_value_BH <- p.adjust(df_OR_having_descendants_by_mut$p_value, method = 'BH')

df_OR_having_descendants_by_pos <- 
  get_OR_having_descendants_having_muts_at_pos_among_enriched_pos_non_syn(df_aa_muts, 
                                                                          df_descendants,
                                                                          enriched_AA_pos)
df_OR_having_descendants_by_pos$p_value_BH <- p.adjust(df_OR_having_descendants_by_pos$p_value, method = 'BH')



# df_OR_having_descendants_by_mut <- 
#   get_OR_having_descendants_having_muts(df_aa_muts, 
#                                         df_descendants,
#                                         threshold_mutation_occurence_inclusion = 0)



df_OR_having_descendants_by_mut %>% filter(OR_lower > 1. | OR_upper < 1.) %>% View()

df_OR_having_descendants_by_pos %>% filter(OR_lower > 1.| OR_upper < 1.) %>% View()

df_aa_muts %>% filter(gene == 'ORF1b', aa_muts_direction == '1000L') %>% 
  left_join(df_descendant_counts, by = 'strain') %>% View()
df_aa_muts %>% filter(gene == 'M', aa_muts_direction == '82T') %>% 
  left_join(df_descendant_counts, by = 'strain') %>% View()
df_aa_muts %>% filter(gene == 'S', aa_muts_direction == '477N') %>% 
  left_join(df_descendant_counts, by = 'strain') %>% View()

df_OR_having_descendants_by_pos %>% 
  arrange(-OR) %>% mutate(rank = 1:n()) %>% 
  ggplot(aes(x = rank, y = OR, colour = gene, alpha = p_value < 0.05)) +
  geom_point() +
  geom_linerange(aes(ymin = OR_lower, ymax = OR_upper)) +
  scale_y_continuous(trans = 'log') +
  scale_colour_manual(values = palette_genes, name = 'Gene') +
  scale_alpha_manual(breaks = c(T, F), values = c(1., 0.5), guide = 'none')

df_OR_having_descendants_by_mut %>% 
  filter(n_with_mut_with_descendants + n_with_mut_without_descendants >= 5 | OR_lower > 1. | OR_upper < 1.) %>% 
  arrange(- OR ) %>% mutate(rank = 1:n()) %>% 
  ggplot(aes(x = rank, y = OR, colour = gene)) +
  geom_point() +
  geom_linerange(aes(ymin = OR_lower, ymax = OR_upper)) +
  scale_y_continuous(trans = 'log') +
  scale_colour_manual(values = palette_genes)
  

df_OR_having_descendants_by_mut %>% 
  filter(n_with_mut_with_descendants + n_with_mut_without_descendants >= 5 | OR_lower > 1. | OR_upper < 1.) %>% 
  arrange(- OR ) %>% mutate(rank = 1:n()) %>% 
  ggplot(aes(x = rank, 
             y = n_with_mut_with_descendants/(n_with_mut_with_descendants + n_with_mut_without_descendants), 
             fill = gene, colour = gene)) +
  geom_bar(stat = 'identity') +
  #geom_linerange(aes(ymin = OR_lower, ymax = OR_upper)) +
  #scale_y_continuous(trans = 'log') +
  scale_colour_manual(values = palette_genes) +
  scale_fill_manual(values = palette_genes)

df_OR_having_descendants_by_mut %>% 
  filter(OR_lower > 1.)
