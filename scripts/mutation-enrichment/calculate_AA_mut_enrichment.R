library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--df_aa_muts', type = 'character', help = 'TSV file with AA muts info')
  parser$add_argument('--name_adaptive_branches', type = 'character', help = 'Text file with name of adaptive branches')
  parser$add_argument('--output_df_enrichment', type = 'character', help = 'TSV file for enrichment OR scores')
  return(parser$parse_args())
}

args <- collect_args()

# Load input data
df_aa_muts <- read_tsv(args$df_aa_muts) %>% 
  rename(strain = node_id) %>% 
  # Keep non-synonymous only
  filter(aa_new != aa_initial) %>%
  # Keep coded ones only
  filter(is_coded) %>% 
  mutate(mut_direction = paste0(aa_muts_position, aa_new))
  
vec_name_node_signif <- read.table(args$name_adaptive_branches) %>% unlist() %>% as.character()

get_OR_AA_mut <- function(df_aa_muts, vec_name_node_signif, threshold_mut_occ_for_inclusion = 1){
  # Get number of mutation by branch type
  df_aa_muts_counts_by_branch_type <- df_aa_muts %>% 
    mutate(is_adaptive = strain %in% vec_name_node_signif) %>% 
    group_by(is_adaptive, mut_direction, gene) %>% 
    summarise(n_mutations = n())
  
  # Get mutated position within adaptive branches
  df_mutated_pos <- df_aa_muts_counts_by_branch_type %>% filter(is_adaptive == T) %>% 
    filter(n_mutations >= threshold_mut_occ_for_inclusion) %>% 
    select(gene, mut_direction)
    
  df_OR_by_mut <- Reduce('bind_rows', lapply(1:nrow(df_mutated_pos), FUN = function(i_row){
    curr_gene <- df_mutated_pos$gene[i_row]
    curr_pos <- df_mutated_pos$mut_direction[i_row]
    
    df_muts_counts <- df_aa_muts_counts_by_branch_type %>% 
      mutate(is_mut = (gene == curr_gene) & (mut_direction == curr_pos)) %>% 
      group_by(is_mut, is_adaptive) %>% summarise(n_mutations = sum(n_mutations))
      
    n_mut_of_interest_adapt_branches <- df_muts_counts$n_mutations[df_muts_counts$is_adaptive == T & df_muts_counts$is_mut == T]
    n_mut_of_interest_normal_branches <- df_muts_counts$n_mutations[df_muts_counts$is_adaptive == F & df_muts_counts$is_mut == T]
    n_other_muts_adapt_branches <- df_muts_counts$n_mutations[df_muts_counts$is_adaptive == T & df_muts_counts$is_mut == F]
    n_other_muts_normal_branches <- df_muts_counts$n_mutations[df_muts_counts$is_adaptive == F & df_muts_counts$is_mut == F]
    
    n_mut_of_interest_adapt_branches <- ifelse(length(n_mut_of_interest_adapt_branches) == 0, 0, n_mut_of_interest_adapt_branches)
    n_mut_of_interest_normal_branches <- ifelse(length(n_mut_of_interest_normal_branches) == 0, 0, n_mut_of_interest_normal_branches)
    n_other_muts_adapt_branches <- ifelse(length(n_other_muts_adapt_branches) == 0, 0, n_other_muts_adapt_branches)
    n_other_muts_normal_branches <- ifelse(length(n_other_muts_normal_branches) == 0, 0, n_other_muts_normal_branches)
    
    mat_fisher <- matrix(c(n_mut_of_interest_adapt_branches, n_mut_of_interest_normal_branches, n_other_muts_adapt_branches, n_other_muts_normal_branches),
                         nrow = 2, ncol = 2, byrow = T)
    
    f_test <- fisher.test(x = mat_fisher)
    
    df_fisher <- tibble(OR = as.numeric(f_test$estimate),
                        OR_lower = as.numeric(f_test$conf.int[1]),
                        OR_upper = as.numeric(f_test$conf.int[2]),
                        p_value = f_test$p.value,
                        n_mut_of_interest_adapt_branches = n_mut_of_interest_adapt_branches,
                        n_mut_of_interest_normal_branches = n_mut_of_interest_normal_branches,
                        n_other_muts_adapt_branches = n_other_muts_adapt_branches,
                        n_other_muts_normal_branches = n_other_muts_normal_branches,
                        gene = curr_gene,
                        mutation = curr_pos)
  }))
  
  
  df_OR_by_mut$p_value_BH <- p.adjust(df_OR_by_pos$p_value, method = 'BH')
  df_OR_by_mut$p_value_bonf <- p.adjust(df_OR_by_pos$p_value, method = 'bonferroni')
  
  return(df_OR_by_mut)
}

df_OR_by_AA_mut_threshold_1 <- get_OR_AA_mut(df_aa_muts, vec_name_node_signif, 1) %>% mutate(threshold_for_inclusion = 1)
df_OR_by_AA_mut_threshold_2 <- get_OR_AA_mut(df_aa_muts, vec_name_node_signif, 2) %>% mutate(threshold_for_inclusion = 2)
df_OR_by_AA_mut_threshold_3 <- get_OR_AA_mut(df_aa_muts, vec_name_node_signif, 3) %>% mutate(threshold_for_inclusion = 3)
df_OR_by_AA_mut_threshold_4 <- get_OR_AA_mut(df_aa_muts, vec_name_node_signif, 4) %>% mutate(threshold_for_inclusion = 4)
df_OR_by_AA_mut_threshold_5 <- get_OR_AA_mut(df_aa_muts, vec_name_node_signif, 5) %>% mutate(threshold_for_inclusion = 5)

df_OR_by_AA_mut_threshold <- bind_rows(df_OR_by_AA_mut_threshold_1, df_OR_by_AA_mut_threshold_2, df_OR_by_AA_mut_threshold_3,
                                       df_OR_by_AA_mut_threshold_4, df_OR_by_AA_mut_threshold_5)

# Save results
write_tsv(df_OR_by_AA_mut_threshold, args$output_df_enrichment)
