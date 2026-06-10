library(tidyverse)

# Define AA windows to look at
char_gene_windows <- tibble(gene = c('E',  'E', 'M', 'M', 'M','M','ORF1a', 'ORF1a', 'ORF1a', 'ORF1a', 'ORF1a'),
                            gene_start = c(5, 14, 2, 2, 3,  4, 1322, 1361, 4085, 4097, 4164),
                            gene_end = c(11, 16, 4, 2, 3, 4, 1324, 1372, 4090, 4102, 4167))

# Load name of saltational branches
name_saltational_branches <- 
  read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
             header = F) %>% unlist() %>% as.character()

# Load AA muts information
df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  rename(strain = node_id) %>% 
  filter(is_coded, aa_initial != aa_new) %>% 
  mutate(is_saltational = strain %in% name_saltational_branches)



Reduce('bind_rows', lapply(1:nrow(char_gene_windows), FUN = function(i_row){
  gene_window <- char_gene_windows$gene[i_row]
  gene_start <- char_gene_windows$gene_start[i_row]
  gene_end <- char_gene_windows$gene_end[i_row]
  
  
  tmp_df_muts <- df_aa_muts %>% 
    mutate(is_within_gene = ((gene == gene_window) & (aa_muts_position >= gene_start)) & (aa_muts_position <= gene_end)) %>% 
    group_by(is_within_gene, is_saltational) %>% 
    summarise(n_nonsyn = n())
  
  mat_fisher <- matrix(c(
    tmp_df_muts$n_nonsyn[tmp_df_muts$is_within_gene & tmp_df_muts$is_saltational],
    tmp_df_muts$n_nonsyn[tmp_df_muts$is_within_gene & !tmp_df_muts$is_saltational],
    tmp_df_muts$n_nonsyn[!tmp_df_muts$is_within_gene & tmp_df_muts$is_saltational],
    tmp_df_muts$n_nonsyn[!tmp_df_muts$is_within_gene & !tmp_df_muts$is_saltational]
  ), ncol = 2, nrow = 2)
  
  f_test <- fisher.test(x = mat_fisher)
  tibble(gene = gene_window,
         gene_start = gene_start,
         gene_end = gene_end,
         OR = f_test$estimate[1],
         OR_lower = f_test$conf.int[1],
         OR_upper = f_test$conf.int[2])
}))
