library(tidyverse)
library(ggrepel)
library(ggpubr)
source('config_plot.R')

# Load OR by AA mutation
df_OR_AA_mut <- read_tsv('../../results/mutation_enrichment/df_AA_mut_enrichment_saltational_CV_1.8_BF_20.tsv') %>% 
  filter(threshold_for_inclusion == 1)
df_OR_AA_mut

# Load Ryan's EPCI mutations
df_chronic_muts <- read_csv('../../data/frequent_mutations/ryan_epci_hqcs_AA_ratio.csv') %>% 
  tidyr::separate(Mutation, into = c("gene", "aa_muts"), sep = ":") %>% 
  group_by(gene, aa_muts) %>% 
  mutate(mutation = substr(aa_muts, start = 2, stop = nchar(aa_muts)),
         aa_initial = substr(aa_muts, start = 1, stop = 1),
         aa_new = substr(aa_muts, start = nchar(aa_muts), stop = nchar(aa_muts))) %>% 
  ungroup() %>% 
  mutate(hqcs_count_estim = case_when(is.infinite(EPCI_to_HQCS_ratio) ~ 0,
                                      T ~ round(EPCI_ct / EPCI_to_HQCS_ratio))) %>% 
  filter(aa_new != '-', aa_initial != '-') %>% 
  group_by(gene, mutation) %>% 
  summarise(n_mut_chronic = sum(EPCI_ct),
            n_mut_high_quality_sequence = sum(hqcs_count_estim))

# Compute modified OR of mutations of occurring in chronic infections
df_OR_chronic <- Reduce('bind_rows', lapply(1:nrow(df_chronic_muts), FUN = function(i_row){
  curr_gene <- df_chronic_muts$gene[i_row]
  curr_mut <- df_chronic_muts$mutation[i_row]
  
  tmp_df <- df_chronic_muts %>% 
    mutate(is_mut = (mutation == curr_mut) & (gene == curr_gene)) %>% 
    group_by(is_mut) %>% 
    summarise(n_chronic = sum(n_mut_chronic),
              n_high_quality_sequences = sum(n_mut_high_quality_sequence))
  
  n_mut_chronic <- tmp_df$n_chronic[tmp_df$is_mut]
  n_mut_background <- tmp_df$n_high_quality_sequences[tmp_df$is_mut]
  
  n_other_muts_chronic <- tmp_df$n_chronic[!tmp_df$is_mut]
  n_other_muts_background <- tmp_df$n_high_quality_sequences[!tmp_df$is_mut]
  
  mat_fisher <- matrix(c(n_mut_chronic, n_other_muts_chronic, n_mut_background, n_other_muts_background),
                       ncol = 2, nrow = 2)
  f_test <- fisher.test(x = mat_fisher)
  
  mat_fisher_modif <- matrix(c(n_mut_chronic + 1, n_other_muts_chronic + 1,
                               n_mut_background + 1, n_other_muts_background + 1),
                             ncol = 2, nrow = 2)
  f_test_modif <- fisher.test(x = mat_fisher_modif)
  
  tibble(gene = curr_gene, 
         mutation = curr_mut,
         OR_chronic = f_test$estimate,
         OR_chronic_lower = f_test$conf.int[1],
         OR_chronic_upper = f_test$conf.int[2],
         OR_chronic_modif = f_test_modif$estimate,
         OR_chronic_modif_lower = f_test_modif$conf.int[1],
         OR_chronic_modif_upper = f_test_modif$conf.int[2],
         p_value_chronic = f_test$p.value)
}))

# Compute modified OR of mutations of occurring in chronic infections
vec_name_node_adapt <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv') %>% unlist() %>% as.character()

df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  filter(is_coded) %>% # Coded AA mut
  filter(aa_new != aa_initial) %>% # Non-synonymous
  mutate(full_mut = paste0(aa_muts_position, aa_new),
         is_saltational = node_id %in% vec_name_node_adapt) %>% 
  group_by(gene, full_mut, is_saltational) %>% 
  summarise(n_muts = n())
  
df_modif_OR_saltational <- Reduce('bind_rows', lapply(1:nrow(df_chronic_muts), FUN = function(i_row){
  curr_gene <- df_chronic_muts$gene[i_row]
  curr_mut <- df_chronic_muts$mutation[i_row]
  
  tmp_df <- df_aa_muts %>% 
    mutate(is_mut = (full_mut == curr_mut) & (gene == curr_gene)) %>% 
    group_by(is_mut) %>% 
    summarise(n_saltational = sum(n_muts[is_saltational]),
              n_others = sum(n_muts[! is_saltational]))
  
  n_mut_saltational <- tmp_df$n_saltational[tmp_df$is_mut]
  n_mut_background <- tmp_df$n_others[tmp_df$is_mut]
  
  n_other_muts_saltational <- tmp_df$n_saltational[!tmp_df$is_mut]
  n_other_muts_background <- tmp_df$n_others[!tmp_df$is_mut]
  
  mat_fisher <- matrix(c(n_mut_saltational, n_other_muts_saltational, n_mut_background, n_other_muts_background),
                       ncol = 2, nrow = 2)
  f_test <- fisher.test(x = mat_fisher)
  
  mat_fisher_modif <- matrix(c(n_mut_saltational + 1, n_other_muts_saltational + 1, 
                               n_mut_background + 1, n_other_muts_background + 1),
                             ncol = 2, nrow = 2)
  f_test_modif <- fisher.test(x = mat_fisher_modif)
  
  tibble(gene = curr_gene, 
         mutation = curr_mut,
         OR_saltational = f_test$estimate,
         OR_saltational_lower = f_test$conf.int[1],
         OR_saltational_upper = f_test$conf.int[2],
         OR_saltational_modif = f_test_modif$estimate,
         OR_saltational_modif_lower = f_test_modif$conf.int[1],
         OR_saltational_modif_upper = f_test_modif$conf.int[2],
         p_value_salt = f_test$p.value)
}))

## 
df_for_plot <- df_chronic_muts %>% 
  left_join(df_OR_chronic, by = c('gene', 'mutation')) %>% 
  left_join(df_modif_OR_saltational, by = c('gene', 'mutation'))

cor(log(df_for_plot$OR_chronic_modif), log(df_for_plot$OR_saltational_modif), method = 'spearman')
cor(log(df_for_plot$OR_chronic_modif), log(df_for_plot$OR_saltational_modif), method = 'pearson')

label_125 <- function(x) {
  m <- signif(x / 10^floor(log10(x)), 12) 
  ifelse(m %in% c(1, 2, 5), x, "")
}

plt_comp_OR_saltational_OR_chronic <- df_for_plot %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  arrange(n_mut_high_quality_sequence) %>% 
  ggplot(aes(x = OR_chronic_modif, y = OR_saltational_modif)) +
  geom_point(aes(colour = gene), alpha = 0.8) +
  stat_cor(method = 'pearson') +
  scale_x_continuous(trans = 'log', name = expression(OR['chronic']),
                     labels = label_125,
                     breaks = c(seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 9e-1, 1e-1),
                                seq(1e0, 9e0, 1e0), seq(1e1, 9e1, 1e1), seq(1e2, 9e2, 1e2))) +
  scale_y_continuous(trans = 'log', name = expression(OR['saltational']),
                     labels = label_125,
                     breaks = c(seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 9e-1, 1e-1),
                                seq(1e0, 9e0, 1e0), seq(1e1, 9e1, 1e1), seq(1e2, 9e2, 1e2))) +
  scale_colour_manual(values = palette_genes, name = 'Gene') +
  theme_bw() +
  coord_fixed() +
  theme(panel.grid.minor = element_blank())

png('../figure_mutation_enrichment/compare_OR_saltational_chronic.png',
    height = 4., width = 5.5, res = 350, units = 'in')
plot(plt_comp_OR_saltational_OR_chronic)
dev.off()
pdf('../figure_mutation_enrichment/compare_OR_saltational_chronic.pdf',
    height = 4., width = 5.5)
plot(plt_comp_OR_saltational_OR_chronic)
dev.off()
