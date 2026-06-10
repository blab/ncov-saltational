library(tidyverse)

# Load name of saltational branches
name_saltational_branches <- 
  read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
             header = F) %>% unlist() %>% as.character()

# Load nucleotide mutation information
df_nuc_mut <- read_tsv('../../results/mutations/df_nuc_mut.tsv')

# Get nucleotide mutation type by gene
df_nuc_mut_type_distrib <- df_nuc_mut %>% 
  group_by(strain, nuc_mutation) %>% 
  mutate(nuc_old = substr(nuc_mutation, 1, 1), 
         nuc_new = substr(nuc_mutation, nchar(nuc_mutation), nchar(nuc_mutation)),
         is_saltational = strain %in% name_saltational_branches,
         mut_type = paste0(nuc_old, nuc_new)) %>% 
  group_by(is_saltational, mut_type) %>% 
  summarise(n_occ = n()) %>% 
  group_by(is_saltational) %>% mutate(prop_occ = n_occ / sum(n_occ))

plt_nuc_mut_type <- df_nuc_mut_type_distrib %>% 
  ggplot(aes(x = mut_type, y = prop_occ, fill = is_saltational)) +
  geom_bar(stat = 'identity', position = position_dodge(0.8), width = 0.8) +
  scale_x_discrete(name = 'Nucleotide mutation type') +
  scale_y_continuous(name = 'Proportion') +
  scale_fill_manual(values = c('orange2', 'darkgrey'), breaks = c(T, F), labels = c('Saltational', 'Others'),
                    name = 'Branch type') +
  theme_classic() +
  theme(legend.position = 'inside', 
        legend.position.inside = c(0.8, 0.8),
        legend.background = element_blank(),
        axis.text = element_text(size = 11))

png('../../figures/figure_mutation_enrichment/plt_nuc_mut_type.png', 
    height = 2.8, width = 4.5, res = 350, units = 'in')
plot(plt_nuc_mut_type)
dev.off()

##################
df_gene_coordinates <- read_csv('../param/proteins_nucleotide_coordinates.csv') %>% 
  filter(! startsWith(gene, 'NSP'), gene != 'RdRp', gene != 'S1', gene != 'S2', gene != 'ORF1ab', gene != 'ORF9b')

df_aa_mut <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  rename(strain = node_id) %>% 
  mutate(aa_muts_position = as.integer(aa_muts_position))


get_gene <- function(vec_positions){
  sapply(vec_positions, FUN = function(position){
    gene_to_return <- df_gene_coordinates %>% 
      filter(start <= position, stop >= position) %>% 
      select(gene) %>% unlist() %>% as.character()
    if(length(gene_to_return) == 0){
      gene_to_return <- NA
    }
    gene_to_return[1]
  }) %>% unlist() %>% return()
}

tmp_df <- df_nuc_mut %>% 
  group_by(strain, nuc_mutation) %>% 
  mutate(nuc_old = substr(nuc_mutation, 1, 1), 
         nuc_new = substr(nuc_mutation, nchar(nuc_mutation), nchar(nuc_mutation)),
         position = as.integer(substr(nuc_mutation, 2, nchar(nuc_mutation) - 1)),
         is_saltational = strain %in% name_saltational_branches,
         mut_type = paste0(nuc_old, nuc_new)) %>% 
  left_join(tibble(position = c(1:30000), gene = get_gene(c(1:30000))), by = 'position') %>% 
  left_join(df_gene_coordinates) %>% 
  mutate(aa_position = as.integer(floor((position - start)/3)) + 1) %>% 
  left_join(df_aa_mut, by = c('strain', 'gene', 'aa_position' = 'aa_muts_position'))

df_nuc_mut_type_distrib_split_syn_nonsyn <- tmp_df %>% 
  mutate(is_synonymous = case_when(is.na(aa_position) ~ T,
                                   ! is_coded ~ T,
                                   aa_initial == aa_new ~ T,
                                   aa_initial != aa_new ~ F)) %>% 
  group_by(is_saltational, is_synonymous, mut_type) %>% 
  summarise(n_occ = n()) %>% 
  group_by(is_saltational, is_synonymous) %>% mutate(prop_occ = n_occ / sum(n_occ))

df_nuc_mut_type_distrib_split_syn_nonsyn_by_gene <- tmp_df %>% 
  mutate(is_synonymous = case_when(is.na(aa_position) ~ T,
                                   ! is_coded ~ T,
                                   aa_initial == aa_new ~ T,
                                   aa_initial != aa_new ~ F)) %>% 
  group_by(is_saltational, is_synonymous, mut_type, gene) %>% 
  summarise(n_occ = n()) %>% 
  group_by(is_saltational, is_synonymous, gene) %>% mutate(prop_occ = n_occ / sum(n_occ))


plt_nuc_mut_type_split_syn_nonsyn <- df_nuc_mut_type_distrib_split_syn_nonsyn %>% 
  mutate(aa_mut_type = ifelse(is_synonymous, 'Synonymous', 'Nonsynonymous')) %>% 
  filter(! is.na(aa_mut_type)) %>% 
  ggplot(aes(x = mut_type, y = prop_occ, fill = is_saltational)) +
  geom_bar(stat = 'identity', position = position_dodge(0.8), width = 0.8) +
  scale_x_discrete(name = 'Nucleotide mutation type') +
  scale_y_continuous(name = 'Proportion') +
  scale_fill_manual(values = c('orange2', 'darkgrey'), breaks = c(T, F), labels = c('Saltational', 'Others'),
                    name = 'Branch type') +
  theme_classic() +
  facet_wrap(. ~ aa_mut_type) +
  theme(legend.position = 'inside', 
        legend.position.inside = c(0.9, 0.8),
        legend.background = element_blank(),
        axis.text = element_text(size = 11))

plt_nuc_mut_type_split_syn_nonsyn_by_gene <- df_nuc_mut_type_distrib_split_syn_nonsyn_by_gene %>% 
  mutate(aa_mut_type = ifelse(is_synonymous, 'Synonymous', 'Nonsynonymous')) %>% 
  filter(! is.na(aa_mut_type)) %>% 
  mutate(gene = case_when(gene %in% c('ORF1a', 'ORF1b', 'S', 'N') ~ gene,
                          is.na(gene) ~ 'Non coding',
                          T ~ 'Other genes')) %>% 
  mutate(type = paste0(gene, ' - ', aa_mut_type)) %>% 
  ggplot(aes(x = mut_type, y = prop_occ, fill = is_saltational)) +
  geom_bar(stat = 'identity', position = position_dodge(0.8), width = 0.8) +
  scale_x_discrete(name = 'Nucleotide mutation type') +
  scale_y_continuous(name = 'Proportion') +
  scale_fill_manual(values = c('orange2', 'darkgrey'), breaks = c(T, F), labels = c('Saltational', 'Others'),
                    name = 'Branch type') +
  theme_classic() +
  facet_wrap(. ~ type) +
  theme(legend.position = 'inside', 
        legend.position.inside = c(0.9, 0.1),
        legend.background = element_blank(),
        axis.text = element_text(size = 11))

png('../../figures/figure_mutation_enrichment/plt_nuc_mut_type_split_syn_nonsyn.png', 
    height = 2.8, width = 9, res = 350, units = 'in')
plot(plt_nuc_mut_type_split_syn_nonsyn)
dev.off()

png('../../figures/figure_mutation_enrichment/plt_nuc_mut_type_split_syn_nonsyn_by_gene.png', 
    height = 5., width = 12, res = 350, units = 'in')
plot(plt_nuc_mut_type_split_syn_nonsyn_by_gene)
dev.off()

