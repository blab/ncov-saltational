library(tidyverse)

# Load saltational branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', header = F) %>%
  unlist() %>% as.character()

# Load reversion information
df_reversion <- read_tsv('../../results/mutation_enrichment/df_aa_reference_sequence.tsv') %>% 
  mutate(reversion = paste0(gene, ':', AA_pos, AA))

# Load AA mutation information
df_OR_AA_muts <- read_tsv('../../results/mutation_enrichment/df_AA_mut_enrichment_saltational_opportunity_CV_1.8_BF_20.tsv') %>% 
  filter(threshold_for_inclusion == 2) %>% 
  mutate(full_mut = paste0(gene, ':', mutation)) %>% 
  mutate(is_reversion = full_mut %in% df_reversion$reversion)

name_enriched_muts <- df_OR_AA_muts %>% filter(p_value_BH < 0.05, OR > 1.) %>% select(full_mut) %>% 
  unlist() %>% as.character()

comp_OR_mut_reversion <- df_OR_AA_muts %>% 
  ggplot(aes(x = is_reversion, y = log((n_mut_of_interest_adapt_branches + 0) * (n_other_muts_normal_branches + 0) / (n_other_muts_adapt_branches + 0) / (n_mut_of_interest_normal_branches + 0)))) +
  geom_jitter(alpha = 0.1, width = 0.2, aes(colour = is_reversion)) +
  geom_boxplot(outliers = F, aes(colour = is_reversion), fill = NA) +
  ggsignif::geom_signif(comparisons = list(c('FALSE', 'TRUE')),
                        map_signif_level = TRUE, 
                        tip_length = 0.005 ) +
  scale_x_discrete(name = 'Mutation is reversion', breaks = c(F, T), labels = c('No', 'Yes'),
                   limits = rev) +
  scale_y_continuous(name = expression(OR['mutation']^{'saltational'}), 
                     trans = 'log',
                     breaks = c(0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5., 10.)) +
  scale_colour_manual(values = c('darkgrey', 'deeppink2'), guide = 'none') +
  scale_fill_manual(values = c('darkgrey', 'deeppink2'), guide = 'none') +
  theme_classic() +
  coord_flip(ylim = c(0.05, NA)) +
  theme(axis.text.x = element_text(size = 11))




# Load descendants information
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(has_descendants = n_terminal_descendants > 1) %>% 
  select(strain, has_descendants)

df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  filter(is_coded) %>% 
  rename(strain = node_id) %>% 
  mutate(is_saltational = strain %in% name_saltational_branches) %>% 
  filter(is_saltational) %>% 
  mutate(is_synonymous = (aa_initial == aa_new),
         aa_muts_direction = paste0(aa_muts_position, aa_new),
         full_mut = paste0(gene, ':', aa_muts_direction)) %>% 
  mutate(is_reversion = full_mut %in% df_reversion$reversion) %>% 
  filter(! is_synonymous) %>% 
  left_join(df_descendants) %>% 
  group_by(full_mut, gene, is_reversion, has_descendants) %>% 
  summarise(n_branches = n()) %>% 
  group_by(full_mut, gene, is_reversion) %>% 
  summarise(n_with_descendants = sum(n_branches[has_descendants]),
            n_with_mut = sum(n_branches))


comp_prop_salt_branches_with_descendants_reversion <- df_aa_muts %>% 
  filter(full_mut %in% name_enriched_muts) %>% 
  group_by(is_reversion) %>% 
  ggplot(aes(x = is_reversion, y = n_with_descendants / n_with_mut)) +
  geom_jitter(alpha = 0.1, width = 0.2, aes(colour = is_reversion)) +
  geom_boxplot(outliers = F, aes(colour = is_reversion), fill = NA) +
  
  ggsignif::geom_signif(comparisons = list(c('FALSE', 'TRUE')),
                        map_signif_level = TRUE, 
                        tip_length = 0.005 ) +
  scale_x_discrete(name = 'Mutation is reversion', breaks = c(F, T), labels = c('No', 'Yes'),
                   limits = rev) +
  scale_y_continuous(name = 'Proportion of saltational branches with descendants') +
  scale_fill_manual(values = c('darkgrey', 'deeppink2'), guide = 'none') +
  scale_colour_manual(values = c('darkgrey', 'deeppink2'), guide = 'none') +
  theme_classic() +
  coord_flip() +
  theme(axis.text.x = element_text(size = 11))


ggpubr::ggarrange(comp_OR_mut_reversion, 
                  comp_prop_salt_branches_with_descendants_reversion, 
                  nrow = 2, labels = 'AUTO')

