library(tidyverse)
library(ggrepel)

# Load frequent mutations in minks / deers
frequent_mutations_minks_deers <- read_csv('../../data/frequent_mutations/frequent_mutations_minks_deers.csv') %>% 
  filter(! is.na(aa_mut_direction)) %>% 
  mutate(aa_muts_tot = paste0(gene, ':', aa_mut_direction)) %>% 
  select(host, aa_muts_tot) %>% distinct() %>%
  filter(! is.na(aa_muts_tot)) %>% 
  rename(type = host)

# Load frequent mutations in persistent infection
frequent_mutations_chronic <- read_csv('../../data/frequent_mutations/frequent_mutations_chronic.csv') %>% 
  filter(source != 'Ghafari et al. Lancet ID. Occurred 4+ times') %>% 
  filter(! is.na(aa_mutation_full)) %>% 
  mutate(aa_muts_tot = paste0(gene, ':', aa_mut_direction)) %>% 
  select(type, aa_muts_tot) %>% distinct()

frequent_mutations_summary <- bind_rows(frequent_mutations_minks_deers, frequent_mutations_chronic) 
mut_profile_persistent <- frequent_mutations_summary %>% filter(type == 'Persistent') %>% select(aa_muts_tot) %>% unlist() %>% as.character()
mut_profile_minks <- frequent_mutations_summary %>% filter(type == 'Minks') %>% select(aa_muts_tot) %>% unlist() %>% as.character()
mut_profile_deers <- frequent_mutations_summary %>% filter(type == 'Deers') %>% select(aa_muts_tot) %>% unlist() %>% as.character()


# Load saltational branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
                                        header = F) %>% unlist() %>% as.character()

# Load AA mut OR
df_OR_AA_mut <- read_tsv('../../results/mutation_enrichment/df_AA_mut_enrichment_saltational_opportunity_CV_1.8_BF_20.tsv') %>% 
  filter(threshold_for_inclusion == 2)  %>% 
  mutate(aa_muts_tot = paste0(gene, ':', mutation)) %>% 
  select(OR, OR_lower, OR_upper, aa_muts_tot, p_value_BH, n_mut_of_interest_adapt_branches) %>% 
  mutate(is_signif = p_value_BH < 0.05) %>% 
  mutate(is_within_profile_persistent = aa_muts_tot %in% mut_profile_persistent,
         is_within_profile_deers = aa_muts_tot %in% mut_profile_deers,
         is_within_profile_minks = aa_muts_tot %in% mut_profile_minks) %>% 
  mutate(n_profiles = as.numeric(is_within_profile_persistent) + as.numeric(is_within_profile_deers) + as.numeric(is_within_profile_minks)) %>% 
  mutate(type = case_when(n_profiles > 1 ~ 'Multiple profiles',
                          is_within_profile_persistent ~ 'Persistent infections-like',
                          is_within_profile_deers ~ 'Deer-like',
                          is_within_profile_minks ~ 'Mink-like',
                          T ~ 'Others')) %>% 
  mutate(type = factor(type, levels = c('Persistent infections-like', 'Deer-like', 'Mink-like', 'Others'))) %>% 
  arrange(desc(type))




df_OR_AA_mut %>% filter(p_value_BH < 0.05) %>% 
  group_by(type) %>% summarise(n_occ = n())

plt_OR_func_mut_frequency_coloured_by_profile_type <- df_OR_AA_mut %>% 
  ggplot(aes(x = n_mut_of_interest_adapt_branches, y = OR, 
             colour = type, alpha = is_signif)) +
  geom_point() +
  geom_text_repel(data = df_OR_AA_mut %>% 
                    filter(is_signif, type != 'Others'),
                  aes(label = aa_muts_tot),
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
  scale_colour_manual(breaks = c('Persistent infections-like', 'Deer-like', 'Mink-like', 'Others'),
                      values = c(col_persistent, col_deers, col_minks, col_others),
                      name = 'Mutation profile') +
  scale_alpha_manual(breaks = c(T, F), values = c(0.9, .01), guide = 'none') +
  theme_classic() +
  theme(axis.text = element_text(size = 11)) +
  coord_cartesian(ylim = c(1, NA)) +
  guides(colour = guide_legend(override.aes = list(label = '')))

pdf('../figure_mutation_enrichment/scatter_plot_OR_func_mut_count_saltational_colour_mut_profile.pdf',
    height = 4., width = 7.)
plot(plt_OR_func_mut_frequency_coloured_by_profile_type)
dev.off()
png('../figure_mutation_enrichment/scatter_plot_OR_func_mut_count_saltational_colour_mut_profile.png',
    height = 4., width = 7., res = 350, units = 'in')
plot(plt_OR_func_mut_frequency_coloured_by_profile_type)
dev.off()
 