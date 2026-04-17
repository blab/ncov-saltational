library(tidyverse)
library(scales)
library(binom)
source('config_plot.R')

# Load saltational branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', header = F) %>%
  unlist() %>% as.character()

# Load mutation information
df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  filter(is_coded) %>% rename(strain = node_id) %>% 
  group_by(strain) %>% summarise(n_non_syn = sum(aa_initial != aa_new))

# Load descendants information
df_prop_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  mutate(has_descendants = n_terminal_descendants > 1,
         is_saltational = strain %in% name_saltational_branches) %>% 
  select(strain, has_descendants, is_saltational) %>% 
  left_join(df_aa_muts, by = 'strain') %>% 
  mutate(n_non_syn = replace_na(n_non_syn, 0)) %>% 
  mutate(has_more_5_NS_mut = n_non_syn >= 5) %>% 
  group_by(is_saltational, has_more_5_NS_mut) %>% 
  summarise(n_with_desc =  sum(has_descendants), 
            n_branches = n(),
            prop_with_desc = n_with_desc/n_branches) %>% 
  group_by(is_saltational, has_more_5_NS_mut) %>% 
  mutate(lower_prop = prop.test(x = n_with_desc, n = n_branches)$conf.int[1],
         upper_prop = prop.test(x = n_with_desc, n = n_branches)$conf.int[2]) %>% 
  ungroup()


# Barplot of proportion of saltational branches with descendants
barplot_prop_branches_with_desc <- df_prop_descendants %>% 
  ggplot(aes(x = interaction(has_more_5_NS_mut, is_saltational), 
             y = prop_with_desc, 
             fill = is_saltational)) +
  geom_bar(stat = 'identity') +
  geom_linerange(aes(ymin = lower_prop, ymax = upper_prop)) +
  scale_x_discrete(breaks = c('TRUE.TRUE', 'TRUE.FALSE', 'FALSE.FALSE'),
                   name = '',
                   labels = c('Saltational', 'Non saltational\nmany muts', 'Non saltational\nfew muts')) +
  scale_y_continuous(name = 'Proportion of branches with descendants', limits = c(0., NA),
                     expand = expansion(mult = c(0., 0.05))) +
  scale_fill_manual(breaks = c(T, F), values = c(col_adaptive, col_non_adaptive),
                    guide = 'none') +
  theme_classic() +
  theme(axis.text = element_text(size = 11)) +
  coord_flip()

png('../figure_descendant_analysis/prop_branches_with_desc_saltational.png', height = 2.3, width = 4.4,
    res = 350, units = 'in')
plot(barplot_prop_branches_with_desc)
dev.off()
pdf('../figure_descendant_analysis/prop_branches_with_desc_saltational.pdf', height = 2.3, width = 4.4,)
plot(barplot_prop_branches_with_desc)
dev.off()

# Number of descendnats per branches
barplot_descendants_counts_dist <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(n_terminal_descendants = ifelse(n_terminal_descendants == 1, 0, n_terminal_descendants),
         n_terminal_descendants_crop = ifelse(n_terminal_descendants > 10, 10, n_terminal_descendants)) %>%
  group_by(n_terminal_descendants_crop) %>% summarise(n_branches = n()) %>% 
  ggplot(aes(x = as.factor(n_terminal_descendants_crop), y = n_branches)) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(breaks = c(0, 2:10), labels = c(0, 2:9, '10+'),
                   name = 'Number of descendants') +
  scale_y_continuous(name = 'Number of saltational branches', limits = c(0., NA),
                     expand = expansion(mult = c(0., 0.05)),
                     breaks = seq(0, 1200, 200)) +
  theme_classic() +
  theme(axis.text = element_text(size = 11))
  
png('../figure_descendant_analysis/barplot_descendants_counts_dist.png', 
    height = 3., width = 4.4, res = 350, units = 'in')
plot(barplot_descendants_counts_dist)
dev.off()
pdf('../figure_descendant_analysis/barplot_descendants_counts_dist.pdf', 
    height = 3., width = 4.4)
plot(barplot_descendants_counts_dist)
dev.off()



#### Rather inferring branches that spread as branches that were detected in multiple countries
inferred_countries <- read_tsv('../../results/bte/inferred_countries.tsv') %>% select(strain, inferred_location) %>% 
  rename(inferred_country = inferred_location)
inferred_subregions <- read_tsv('../../results/bte/inferred_subregions.tsv') %>% select(strain, inferred_location) %>% 
  rename(inferred_subregion = inferred_location)
inferred_continents <- read_tsv('../../results/bte/inferred_continents.tsv') %>% select(strain, inferred_location) %>% 
  rename(inferred_continent = inferred_location)

plt_prop_with_descendants_across_scales <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  mutate(bool_spread_to_multiple_people = n_terminal_descendants > 1,
         is_saltational = strain %in% name_saltational_branches) %>% 
  filter(is_saltational) %>% 
  select(strain, bool_spread_to_multiple_people, is_saltational) %>% 
  left_join(inferred_countries) %>% left_join(inferred_subregions) %>% left_join(inferred_continents) %>% 
  mutate(bool_spread_to_multiple_countries = (inferred_country == 'Multiple descendant locations'),
         bool_spread_to_multiple_subregions = (inferred_subregion == 'Multiple descendant locations'),
         bool_spread_to_multiple_continents = (inferred_continent == 'Multiple descendant locations')) %>% 
  select(strain, bool_spread_to_multiple_people, bool_spread_to_multiple_countries, bool_spread_to_multiple_subregions, bool_spread_to_multiple_continents) %>% 
  pivot_longer(cols = 2:5, values_to = 'bool_spread', names_prefix = 'bool_spread_to_multiple_', names_to = 'scale') %>% 
  group_by(scale) %>% summarise(n_spread = sum(bool_spread), n_tot = n()) %>% 
  ggplot(aes(x = reorder(scale, -n_spread), y = n_spread/n_tot)) +
  geom_bar(stat = 'identity', fill = 'gray22') +
  geom_text(aes(label = n_spread), vjust = - .5) +
  scale_x_discrete(name = '', 
                   breaks = c('people', 'countries', 'subregions', 'continents'),
                   labels = c('Overall', 'Different\ncountries', 'Different\nsubregions', 'Different\ncontinents')) +
  scale_y_continuous(name = 'Proportion of branches with descendants', limits = c(0., NA),
                     expand = expansion(mult = c(0., 0.1))) +
  theme_classic() +
  theme(axis.text = element_text(size = 11))
  
  
png('../figure_descendant_analysis/prop_branches_with_desc_across_scales.png', height = 3.5, width = 4.,
    res = 350, units = 'in')
plot(plt_prop_with_descendants_across_scales)
dev.off()
pdf('../figure_descendant_analysis/prop_branches_with_desc_across_scales.pdf', height = 3.5, width = 4.,)
plot(plt_prop_with_descendants_across_scales)
dev.off()