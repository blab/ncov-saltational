### This script reproduces figures illustrating how the 
### frequency at which saltational branches have descendants
### (Figure 5A and 5B)

library(tidyverse)
library(scales)
source('config_plot.R')

#############################
### LOADING MODEL RESULTS ###
#############################

# Load name of saltational branches
name_saltational_branches <- 
  read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
             header = F) %>% unlist() %>% as.character()

# Load number of descendants by node
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv')

# Load inferred geographic locations
inferred_countries <- read_tsv('../../results/bte/inferred_countries.tsv') %>% select(strain, inferred_location) %>% 
  rename(inferred_country = inferred_location)
inferred_subregions <- read_tsv('../../results/bte/inferred_subregions.tsv') %>% select(strain, inferred_location) %>% 
  rename(inferred_subregion = inferred_location)
inferred_continents <- read_tsv('../../results/bte/inferred_continents.tsv') %>% select(strain, inferred_location) %>% 
  rename(inferred_continent = inferred_location)

# Load metadata 
metadata_taxonium <- read_tsv('../../results/saltational_branches/metadata_taxonium_CV_1.8_BF_20.tsv') %>% 
  mutate(subdivision_1 = ifelse(subdivision_1 == 'nan', NA, subdivision_1),
         subdivision_1 = ifelse(subdivision_1 == 'Freestate', 'Free State', subdivision_1),
         subdivision_2 = ifelse(subdivision_2 == 'nan', NA, subdivision_2))

# Load child parent relationships (list of all nodes)
df_child_parent <- read_tsv('../../results/bte/df_parent_child_relationship.tsv') %>% 
  mutate(bool_child_is_terminal = ! startsWith(strain, 'node_'))

# Among saltational branches with descendants that only spread to one country, get the ones with descendants with distinct metadata
df_descendants_single_country <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  mutate(bool_spread_to_multiple_people = n_terminal_descendants > 1,
         is_saltational = strain %in% name_saltational_branches) %>% 
  filter(is_saltational, bool_spread_to_multiple_people) %>% 
  left_join(inferred_countries) %>% filter(inferred_country != 'Multiple descendant locations')

df_multiple_distinct_metadata <- Reduce('bind_rows', lapply(1:nrow(df_descendants_single_country), FUN = function(i_row){
  start_node <- df_descendants_single_country$strain[i_row]
  list_descendants_from_start_node <- c(start_node)
  while(length(list_descendants_from_start_node) != (df_descendants_single_country$n_descendants[i_row] + 1)){
    list_descendants_from_start_node <- c(list_descendants_from_start_node,
                                          df_child_parent %>% filter(parent_name %in% list_descendants_from_start_node) %>% 
                                            select(strain) %>% unlist() %>% as.character()
    )
    list_descendants_from_start_node <- unique(list_descendants_from_start_node)
  }
  metadata_taxonium %>% filter(strain %in% list_descendants_from_start_node) %>% 
    summarise(n_distinct_ages = length(unique(age[! is.na(age)])),
              n_distinct_sexes = length(unique(sex[! is.na(sex)])),
              n_distinct_subdivision_1 =  length(unique(subdivision_1[! is.na(subdivision_1)])),
              n_distinct_subdivision_2 = length(unique(subdivision_1[! is.na(subdivision_1)]))) %>% 
    mutate(i_row = i_row, 
           start_node = start_node)
})) %>% 
  mutate(has_descendant_distinct_metadata = (n_distinct_ages > 1 | n_distinct_sexes > 1 | n_distinct_subdivision_1 > 1 | n_distinct_subdivision_2 > 1)) %>% 
  select(start_node, has_descendant_distinct_metadata)


####################
### MAKE FIGURES ###
####################
# Color used for the plot
col_plot <- 'gray22'

# Figure 5A - Distribution of the number of descendants per saltational branches
barplot_descendants_counts_dist <- df_descendants %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(n_terminal_descendants = ifelse(n_terminal_descendants == 1, 0, n_terminal_descendants),
         n_terminal_descendants_crop = ifelse(n_terminal_descendants > 10, 10, n_terminal_descendants)) %>%
  group_by(n_terminal_descendants_crop) %>% summarise(n_branches = n()) %>% 
  ggplot(aes(x = as.factor(n_terminal_descendants_crop), y = n_branches)) +
  geom_bar(stat = 'identity', fill = col_plot) +
  scale_x_discrete(breaks = c(0, 2:10), labels = c(0, 2:9, '10+'),
                   name = 'Number of descendants') +
  scale_y_continuous(name = 'Number of saltational branches', limits = c(0., NA),
                     expand = expansion(mult = c(0., 0.05)),
                     breaks = seq(0, 1200, 200)) +
  theme_classic() +
  theme(axis.text = element_text(size = 11))

png('../figure_descendant_analysis/barplot_descendants_counts_dist.png', 
    height = 3., width = 3.2, res = 350, units = 'in')
plot(barplot_descendants_counts_dist)
dev.off()
pdf('../figure_descendant_analysis/barplot_descendants_counts_dist.pdf', 
    height = 3., width = 3.2)
plot(barplot_descendants_counts_dist)
dev.off()

# Get proportion of branches without descendants
df_descendants %>%
  filter(strain %in% name_saltational_branches) %>%
  mutate(n_terminal_descendants = ifelse(n_terminal_descendants == 1, 0, n_terminal_descendants),
         n_terminal_descendants_crop = ifelse(n_terminal_descendants > 10, 10, n_terminal_descendants)) %>%
  group_by(n_terminal_descendants_crop) %>% summarise(n_branches = n()) %>%
  summarise(prop_with_desc = sum(n_branches[n_terminal_descendants_crop == 0])/sum(n_branches))


# Figure 5B - Proportion of branches with descendants with different metrics 
# of "cluster spread" (looking at clusters with descendants observed in different geographic regions)
df_desc_across_scales <- df_descendants %>% 
  mutate(bool_spread_to_multiple_people = n_terminal_descendants > 1,
         is_saltational = strain %in% name_saltational_branches) %>% 
  filter(is_saltational) %>% 
  select(strain, bool_spread_to_multiple_people, is_saltational) %>% 
  left_join(inferred_countries) %>% left_join(inferred_subregions) %>% left_join(inferred_continents) %>% 
  mutate(bool_spread_to_multiple_countries = ((! is.na(inferred_country)) & (inferred_country == 'Multiple descendant locations')),
         bool_spread_to_multiple_subregions = ((! is.na(inferred_subregion)) & (inferred_subregion == 'Multiple descendant locations')),
         bool_spread_to_multiple_continents = ((! is.na(inferred_continent)) & (inferred_continent == 'Multiple descendant locations'))) %>% 
  left_join(df_multiple_distinct_metadata, by = c('strain' = 'start_node')) %>% 
  mutate(bool_spread_to_multiple_distinct_people = case_when(!bool_spread_to_multiple_people ~ F,
                                                             bool_spread_to_multiple_countries ~ T,
                                                             has_descendant_distinct_metadata ~ T,
                                                             !has_descendant_distinct_metadata ~ F)) %>% 
  select(- has_descendant_distinct_metadata) %>% 
  select(strain, bool_spread_to_multiple_people, bool_spread_to_multiple_distinct_people, bool_spread_to_multiple_countries, bool_spread_to_multiple_subregions, bool_spread_to_multiple_continents)

#write_tsv(df_desc_across_scales, '../figure_descendant_analysis/df_desc_across_scales.tsv')

# Get proportion of descendants across scales
df_prop_desc_across_scales <- df_desc_across_scales %>% 
  pivot_longer(cols = 2:6, values_to = 'bool_spread', names_prefix = 'bool_spread_to_multiple_', names_to = 'scale') %>% 
  group_by(scale) %>% summarise(n_spread = sum(bool_spread), n_tot = n()) %>% 
  group_by(scale) %>% 
  mutate(lower_prop = prop.test(x = n_spread, n = n_tot)$conf.int[1],
         upper_prop = prop.test(x = n_spread, n = n_tot)$conf.int[2],
         prop = n_spread/n_tot)

# Make Figure 5B
plt_prop_with_descendants_across_scales <- df_prop_desc_across_scales %>% 
  ggplot(aes(x = reorder(scale, n_spread), y = n_spread/n_tot)) +
  geom_point(color = col_plot) +
  geom_linerange(aes(ymin = lower_prop, ymax = upper_prop), color = col_plot) +
  geom_text(aes(label = n_spread), vjust = - .5, color = col_plot) +
  coord_flip() +
  scale_x_discrete(name = '', 
                   breaks = c('people', 'distinct_people', 'countries', 'subregions', 'continents'),
                   labels = c('Any\nsamples', 'Distinct\nmetadata', 'In different\ncountries', 'In different\nsubregions', 'In different\ncontinents')) +
  scale_y_continuous(name = 'Proportion of saltational\nbranches with descendants', limits = c(0., NA),
                     expand = expansion(mult = c(0., 0.1))) +
  theme_classic() +
  theme(axis.text = element_text(size = 11),
        axis.title.y = element_blank())

png('../figure_descendant_analysis/prop_branches_with_desc_across_scales.png', 
    height = 3., width = 3., res = 350, units = 'in')
plot(plt_prop_with_descendants_across_scales)
dev.off()
pdf('../figure_descendant_analysis/prop_branches_with_desc_across_scales.pdf', 
    height = 3., width = 3.)
plot(plt_prop_with_descendants_across_scales)
dev.off()
