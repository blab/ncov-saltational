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
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv')%>% 
  mutate(n_terminal_descendants = ifelse(n_terminal_descendants == 1, 0, n_terminal_descendants),
         n_terminal_descendants_crop = ifelse(n_terminal_descendants > 10, 10, n_terminal_descendants)) 

# Load inferred geographic locations
inferred_countries <- read_tsv('../../results/bte/inferred_countries.tsv') %>% 
  dplyr::select(strain, inferred_location) %>% 
  rename(inferred_country = inferred_location)
inferred_subregions <- read_tsv('../../results/bte/inferred_subregions.tsv') %>% 
  dplyr::select(strain, inferred_location) %>% 
  rename(inferred_subregion = inferred_location)
inferred_continents <- read_tsv('../../results/bte/inferred_continents.tsv') %>% 
  dplyr::select(strain, inferred_location) %>% 
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
                                            dplyr::select(strain) %>% unlist() %>% as.character()
    )
    list_descendants_from_start_node <- unique(list_descendants_from_start_node)
  }
  metadata_taxonium %>% filter(strain %in% list_descendants_from_start_node) %>% 
    summarise(n_distinct_ages = length(unique(age[! is.na(age)])),
              n_distinct_sexes = length(unique(sex[! is.na(sex)])),
              n_distinct_subdivision_1 =  length(unique(subdivision_1[! is.na(subdivision_1)])),
              n_distinct_subdivision_2 = length(unique(subdivision_1[! is.na(subdivision_1)])),
              n_descendants_with_age_info = sum(! is.na(age)),,
              n_descendants_with_sex_info = sum(! is.na(sex)),
              n_descendants_with_subdivision_1_info = sum(! is.na(subdivision_1)),
              n_descendants_with_subdivision_2_info = sum(! is.na(subdivision_2))) %>% 
    mutate(i_row = i_row, 
           start_node = start_node,
           n_descendants = length(list_descendants_from_start_node))
})) %>% 
  mutate(has_descendant_distinct_metadata = (n_distinct_ages > 1 | n_distinct_sexes > 1 | n_distinct_subdivision_1 > 1 | n_distinct_subdivision_2 > 1)) %>% 
  dplyr::select(start_node, has_descendant_distinct_metadata)


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
  summarise(prop_with_desc = sum(n_branches[n_terminal_descendants_crop == 0])/sum(n_branches),
            n_tot_branches = sum(n_branches),
            n_branches_with_desc = sum(n_branches[n_terminal_descendants_crop == 0]))
prop.test(x = 1079, n = 1535, conf.level = 0.95)

# Figure 5B - Proportion of branches with descendants with different metrics 
# of "cluster spread" (looking at clusters with descendants observed in different geographic regions)
df_desc_across_scales <- df_descendants %>% 
  mutate(bool_spread_to_multiple_people = n_terminal_descendants > 1,
         is_saltational = strain %in% name_saltational_branches) %>% 
  filter(is_saltational) %>% 
  dplyr::select(strain, bool_spread_to_multiple_people, is_saltational) %>% 
  left_join(inferred_countries) %>% left_join(inferred_subregions) %>% left_join(inferred_continents) %>% 
  mutate(bool_spread_to_multiple_countries = ((! is.na(inferred_country)) & (inferred_country == 'Multiple descendant locations')),
         bool_spread_to_multiple_subregions = ((! is.na(inferred_subregion)) & (inferred_subregion == 'Multiple descendant locations')),
         bool_spread_to_multiple_continents = ((! is.na(inferred_continent)) & (inferred_continent == 'Multiple descendant locations'))) %>% 
  left_join(df_multiple_distinct_metadata, by = c('strain' = 'start_node')) %>% 
  mutate(bool_spread_to_multiple_distinct_people = case_when(!bool_spread_to_multiple_people ~ F,
                                                             bool_spread_to_multiple_countries ~ T,
                                                             has_descendant_distinct_metadata ~ T,
                                                             !has_descendant_distinct_metadata ~ F)) %>% 
  dplyr::select(- has_descendant_distinct_metadata) %>% 
  dplyr::select(strain, bool_spread_to_multiple_people, bool_spread_to_multiple_distinct_people, bool_spread_to_multiple_countries, bool_spread_to_multiple_subregions, bool_spread_to_multiple_continents)

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

###########
metadata_taxonium |> 
  summarise(n_with_sex = sum(! is.na(sex)),
            n_with_age = sum(! is.na(age)),
            n_with_country = sum(! is.na(country)),
            n_with_subdivision = sum(! is.na(subdivision_full)),
            n_tot_branches = n()) |> 
  mutate(prop_with_sex = n_with_sex/n_tot_branches,
         prop_with_age = n_with_age/n_tot_branches,
         prop_with_country = n_with_country/n_tot_branches,
         prop_with_subdivision = n_with_subdivision/n_tot_branches)

df_assess_missing_metadata <- Reduce('bind_rows', lapply(1:nrow(df_descendants_single_country), FUN = function(i_row){
  start_node <- df_descendants_single_country$strain[i_row]
  list_descendants_from_start_node <- c(start_node)
  while(length(list_descendants_from_start_node) != (df_descendants_single_country$n_descendants[i_row] + 1)){
    list_descendants_from_start_node <- c(list_descendants_from_start_node,
                                          df_child_parent %>% filter(parent_name %in% list_descendants_from_start_node) %>% 
                                            dplyr::select(strain) %>% unlist() %>% as.character()
    )
    list_descendants_from_start_node <- unique(list_descendants_from_start_node)
  }
  metadata_taxonium %>% filter(strain %in% list_descendants_from_start_node) %>% 
    summarise(n_distinct_ages = length(unique(age[! is.na(age)])),
              n_distinct_sexes = length(unique(sex[! is.na(sex)])),
              n_distinct_subdivision_1 =  length(unique(subdivision_1[! is.na(subdivision_1)])),
              n_distinct_subdivision_2 = length(unique(subdivision_1[! is.na(subdivision_1)])),
              n_descendants_with_age_info = sum(! is.na(age)),,
              n_descendants_with_sex_info = sum(! is.na(sex)),
              n_descendants_with_subdivision_1_info = sum(! is.na(subdivision_1)),
              n_descendants_with_subdivision_2_info = sum(! is.na(subdivision_2))) %>% 
    mutate(i_row = i_row, 
           start_node = start_node,
           n_descendants = length(list_descendants_from_start_node))
})) %>% 
  mutate(has_descendant_distinct_metadata = (n_distinct_ages > 1 | n_distinct_sexes > 1 | n_distinct_subdivision_1 > 1 | n_distinct_subdivision_2 > 1))


df_assess_missing_metadata |> 
  filter(n_descendants_with_age_info <= 1 &
           n_descendants_with_sex_info <= 1 &
           n_descendants_with_subdivision_1_info <= 1 &
           n_descendants_with_subdivision_2_info <= 1) |> View()
  
df_assess_missing_metadata |> 
  filter(n_descendants_with_age_info == n_descendants |
           n_descendants_with_sex_info == n_descendants |
           n_descendants_with_subdivision_1_info == n_descendants |
           n_descendants_with_subdivision_2_info == n_descendants) 

df_assess_missing_metadata |>
  ggplot(aes(x = n_distinct_ages > 1, y = n_descendants_with_age_info)) +
  geom_boxplot(aes(group = n_distinct_ages > 1))


plt_proba_distinct_age_metadata <- df_assess_missing_metadata |> 
  ggplot(aes(x = n_descendants_with_age_info + 1, 
             y = as.integer(n_distinct_ages > 1))) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              color = 'darkcyan', fill = 'darkcyan') +
  geom_point(data = df_assess_missing_metadata |> 
               group_by(n_descendants_with_age_info, n_distinct_ages) |> 
               summarise(n_branches = n()),
             aes(size = n_branches)) +
  scale_x_continuous(trans = 'log',
                     name = 'Number of descendants with age information + 1',
                     breaks = c(seq(1, 9, 1), seq(10, 90, 10), 100),
                     labels = c(1, 2, '', '', 5, rep('', 4),
                                10, 20, '', '', 50, rep('', 4),
                                100)) +
  scale_y_continuous(breaks = seq(0., 1., 0.2),
                     name = 'Proportion with descendants\nhaving distinct age metadata') +
  scale_size_continuous(name = 'Number of branches') +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

plt_proba_distinct_sex_metadata <- df_assess_missing_metadata |> 
  ggplot(aes(x = n_descendants_with_sex_info + 1, 
             y = as.integer(n_distinct_sexes > 1))) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              color = 'darkcyan', fill = 'darkcyan') +
  geom_point(data = df_assess_missing_metadata |> 
               group_by(n_descendants_with_sex_info, n_distinct_sexes) |> 
               summarise(n_branches = n()),
             aes(size = n_branches)) +
  scale_x_continuous(trans = 'log',
                     name = 'Number of descendants with sex information + 1',
                     breaks = c(seq(1, 9, 1), seq(10, 90, 10), 100),
                     labels = c(1, 2, '', '', 5, rep('', 4),
                                10, 20, '', '', 50, rep('', 4),
                                100)) +
  scale_y_continuous(breaks = seq(0., 1., 0.2),
                     name = 'Proportion with descendants\nhaving distinct sex metadata') +
  scale_size_continuous(name = 'Number of branches') +
  theme_bw() +
  theme(panel.grid.minor = element_blank())


panel_proba_descendants <- plt_proba_distinct_age_metadata / plt_proba_distinct_sex_metadata +
  plot_annotation(tag_levels = c('A')) & 
  theme(plot.tag = element_text(face = 'bold', size = 15))


png('../figure_descendant_analysis/panel_prop_descendants_with_distinct_metadata_func_availability.png',
    height = 5., width = 6., res = 350, units = 'in')
plot(panel_proba_descendants)
dev.off()
