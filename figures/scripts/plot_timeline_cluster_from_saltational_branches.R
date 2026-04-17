library(tidyverse)
library(viridis)
# Load name of saltational branches
vec_name_node_saltational <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv') %>% unlist() %>% as.character()

# Load number of descendants
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% vec_name_node_saltational) %>% 
  mutate(earliest_date_descendants = as.Date(earliest_date_descendants),
         latest_date_descendant = as.Date(latest_date_descendant))

# Load inferred countries information
inferred_countries <- read_tsv('../../results/bte/inferred_countries.tsv') %>% select(strain, inferred_location) %>% 
  rename(inferred_country = inferred_location)
inferred_subregions <- read_tsv('../../results/bte/inferred_subregions.tsv') %>% select(strain, inferred_location) %>% 
  rename(inferred_subregion = inferred_location)
inferred_continents <- read_tsv('../../results/bte/inferred_continents.tsv') %>% select(strain, inferred_location) %>% 
  rename(inferred_continent = inferred_location)

# Load AA muts information
df_non_syn_counts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  filter(is_coded, node_id %in% vec_name_node_saltational) %>% 
  group_by(node_id) %>% summarise(n_non_syn = sum(aa_initial != aa_new), n_muts = n()) %>% 
  rename(strain = node_id)

df_descendants %>% 
  left_join(inferred_countries) %>%
  left_join(inferred_subregions) %>% left_join(inferred_continents) %>%  
  left_join(df_non_syn_counts) %>% 
  filter(! is.na(earliest_date_descendants), ! is.na(latest_date_descendant)) %>%
  filter(earliest_date_descendants < as.Date('2021-01-01')) %>% 
  mutate(cluster_duration = as.numeric(latest_date_descendant - earliest_date_descendants, units = 'days')) %>% 
  filter(n_terminal_descendants > 1) %>% View()
  # summarise(median_cluster_duration = median(cluster_duration, na.rm = T),
  #           mean_cluster_duration = mean(cluster_duration, na.rm = T),
  #           lower_IQR_cluster_duration = quantile(cluster_duration, probs = 0.25, na.rm = T),
  #           upper_IQR_cluster_duration = quantile(cluster_duration, probs = 0.75, na.rm = T))
  ggplot(aes(y = reorder(strain, earliest_date_descendants), 
             colour = n_non_syn,
             #colour = n_terminal_descendants
             )) +
  geom_segment(aes(x = earliest_date_descendants, 
                   xend = latest_date_descendant,
                   yend = strain)) +
  geom_point(aes(x = earliest_date_descendants)) +
  scale_y_discrete(labels = NULL, name = 'Cluster ID') +
  scale_colour_viridis_b() 
  facet_wrap(. ~ inferred_continent, scales = 'free_y')
