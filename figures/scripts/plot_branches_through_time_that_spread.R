library(tidyverse)
library(binom)

# Load adaptive branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.7_BF_20.tsv', 
                                        header = F) %>% unlist() %>% as.character()

# Load name of all branches
df_child_parent <- read_tsv('../../results/bte/df_parent_child_relationship.tsv') %>% 
  filter(parent_name != 'None') %>% 
  mutate(is_adaptive = strain %in% name_saltational_branches) %>% 
  filter(is_adaptive)

# Load chronumental data for tips with missing year and internal nodes
chronumental_internal_nodes <- read_tsv('../../results/chronumental/chronumental_inferred_dates.tsv') %>% 
  filter(strain %in% name_saltational_branches, startsWith(strain, 'node_')) %>% 
  mutate(date_chronumental = as.Date(predicted_date),
         year_chronumental = as.Date(paste0(substr(date_chronumental, start = 1, stop = 4), '-01-01'))) %>% 
  select(- predicted_date)

# Use SRA dates if available and chronumental dates for internal nodes and tips with missing dates
df_branches_with_dates <- df_child_parent %>% 
  left_join(metadata %>% select(strain, year, bool_missing_year), by = 'strain') %>% 
  left_join(chronumental %>% select(strain, date_chronumental, year_chronumental), by = 'strain') %>% 
  mutate(bool_missing_year = replace_na(bool_missing_year, T),
         year_final = as.Date(ifelse(bool_missing_year, as.character(year_chronumental), as.character(year)))) %>% 
  mutate(is_adaptive = strain %in% name_saltational_branches,
         is_internal = startsWith(strain, 'node_'))

# Load descendants information
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches)

df_child_parent %>% 
  filter(startsWith(strain, 'node_')) %>% 
  left_join(df_descendants) %>% 
  select(strain, branch_length, n_terminal_descendants) %>% 
  left_join(chronumental_internal_nodes) %>% 
  arrange(date_chronumental) %>% 
  filter(n_terminal_descendants >= 20, branch_length >= 7, branch_length <= 10) %>% View()


# Interesting examples:
# node_861367: With descendants first detected at the end of 2020
# node_9379
# node_38879
# node_176
