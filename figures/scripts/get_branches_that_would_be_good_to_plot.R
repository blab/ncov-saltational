library(tidyverse)

# Load child parent relationships
df_child_parent <- read_tsv('../../results/bte/df_parent_child_relationship.tsv')

# Load adaptive branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
                                        header = F) %>% unlist() %>% as.character()

# Load chronic-like infections
name_muts_chronic_like <- read_csv('../../data/frequent_mutations/frequent_mutations_chronic.csv') %>% 
  filter(source != 'Ghafari et al. Lancet ID. Occurred 4+ times') %>% 
  mutate(aa_mutation_full = paste0(gene, ':', aa_mut_direction)) %>% 
  select(aa_mutation_full) %>% unlist() %>% as.character()

# Get characteristics of saltational branches
df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  group_by(node_id, gene, aa_muts) %>% 
  mutate(aa_mutation_full = paste0(gene, ':', aa_muts_position, aa_new),
         is_chronic_like = aa_mutation_full %in% name_muts_chronic_like) %>% 
  filter(is_coded, node_id %in% name_saltational_branches) %>% 
  group_by(node_id) %>% 
  summarise(n_muts = n(), n_non_syn = sum(aa_initial != aa_new), 
            n_muts_S = sum(gene == 'S'), n_non_syn_S = sum((aa_initial != aa_new) & (gene == 'S')),
            n_non_syn_chronic = sum((aa_initial != aa_new) & (is_chronic_like))) %>% 
  rename(strain = node_id)

# Get descendants information
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches)

# Load inferred location information
df_inferred_location <- read_tsv('../../results/bte/inferred_countries.tsv') %>% 
  filter(strain %in% name_saltational_branches) %>% 
  select(strain, inferred_location) %>% rename(country = inferred_location)
df_inferred_continent <- read_tsv('../../results/bte/inferred_continents.tsv') %>% 
  filter(strain %in% name_saltational_branches) %>% 
  select(strain, inferred_location) %>% rename(continent = inferred_location)


## Identify nodes that would be nice to plot
name_nodes_to_plot <- c('node_176', 'node_2803', 'node_8768', 'node_9379', 'node_38879', 'node_39826', 'node_812479', 'node_810966', 'node_841663', 'node_861367')

df_child_parent %>% filter(strain %in% name_nodes_to_plot) %>% 
  select(strain, parent_name, branch_length)

## Objectives - Get one for each year
## Illustrate spread within a country - Focus on US with many states
## Illustrate spread within a region (For Europe: node_655327)
## Illustrate spread within a worldwide
## Illustrate that many chronic-like infection on some of the branches leading to it
## 

## 2020

## 2021
# node_2803 - US small cluster that spread to 5 states - 1 chronic like mut

## 2022
# node_812479 - lineage defining for BS.1 - spread in the USA, UK, Denmark, Thaildand - 5 chronic like infections, 7 NS in S, 18 NS
# Spread for 5 months - closest related sample in the UK from March 2022 (6 months before) - 26 terminal descendants

## 2023

## 2024



# Little cluster with many relevant mutations observed in Angola + UK (node_39826)



# node_176 - 16 muts, 12 NS, 5 NS in Spike (1 chronic. Spread over 3 months in US, USA, Austria, Germany and France)


# node_39264 - 1 chronic like muts, 10 NS Spike muts, spread in USA + Mexico. 

# node_810966 - 4 chronic like muts! spread in USA, Denmark + UK (countries with the most samples). Many many descendants.
# Closest other sequence - sequence for the USA from May 2022 (6 months before first descendant is observed)
# Lineage defining for BA.2.3.20

# node_8768 - multiple spillover?


# node_841663 - chronic infections from Germany (multiple from the same male of 48)





# 


df_aa_muts %>% 
  left_join(df_descendants) %>% 
  left_join(df_inferred_location) %>% left_join(df_inferred_continent) %>% 
  filter(n_terminal_descendants > 1) %>% View()

df_aa_muts %>% 
  left_join(df_descendants) %>% 
  left_join(df_inferred_location) %>% left_join(df_inferred_continent) %>% filter(country == 'USA', n_terminal_descendants > 1) %>% View()

df_aa_muts %>% 
  left_join(df_descendants) %>% 
  left_join(df_inferred_location) %>% 
  filter(n_terminal_descendants >= 30)
