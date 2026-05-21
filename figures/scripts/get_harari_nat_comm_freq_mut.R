library(tidyverse)

df_frequent_mut_metadata_cluster <- read_csv('../../data/chronic-associated-mutations/harari-top-mutation-fitness.csv') %>% 
  filter(gene != 'Non-coding', startsWith(full_mutation, gene))

df_frequent_mut_metadata_cluster %>% 
  group_by(full_mutation, gene) %>% 
  summarise(n_samples = sum(n_samples), 
            mean_explainability = mean(explainability), 
            mean_reliability = mean(reliability)) %>% 
  filter(n_samples >= 3) %>% arrange(- n_samples) %>% View()
