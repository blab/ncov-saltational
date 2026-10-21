library(tidyverse)
library(scales)
source('config_plot.R')

#############################
### LOADING MODEL RESULTS ###
#############################
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

# Among branches with descendants that only spread to one country, get the ones with descendants with distinct metadata
df_descendants_single_country <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  mutate(bool_spread_to_multiple_people = n_terminal_descendants > 1) %>% 
  filter(bool_spread_to_multiple_people) %>% 
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
              n_distinct_subdivision_2 = length(unique(subdivision_1[! is.na(subdivision_1)]))) %>% 
    mutate(i_row = i_row, 
           start_node = start_node)
})) %>% 
  mutate(has_descendant_distinct_metadata = (n_distinct_ages > 1 | n_distinct_sexes > 1 | n_distinct_subdivision_1 > 1 | n_distinct_subdivision_2 > 1)) %>% 
  dplyr::select(start_node, has_descendant_distinct_metadata)


##############
vec_BF <- c('5', '10', '20', '40', '80')

df_descendants_char_across_threshold <- bind_rows(
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    name_salt_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_', curr_BF, '.tsv'), 
                                     header = F) %>% unlist() %>% as.character()
    
    
    df_descendants |> 
      filter(strain %in% name_salt_branches) |> 
      left_join(inferred_countries) |> 
      mutate(bool_spread_to_multiple_people = n_terminal_descendants > 1,
             bool_spread_to_multiple_countries = (inferred_country == 'Multiple descendant locations'),
             bool_spread_to_multiple_countries = replace_na(bool_spread_to_multiple_countries, F)) |> 
      summarise(n_with_desc = sum(bool_spread_to_multiple_people),
                n_with_desc_multiple_countries = sum(bool_spread_to_multiple_countries),
                n_tot = n()) |> 
      mutate(lower_prop = prop.test(x = n_with_desc, n = n_tot)$conf.int[1],
             upper_prop = prop.test(x = n_with_desc, n = n_tot)$conf.int[2],
             lower_prop_multiple_countries = prop.test(x = n_with_desc_multiple_countries, n = n_tot)$conf.int[1],
             upper_prop_multiple_countries = prop.test(x = n_with_desc_multiple_countries, n = n_tot)$conf.int[2]
             ) |> 
      mutate(BF = curr_BF,
             nonsyn_count_percentile = 0.995)
    
  })),
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    name_salt_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_perc_99.9_CV_1.8_BF_', curr_BF, '.tsv'), 
                                     header = F) %>% unlist() %>% as.character()
    
    
    df_descendants |> 
      filter(strain %in% name_salt_branches) |> 
      left_join(inferred_countries) |> 
      mutate(bool_spread_to_multiple_people = n_terminal_descendants > 1,
             bool_spread_to_multiple_countries = (inferred_country == 'Multiple descendant locations'),
             bool_spread_to_multiple_countries = replace_na(bool_spread_to_multiple_countries, F)) |> 
      summarise(n_with_desc = sum(bool_spread_to_multiple_people),
                n_with_desc_multiple_countries = sum(bool_spread_to_multiple_countries),
                n_tot = n()) |> 
      mutate(lower_prop = prop.test(x = n_with_desc, n = n_tot)$conf.int[1],
             upper_prop = prop.test(x = n_with_desc, n = n_tot)$conf.int[2],
             lower_prop_multiple_countries = prop.test(x = n_with_desc_multiple_countries, n = n_tot)$conf.int[1],
             upper_prop_multiple_countries = prop.test(x = n_with_desc_multiple_countries, n = n_tot)$conf.int[2]
      ) |> 
      mutate(BF = curr_BF,
             nonsyn_count_percentile = 0.999)
    
  })),
  Reduce('bind_rows', lapply(vec_BF, FUN = function(curr_BF){
    name_salt_branches <- read.table(paste0('../../results/saltational_branches/name_saltational_branches_perc_99_CV_1.8_BF_', curr_BF, '.tsv'), 
                                     header = F) %>% unlist() %>% as.character()
    
    
    df_descendants |> 
      filter(strain %in% name_salt_branches) |> 
      left_join(inferred_countries) |> 
      mutate(bool_spread_to_multiple_people = n_terminal_descendants > 1,
             bool_spread_to_multiple_countries = (inferred_country == 'Multiple descendant locations'),
             bool_spread_to_multiple_countries = replace_na(bool_spread_to_multiple_countries, F)) |> 
      summarise(n_with_desc = sum(bool_spread_to_multiple_people),
                n_with_desc_multiple_countries = sum(bool_spread_to_multiple_countries),
                n_tot = n()) |> 
      mutate(lower_prop = prop.test(x = n_with_desc, n = n_tot)$conf.int[1],
             upper_prop = prop.test(x = n_with_desc, n = n_tot)$conf.int[2],
             lower_prop_multiple_countries = prop.test(x = n_with_desc_multiple_countries, n = n_tot)$conf.int[1],
             upper_prop_multiple_countries = prop.test(x = n_with_desc_multiple_countries, n = n_tot)$conf.int[2]
      )|> 
      mutate(BF = curr_BF,
             nonsyn_count_percentile = 0.99)
    
  }))
)


## Display results
highlight_labels <- function(x, pat, color = "black", family = "") {
  ifelse(grepl(pat, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}

impact_thresholds_prop_branches_with_desc <- df_descendants_char_across_threshold %>% 
  mutate(BF_threshold_label = case_when(
    (nonsyn_count_percentile == 0.995 & BF == "20") ~ "**20**",
    TRUE ~ as.character(BF)),
    BF_threshold_label = factor(BF_threshold_label, c('5', '10', '20', '**20**', '40', '80'))) %>%
  ggplot(aes(x = BF_threshold_label)) +
  geom_point(aes(y = n_with_desc/n_tot, shape = 'Any detected descendants')) +
  geom_point(aes(y = n_with_desc_multiple_countries/n_tot, shape = 'Descendants in multiple countries')) +
  geom_linerange(aes(ymin = lower_prop, ymax = upper_prop), alpha = 0.5) +
  geom_linerange(aes(ymin = lower_prop_multiple_countries, ymax = upper_prop_multiple_countries), alpha = 0.5) +
  facet_wrap(.~ nonsyn_count_percentile, 
             scales = 'free_x',
             labeller = as_labeller(function(x) highlight_labels(paste0('Percentile nonsyn count: ', scales::percent(as.numeric(x))),
                                                                 'Percentile nonsyn count: 99.5%'))) +
  scale_x_discrete(name = 'Bayes Factor threshold') +
  scale_y_continuous(name = 'Proportion of saltational branches\nwith detected descendants',
                     breaks = seq(0., 1., 0.1),
                     limits = c(0., 0.6)
                     ) +
  scale_shape_manual(name = '', values = c(16, 17)) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        strip.text.x = element_markdown(),
        axis.text.x = element_markdown(),
        legend.position = 'bottom')

impact_thresholds_nb_branches_with_desc <- df_descendants_char_across_threshold %>% 
  mutate(BF_threshold_label = case_when(
    (nonsyn_count_percentile == 0.995 & BF == "20") ~ "**20**",
    TRUE ~ as.character(BF)),
    BF_threshold_label = factor(BF_threshold_label, c('5', '10', '20', '**20**', '40', '80'))) %>%
  ggplot(aes(x = BF_threshold_label)) +
  geom_bar(aes(y = n_with_desc),
           stat = 'identity', fill = 'darkcyan') +
  geom_text(aes(y = n_with_desc, label = n_with_desc),
            vjust = -.1) +
  facet_wrap(.~ nonsyn_count_percentile, 
             scales = 'free_x',
             labeller = as_labeller(function(x) highlight_labels(paste0('Percentile nonsyn count: ', scales::percent(as.numeric(x))),
                                                                 'Percentile nonsyn count: 99.5%'))) +
  scale_x_discrete(name = 'Bayes Factor threshold') +
  scale_y_continuous(name = 'Number of saltational branches\nwith detected descendants',
                     expand = expansion(mult = c(0.01, 0.1)),
                     breaks = seq(0., 8000, 2000)) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        strip.text.x = element_markdown(),
        axis.text.x = element_markdown(),
        legend.position = 'top')

impact_thresholds_nb_branches_with_desc_multiple_countries <- df_descendants_char_across_threshold %>% 
  mutate(BF_threshold_label = case_when(
    (nonsyn_count_percentile == 0.995 & BF == "20") ~ "**20**",
    TRUE ~ as.character(BF)),
    BF_threshold_label = factor(BF_threshold_label, c('5', '10', '20', '**20**', '40', '80'))) %>%
  ggplot(aes(x = BF_threshold_label)) +
  geom_bar(aes(y = n_with_desc_multiple_countries),
           stat = 'identity', fill = 'darkcyan') +
  geom_text(aes(y = n_with_desc_multiple_countries, label = n_with_desc_multiple_countries),
            vjust = -.1) +
  facet_wrap(.~ nonsyn_count_percentile, 
             scales = 'free_x',
             labeller = as_labeller(function(x) highlight_labels(paste0('Percentile nonsyn count: ', scales::percent(as.numeric(x))),
                                                                 'Percentile nonsyn count: 99.5%'))) +
  scale_x_discrete(name = 'Bayes Factor threshold') +
  scale_y_continuous(name = 'Number of saltational branches\nwith detected descendants\nin multiple countries',
                     expand = expansion(mult = c(0.01, 0.1)),
                     #breaks = seq(0., 8000, 2000)
                     ) +
  scale_shape_manual(name = '', values = c(16, 17)) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        strip.text.x = element_markdown(),
        axis.text.x = element_markdown(),
        legend.position = 'top')

panel_impact_threshold_descendants <-
  (impact_thresholds_nb_branches_with_desc / impact_thresholds_nb_branches_with_desc_multiple_countries / impact_thresholds_prop_branches_with_desc) &
  theme(plot.tag = element_text(face = "bold", size = 16)) & 
  plot_annotation(tag_levels = 'A') &
  theme(plot.margin = unit(c(.01,.01,.01,.01), "cm"))

png('../figure_model_illustration/impact_threshold_char_descendants.png',
    height = 7.9, width = 6.9, res = 350, units = 'in')
plot(panel_impact_threshold_descendants)
dev.off()
