library(tidyverse)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(binom)
source('utils_plot_chronic.R')
source('config_plot.R')

#############################
### LOADING MODEL RESULTS ###
#############################

# Load recurrent mutations in deers / minks / persistent infections
frequent_mutations_summary <- get_frequent_mutations_summary('../../data/frequent_mutations/frequent_mutations_minks_deers.csv',
                                                             '../../data/frequent_mutations/frequent_mutations_chronic.csv')

# Load saltational branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
                                        header = F) %>% unlist() %>% as.character()

# Load enrichment saltational branches
df_OR_aa_muts <- read_tsv('../../results/mutation_enrichment/df_AA_mut_enrichment_saltational_opportunity_CV_1.8_BF_20.tsv') |>
  mutate(full_mut = paste0(gene, ':', mutation)) |> 
  filter(threshold_for_inclusion == 2) |> 
  select(full_mut, OR, OR_lower, OR_upper, p_value, p_value_BH, n_mut_of_interest_adapt_branches)
  

df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') |> 
  filter(is_coded, aa_initial != aa_new) |> 
  mutate(aa_muts_tot = paste0(gene, ':', aa_muts_position, aa_new)) |> 
  left_join(frequent_mutations_summary) |> 
  filter(! is.na(type)) |> 
  mutate(is_salt = node_id %in% name_saltational_branches) |> 
  filter(is_salt) |> 
  group_by(aa_muts_tot) |> 
  summarise(n_occ_salt = n()) |> 
  full_join(frequent_mutations_summary) |> 
  mutate(n_occ_salt = replace_na(n_occ_salt, 0),
         occur_in_salt = n_occ_salt > 1) |> 
  group_by(type) |> 
  summarise(n_occ_in_salt = sum(occur_in_salt),
            n_muts = n()) |> 
  mutate(prop_occ_in_salt = n_occ_in_salt/n_muts)

df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') |> 
  filter(is_coded, aa_initial != aa_new) |> 
  mutate(aa_muts_tot = paste0(gene, ':', aa_muts_position, aa_new)) |> 
  left_join(frequent_mutations_summary) |> 
  mutate(type = replace_na(type, '')) |> 
  mutate(is_salt = node_id %in% name_saltational_branches)

vec_types <- c('Persistent', 'Deers', 'Minks')

df_OR_by_type <- Reduce('bind_rows', lapply(vec_types, FUN = function(curr_type){
  
  # Get OR of muts from profile in salt vs non-salt branches
  tmp_df <- df_aa_muts |> 
    mutate(is_type = type == curr_type) |> 
    group_by(is_type, is_salt) |> summarise(n_muts = n()) |> 
    ungroup()
  
  n_muts_type_salt <- tmp_df$n_muts[tmp_df$is_salt & tmp_df$is_type]
  n_muts_type_nonsalt <- tmp_df$n_muts[!tmp_df$is_salt & tmp_df$is_type]
  n_muts_not_type_salt <- tmp_df$n_muts[tmp_df$is_salt & !tmp_df$is_type]
  n_muts_not_type_nonsalt <- tmp_df$n_muts[!tmp_df$is_salt & !tmp_df$is_type]
  
  mat_fisher <- matrix(
    c(n_muts_type_salt, n_muts_type_nonsalt,
      n_muts_not_type_salt, n_muts_not_type_nonsalt),
    ncol = 2, nrow = 2
    )
  f_test <- fisher.test(mat_fisher)

  # Get number of saltational branches that have at least
  # one mutation from profile
  df_nb_salt_branches_with_type <- df_aa_muts |> 
    mutate(is_type = type == curr_type) |> 
    filter(is_salt) |> 
    group_by(node_id) |> summarise(has_type = sum(is_type) > 0) |> 
    group_by(has_type) |> summarise(n_salt_branches = n())
    
  # Get number of mutations from profile across all saltational branches
  df_n_muts_type_per_salt_branch <- df_aa_muts |> 
    mutate(is_type = type == curr_type) |> 
    filter(is_salt) |> 
    group_by(node_id) |> 
    summarise(n_muts_type = sum(is_type))
    
  # Get number of mutations from profile that are observed at least once
  df_summary_mut_occurence <- frequent_mutations_summary |> 
    filter(type == curr_type) |> 
    left_join(df_OR_aa_muts, by = c('aa_muts_tot' = 'full_mut')) |> 
    mutate(is_enriched = case_when(is.na(p_value_BH) ~ F,
                                   T ~ p_value_BH < 0.05)) |> 
    select(aa_muts_tot, is_enriched) |> 
    left_join(df_aa_muts |> mutate(is_type = type == curr_type) |> 
                 filter(is_salt, is_type) |> group_by(aa_muts_tot) |> 
                summarise(n_occ = n())) |> 
    mutate(n_occ = replace_na(n_occ, 0)) |> 
    summarise(n_enriched = sum(is_enriched),
              n_with_occ = sum(n_occ > 0),
              n_with_occ_non_enriched = sum(n_occ > 0 & !is_enriched),
              n_no_occ = sum(n_occ == 0))
  
  tibble(type = curr_type,
         OR = as.numeric(f_test$estimate),
         OR_lower = as.numeric(f_test$conf.int[1]),
         OR_upper = as.numeric(f_test$conf.int[2]),
         n_muts_salt = n_muts_type_salt,
         n_muts_nonsalt = n_muts_type_nonsalt,
         n_salt_branch_with_mut_from_type = df_nb_salt_branches_with_type$n_salt_branches[df_nb_salt_branches_with_type$has_type],
         n_salt_branch_without_mut_from_type = df_nb_salt_branches_with_type$n_salt_branches[!df_nb_salt_branches_with_type$has_type],
         n_muts_enriched = df_summary_mut_occurence$n_enriched,
         n_muts_with_occ = df_summary_mut_occurence$n_with_occ,
         n_muts_with_occ_non_enriched = df_summary_mut_occurence$n_with_occ_non_enriched,
         n_no_occ = df_summary_mut_occurence$n_no_occ)
}))


panel_overlap_mutations_across_different_profiles <- (df_OR_by_type|> 
  ggplot(aes(x = type, y = n_muts_salt, fill = type)) +
  geom_bar(stat = 'identity')  +
  geom_text(aes(label = n_muts_salt), vjust = -0.5) +
  scale_x_discrete(name = '', breaks = c('Deers', 'Minks', 'Persistent'),
                   labels = c('Deer', 'Mink', 'Persistent')) +
   scale_y_continuous(name = 'Number of mutations from profile\nacross saltational branches',
                      expand = expansion(mult = c(0., 0.1))) +
  scale_fill_manual(values = c(col_deers, col_minks, col_persistent), guide = 'none') +
  theme_classic() +

  df_OR_by_type|> 
  ggplot(aes(x = type, y = n_salt_branch_with_mut_from_type, fill = type)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = n_salt_branch_with_mut_from_type), vjust = -0.5) +
  scale_x_discrete(name = '', breaks = c('Deers', 'Minks', 'Persistent'),
                   labels = c('Deer', 'Mink', 'Persistent'))+
  scale_y_continuous(name = 'Number of saltational branches\nwith at least 1 mutation from profile',
                     expand = expansion(mult = c(0., 0.1))) +
  scale_fill_manual(values = c(col_deers, col_minks, col_persistent), guide = 'none') +
  theme_classic()) +
  
df_OR_by_type |> 
  select(type, n_muts_enriched, n_muts_with_occ_non_enriched, n_no_occ) |> 
  group_by(type) |> 
  pivot_longer(cols = 2:4, values_to = 'count') |> 
  group_by(type) |> mutate(prop = count/sum(count)) |> 
  mutate(name = factor(name, levels = rev(c('n_muts_enriched', 'n_muts_with_occ_non_enriched', 'n_no_occ')))) |> 
  ggplot(aes(x = type, y = prop, fill = name)) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(name = '', breaks = c('Deers', 'Minks', 'Persistent'),
                   labels = c('Deer', 'Mink', 'Persistent')) +
  scale_y_continuous(name = 'Proportion of recurrent mutations',
                     expand = expansion(mult = c(0., 0.05))) +
  scale_fill_manual(values = c('orange2', 'red3', 'gray22'),
                    breaks = c('n_muts_enriched', 'n_muts_with_occ_non_enriched', 'n_no_occ'),
                    labels = c('Enriched on\nsaltational branches',
                               'Not enriched but occurs\non saltational branches',
                               "Doesn't occur on\nsaltational branches"),
                    name = 'Recurrent mutation\ntype') +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.title.position = 'top') +
  plot_annotation(tag_levels = c('A')) & 
  theme(plot.tag = element_text(face = "bold", size = 15))

# Save results
png('../figure_mutation_enrichment/panel_overlap_recurrent_mutations_saltational_by_profile.png',
    height = 4.5, width = 7, res = 350, units = 'in')
plot(panel_overlap_mutations_across_different_profiles)
dev.off()
