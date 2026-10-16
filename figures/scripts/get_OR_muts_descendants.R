#############################
### LOADING MODEL RESULTS ###
#############################

# Load name of saltational branches
name_saltational_branches <- 
  read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
             header = F) %>% unlist() %>% as.character()

# Get mutation enriched on saltational branches among the ones that occur at least twice
enriched_AA_mut <- read_tsv('../../results/mutation_enrichment/df_AA_mut_enrichment_saltational_opportunity_CV_1.8_BF_20.tsv') %>% 
  filter(threshold_for_inclusion == 2) %>% filter(p_value_BH < 0.05) %>% select(gene, mutation, OR) %>% 
  mutate(mut_full = paste0(gene, ':', mutation))

# Only retain descendant information from saltational branches
df_descendants <-read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(has_descendants = n_terminal_descendants > 1) %>% 
  select(strain, has_descendants, n_terminal_descendants)

# Load AA mutation information
df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  rename(strain = node_id) %>% 
  filter(is_coded, strain %in% name_saltational_branches) %>% 
  mutate(is_synonymous = (aa_initial == aa_new),
         aa_muts_direction = paste0(aa_muts_position, aa_new)) %>% 
  filter(! is_synonymous) %>% 
  left_join(df_descendants) |> 
  mutate(is_saltational = strain %in% name_saltational_branches,
         mut_full = paste0(gene, ':', aa_muts_direction))
    
# List of enriched mutations     
vec_enriched_muts <- unique(enriched_AA_mut$mut_full)

# Get OR of a saltational branch having descendants depending on whether that branch has a mutation of interest
df_OR_descendants <- Reduce('bind_rows', lapply(vec_enriched_muts, FUN = function(curr_mut){
  curr_df <- df_aa_muts |> 
    mutate(is_mut = mut_full == curr_mut) |> 
    group_by(strain, has_descendants, is_saltational) |> 
    summarise(has_mut = sum(is_mut)) |> 
    group_by(has_mut, has_descendants) |> 
    summarise(n_branches = n())
  
  first_or_zero <- function(x) {
    if (length(x) == 0L) 0 else x[[1L]] + 0
  }
  
  ## Fisher test
  mat_fisher <- matrix(
    c(first_or_zero(curr_df$n_branches[curr_df$has_descendants & curr_df$has_mut]),
      first_or_zero(curr_df$n_branches[!curr_df$has_descendants & curr_df$has_mut]),
      first_or_zero(curr_df$n_branches[curr_df$has_descendants & !curr_df$has_mut]),
      first_or_zero(curr_df$n_branches[!curr_df$has_descendants & !curr_df$has_mut])),
    nrow = 2, ncol = 2, byrow = T
  )
  f_test <- fisher.test(mat_fisher)
  
  ## Mann-Whitney U test
  curr_df_wilcox <- df_descendants |> 
    left_join(df_aa_muts |> 
                filter(mut_full == curr_mut) |> 
                select(strain) |> 
                mutate(has_mut = T)) |> 
    mutate(has_mut = replace_na(has_mut, F)) |> 
    mutate(n_terminal_descendants = ifelse(n_terminal_descendants == 1, 0, n_terminal_descendants))
  
  w_test <- wilcox.test(data = curr_df_wilcox, n_terminal_descendants ~ has_mut)
  
  tibble(mut = curr_mut, 
         OR_fisher = f_test$estimate,
         OR_fisher_lower = as.numeric(f_test$conf.int[1]),
         OR_fisher_upper = as.numeric(f_test$conf.int[2]),
         p_value_fisher = f_test$p.value,
         statistic_wilcox = w_test$statistic,
         p_value_wilcox = w_test$p.value,
         n_salt_with_mut_with_desc = first_or_zero(curr_df$n_branches[curr_df$has_descendants & curr_df$has_mut]),
         n_salt_with_mut_without_desc = first_or_zero(curr_df$n_branches[!curr_df$has_descendants & curr_df$has_mut]),
         n_salt_without_mut_with_desc = first_or_zero(curr_df$n_branches[curr_df$has_descendants & !curr_df$has_mut]),
         n_salt_without_mut_without_desc = first_or_zero(curr_df$n_branches[!curr_df$has_descendants & !curr_df$has_mut])
         )
}))
df_OR_descendants$p_value_fisher_bonf <- p.adjust(df_OR_descendants$p_value_fisher, method = 'bonferroni')
df_OR_descendants$p_value_fisher_BH <- p.adjust(df_OR_descendants$p_value_fisher, method = 'BH')
df_OR_descendants$p_value_wilcox_bonf <- p.adjust(df_OR_descendants$p_value_wilcox, method = 'bonferroni')
df_OR_descendants$p_value_wilcox_BH <- p.adjust(df_OR_descendants$p_value_wilcox, method = 'BH')

df_OR_descendants |> filter(mut %in% c('ORF1a:3255I', 'S:452R', 'S:681H', 'S:484K', 'E:30I', 'S:376T')) |> View()

df_OR_descendants |> 
  arrange(- n_salt_with_mut_with_desc) |> 
  mutate(OR_for_plot = paste0(round(OR_fisher, 1), ' (', round(OR_fisher_lower, 1), '-', round(OR_fisher_upper, 2), ')')) |> 
  select(mut, n_salt_with_mut_with_desc, n_salt_with_mut_without_desc, OR_for_plot) |> 
  write_csv('../figure_descendant_analysis/df_OR_descendants_muts.csv')

