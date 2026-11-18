library(tidyverse)
library(argparse)
library(parallel)

# Enrichment of AA mutations on saltational branches, computed separately for saltational
# events occurring on internal and on terminal branches (reviewer comment 3).
#
# The background is always all non-saltational branches, as in the published analysis: only
# the saltational row of the Fisher table is restricted to a branch type. The denominators
# come from the branch-type split of the mutation opportunity counts.
#
# The mutations tested are those passing the inclusion threshold over ALL saltational
# branches, so that the same mutations are tested in both strata and the odds ratios are
# directly comparable.

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--df_aa_muts', type = 'character', help = 'TSV file with AA muts info')
  parser$add_argument('--name_adaptive_branches', type = 'character', help = 'Text file with name of adaptive branches')
  parser$add_argument('--df_mutation_opportunity', type = 'character', help = 'TSV file with mutation opportunity counts split by branch type')
  parser$add_argument('--output_df_enrichment', type = 'character', help = 'TSV file for enrichment OR scores')
  return(parser$parse_args())
}

args <- collect_args()

vec_branch_types <- c('all', 'internal', 'terminal')
vec_thresholds_for_inclusion <- 1:5

# Load input data
df_aa_muts <- read_tsv(args$df_aa_muts) %>%
  rename(strain = node_id) %>%
  # Keep non-synonymous only
  filter(aa_new != aa_initial) %>%
  # Keep coded ones only
  filter(is_coded) %>%
  mutate(mut_direction = paste0(aa_muts_position, aa_new))

vec_name_node_signif <- read.table(args$name_adaptive_branches) %>% unlist() %>% as.character()

df_mut_count_accounting_opportunity <- read_tsv(args$df_mutation_opportunity)

# Internal nodes are named node_XXX, terminal ones carry a sample name
is_internal_branch <- function(strain) grepl('^node_', strain)

get_OR_AA_mut_by_branch_type <- function(df_aa_muts, df_mut_count_accounting_opportunity, vec_name_node_signif,
                                         threshold_mut_occ_for_inclusion = 1, branch_type = 'all'){

  # Restrict the saltational branches to the branch type of interest. The background
  # (is_adaptive == FALSE) keeps all non-saltational branches whatever the branch type.
  vec_name_node_signif_branch_type <- switch(branch_type,
                                             'all' = vec_name_node_signif,
                                             'internal' = vec_name_node_signif[is_internal_branch(vec_name_node_signif)],
                                             'terminal' = vec_name_node_signif[! is_internal_branch(vec_name_node_signif)])

  # Opportunity column matching the branch type
  col_saltational_opportunity <- switch(branch_type,
                                        'all' = 'n_saltational_opportunity',
                                        'internal' = 'n_saltational_opportunity_internal',
                                        'terminal' = 'n_saltational_opportunity_terminal')

  # Get number of mutation by branch type. Mutations on saltational branches of another
  # branch type are dropped entirely, so they count neither as saltational nor as background.
  df_aa_muts_branch_type <- df_aa_muts %>%
    filter(! strain %in% setdiff(vec_name_node_signif, vec_name_node_signif_branch_type))

  df_aa_muts_counts_by_branch_type <- df_aa_muts_branch_type %>%
    mutate(is_adaptive = strain %in% vec_name_node_signif_branch_type) %>%
    group_by(is_adaptive, mut_direction, gene) %>%
    summarise(n_mutations = n())

  # Mutations to test: those passing the threshold over ALL saltational branches, so that
  # the same set is tested in every stratum
  df_mutated_pos <- df_aa_muts %>%
    filter(strain %in% vec_name_node_signif) %>%
    group_by(gene, mut_direction) %>%
    summarise(n_mutations = n(), .groups = 'drop') %>%
    filter(n_mutations >= threshold_mut_occ_for_inclusion) %>%
    select(gene, mut_direction) %>%
    mutate(mut_all = paste0(gene, ':', mut_direction)) %>%
    left_join(df_mut_count_accounting_opportunity %>%
                select(mutation, all_of(col_saltational_opportunity), n_non_saltational_opportunity),
              by = c('mut_all' = 'mutation')) %>%
    rename(n_saltational_opportunity_branch_type = all_of(col_saltational_opportunity))

  df_OR_by_mut <- Reduce('bind_rows', mclapply(1:nrow(df_mutated_pos), FUN = function(i_row){
    curr_gene <- df_mutated_pos$gene[i_row]
    curr_pos <- df_mutated_pos$mut_direction[i_row]
    curr_n_saltational_opportunity <- df_mutated_pos$n_saltational_opportunity_branch_type[i_row]
    curr_n_non_saltational_opportunity <- df_mutated_pos$n_non_saltational_opportunity[i_row]

    df_muts_counts <- df_aa_muts_counts_by_branch_type %>%
      mutate(is_mut = (gene == curr_gene) & (mut_direction == curr_pos)) %>%
      group_by(is_mut, is_adaptive) %>% summarise(n_mutations = sum(n_mutations))

    n_mut_of_interest_adapt_branches <- df_muts_counts$n_mutations[df_muts_counts$is_adaptive == T & df_muts_counts$is_mut == T]
    n_mut_of_interest_normal_branches <- df_muts_counts$n_mutations[df_muts_counts$is_adaptive == F & df_muts_counts$is_mut == T]

    n_mut_of_interest_adapt_branches <- ifelse(length(n_mut_of_interest_adapt_branches) == 0, 0, n_mut_of_interest_adapt_branches)
    n_mut_of_interest_normal_branches <- ifelse(length(n_mut_of_interest_normal_branches) == 0, 0, n_mut_of_interest_normal_branches)

    n_other_muts_adapt_branches <- curr_n_saltational_opportunity - n_mut_of_interest_adapt_branches
    n_other_muts_normal_branches <- curr_n_non_saltational_opportunity - n_mut_of_interest_normal_branches

    n_other_muts_adapt_branches <- ifelse(length(n_other_muts_adapt_branches) == 0, 0, n_other_muts_adapt_branches)
    n_other_muts_normal_branches <- ifelse(length(n_other_muts_normal_branches) == 0, 0, n_other_muts_normal_branches)

    mat_fisher <- matrix(c(n_mut_of_interest_adapt_branches, n_mut_of_interest_normal_branches, n_other_muts_adapt_branches, n_other_muts_normal_branches),
                         nrow = 2, ncol = 2, byrow = T)

    f_test <- fisher.test(x = mat_fisher)

    df_fisher <- tibble(OR = as.numeric(f_test$estimate),
                        OR_lower = as.numeric(f_test$conf.int[1]),
                        OR_upper = as.numeric(f_test$conf.int[2]),
                        p_value = f_test$p.value,
                        n_mut_of_interest_adapt_branches = n_mut_of_interest_adapt_branches,
                        n_mut_of_interest_normal_branches = n_mut_of_interest_normal_branches,
                        n_other_muts_adapt_branches = n_other_muts_adapt_branches,
                        n_other_muts_normal_branches = n_other_muts_normal_branches,
                        gene = curr_gene,
                        mutation = curr_pos)
  }, mc.cores = detectCores() - 1))

  df_OR_by_mut$p_value_BH <- p.adjust(df_OR_by_mut$p_value, method = 'BH')
  df_OR_by_mut$p_value_bonf <- p.adjust(df_OR_by_mut$p_value, method = 'bonferroni')

  return(df_OR_by_mut)
}

# Run for each branch type and each inclusion threshold
df_OR_by_AA_mut <- Reduce('bind_rows', lapply(vec_branch_types, FUN = function(curr_branch_type){
  Reduce('bind_rows', lapply(vec_thresholds_for_inclusion, FUN = function(curr_threshold){
    get_OR_AA_mut_by_branch_type(df_aa_muts, df_mut_count_accounting_opportunity, vec_name_node_signif,
                                 curr_threshold, curr_branch_type) %>%
      mutate(threshold_for_inclusion = curr_threshold,
             branch_type = curr_branch_type)
  }))
}))

# Test whether the OR differs between internal and terminal saltational branches (Woolf test
# on the difference of log ORs). With a shared background this is a direct test of whether the
# mutation occurs at a different rate on the two sets of saltational branches.
get_log_OR_and_variance <- function(df){
  df %>%
    # Haldane-Anscombe correction, so that the test is defined when a cell is empty
    mutate(a = n_mut_of_interest_adapt_branches + 0.5,
           b = n_other_muts_adapt_branches + 0.5,
           c = n_mut_of_interest_normal_branches + 0.5,
           d = n_other_muts_normal_branches + 0.5,
           log_OR = log(a * d / (b * c)),
           var_log_OR = 1/a + 1/b + 1/c + 1/d) %>%
    select(gene, mutation, threshold_for_inclusion, log_OR, var_log_OR)
}

df_heterogeneity <- full_join(
  get_log_OR_and_variance(df_OR_by_AA_mut %>% filter(branch_type == 'internal')) %>%
    rename(log_OR_internal = log_OR, var_log_OR_internal = var_log_OR),
  get_log_OR_and_variance(df_OR_by_AA_mut %>% filter(branch_type == 'terminal')) %>%
    rename(log_OR_terminal = log_OR, var_log_OR_terminal = var_log_OR),
  by = c('gene', 'mutation', 'threshold_for_inclusion')
) %>%
  mutate(log_OR_ratio_internal_over_terminal = log_OR_internal - log_OR_terminal,
         se_log_OR_ratio = sqrt(var_log_OR_internal + var_log_OR_terminal),
         p_value_internal_vs_terminal = 2 * pnorm(-abs(log_OR_ratio_internal_over_terminal / se_log_OR_ratio))) %>%
  group_by(threshold_for_inclusion) %>%
  mutate(p_value_internal_vs_terminal_BH = p.adjust(p_value_internal_vs_terminal, method = 'BH')) %>%
  ungroup() %>%
  select(gene, mutation, threshold_for_inclusion, log_OR_ratio_internal_over_terminal,
         se_log_OR_ratio, p_value_internal_vs_terminal, p_value_internal_vs_terminal_BH)

df_OR_by_AA_mut <- df_OR_by_AA_mut %>%
  left_join(df_heterogeneity, by = c('gene', 'mutation', 'threshold_for_inclusion'))

# Save results
write_tsv(df_OR_by_AA_mut, args$output_df_enrichment)
