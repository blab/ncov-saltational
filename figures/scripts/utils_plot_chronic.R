# Loading and processing data
get_frequent_mutations_summary <- function(file_freq_mut_mink_deers, file_freq_mut_chronic){
  # Load frequent mutations in minks / deers
  frequent_mutations_minks_deers <- read_csv(file_freq_mut_mink_deers) %>% 
    filter(! is.na(aa_mut_direction)) %>% 
    mutate(aa_muts_tot = paste0(gene, ':', aa_mut_direction)) %>% 
    select(host, aa_muts_tot) %>% distinct() %>%
    filter(! is.na(aa_muts_tot)) %>% 
    rename(type = host)
  
  # Load frequent mutations in persistent infection
  frequent_mutations_chronic <- read_csv(file_freq_mut_chronic) %>% 
    filter(! is.na(aa_mutation_full)) %>% 
    mutate(aa_muts_tot = paste0(gene, ':', aa_mut_direction)) %>% 
    select(type, aa_muts_tot) %>% distinct()
  
  frequent_mutations_summary <- bind_rows(frequent_mutations_minks_deers, frequent_mutations_chronic) 
 
  return(frequent_mutations_summary)
}

get_OR_AA_mut <- function(file_OR, threshold_inclusion, frequent_mutations_summary){
  mut_profile_persistent <- frequent_mutations_summary %>% filter(type == 'Persistent') %>% select(aa_muts_tot) %>% unlist() %>% as.character()
  mut_profile_minks <- frequent_mutations_summary %>% filter(type == 'Minks') %>% select(aa_muts_tot) %>% unlist() %>% as.character()
  mut_profile_deers <- frequent_mutations_summary %>% filter(type == 'Deers') %>% select(aa_muts_tot) %>% unlist() %>% as.character()
  
  df_OR_AA_mut <- read_tsv(file_OR) %>% 
    filter(threshold_for_inclusion == threshold_inclusion)  %>% 
    mutate(aa_muts_tot = paste0(gene, ':', mutation)) %>% 
    select(OR, OR_lower, OR_upper, aa_muts_tot, p_value_BH, n_mut_of_interest_adapt_branches) %>% 
    mutate(is_signif = p_value_BH < 0.05) %>% 
    mutate(is_within_profile_persistent = aa_muts_tot %in% mut_profile_persistent,
           is_within_profile_deers = aa_muts_tot %in% mut_profile_deers,
           is_within_profile_minks = aa_muts_tot %in% mut_profile_minks) %>% 
    mutate(n_profiles = as.numeric(is_within_profile_persistent) + as.numeric(is_within_profile_deers) + as.numeric(is_within_profile_minks)) %>% 
    mutate(type = case_when(n_profiles > 1 ~ 'Multiple profiles',
                            is_within_profile_persistent ~ 'Persistent infections-like',
                            is_within_profile_deers ~ 'Deer-like',
                            is_within_profile_minks ~ 'Mink-like',
                            T ~ 'Others')) %>% 
    mutate(type = factor(type, levels = c('Persistent infections-like', 'Deer-like', 'Mink-like', 'Others'))) %>% 
    arrange(desc(type))
  
  return(df_OR_AA_mut)
}

get_chronic_muts_count <- function(file_large_chronic_mut_dataset){
  df_chronic_muts <- read_csv(file_large_chronic_mut_dataset) %>% 
    tidyr::separate(aa_mut, into = c("gene", "aa_muts"), sep = ":") %>% 
    group_by(gene, aa_muts) %>% 
    mutate(mutation = substr(aa_muts, start = 2, stop = nchar(aa_muts)),
           aa_initial = substr(aa_muts, start = 1, stop = 1),
           aa_new = substr(aa_muts, start = nchar(aa_muts), stop = nchar(aa_muts))) %>% 
    ungroup() %>% 
    group_by(gene, mutation) %>% 
    summarise(n_mut_chronic = sum(EPCI_count),
              n_mut_high_quality_sequence = sum(HQCS_count)) %>% 
    ungroup()
  
}

compute_OR_muts_chronic_infections_saltational <- function(file_large_chronic_mut_dataset,
                                                           file_aa_muts, 
                                                           name_saltational_branches){
  # Load chronic mutation counts
  df_chronic_muts <- get_chronic_muts_count(file_large_chronic_mut_dataset)
  
  # Get mutations in viridian tree
  df_aa_muts <- read_tsv(file_aa_muts) %>% 
    filter(is_coded) %>% # Coded AA mut
    filter(aa_new != aa_initial) %>% # Non-synonymous
    mutate(full_mut = paste0(aa_muts_position, aa_new),
           is_saltational = node_id %in% name_saltational_branches) %>% 
    group_by(gene, full_mut, is_saltational) %>% 
    summarise(n_muts = n()) %>% 
    group_by(is_saltational) %>% 
    mutate(n_muts_tot = sum(n_muts)) %>% 
    ungroup()
  
  # Get name of mutations observed in viridian tree
  df_muts_in_persistent_inf <- df_aa_muts %>% 
    group_by(gene, full_mut) %>% 
    summarise(is_within_saltational = sum(is_saltational) > 0)
  
  # Get mutations observed in persistent infections dataset or in saltational branches
  # but removing the ones that are never observed in any of the datasets
  df_chronic_muts <- df_chronic_muts %>% 
    full_join(df_muts_in_persistent_inf, by = c('gene', 'mutation' = 'full_mut')) %>% 
    filter(! is.na(is_within_saltational)) %>%   # Remove mutations not observed in the viridian tree
    filter(! is.na(n_mut_chronic)) # Remove mutations not observed in the EPCI dataset
  
  df_chronic_muts_for_OR_comput <- df_chronic_muts %>% 
    filter(is_within_saltational | n_mut_chronic > 0)
  
  # Get OR of mut in saltational branches
  df_modif_OR_saltational <- Reduce('bind_rows', lapply(1:nrow(df_chronic_muts_for_OR_comput), FUN = function(i_row){
    curr_gene <- df_chronic_muts_for_OR_comput$gene[i_row]
    curr_mut <- df_chronic_muts_for_OR_comput$mutation[i_row]
    
    tmp_df <- df_aa_muts %>% 
      filter((full_mut == curr_mut) & (gene == curr_gene))
    
    n_mut_saltational <- tmp_df$n_muts[tmp_df$is_saltational]
    n_mut_saltational <- ifelse(length(n_mut_saltational) == 0, 0, n_mut_saltational)
    
    n_mut_background <- tmp_df$n_muts[(!tmp_df$is_saltational)]
    n_mut_background <- ifelse(length(n_mut_background) == 0, 0, n_mut_background)
    
    n_other_muts_saltational <- tmp_df$n_muts_tot[tmp_df$is_saltational] - n_mut_saltational
    n_other_muts_background <-  tmp_df$n_muts_tot[! tmp_df$is_saltational] - n_mut_background
    
    n_other_muts_saltational <- ifelse(length(n_other_muts_saltational) == 0, df_aa_muts$n_muts_tot[df_aa_muts$is_saltational][1], n_other_muts_saltational)
    n_other_muts_background <-  ifelse(length(n_other_muts_background) == 0, df_aa_muts$n_muts_tot[!df_aa_muts$is_saltational][1], n_other_muts_background)
    
    mat_fisher <- matrix(c(n_mut_saltational, n_other_muts_saltational, n_mut_background, n_other_muts_background),
                         ncol = 2, nrow = 2)
    f_test <- fisher.test(x = mat_fisher)
    
    mat_fisher_modif <- matrix(c(n_mut_saltational + 1, n_other_muts_saltational + 1, 
                                 n_mut_background + 1, n_other_muts_background + 1),
                               ncol = 2, nrow = 2)
    f_test_modif <- fisher.test(x = mat_fisher_modif)
    
    tibble(gene = curr_gene, 
           mutation = curr_mut,
           OR_saltational = f_test$estimate,
           OR_saltational_lower = f_test$conf.int[1],
           OR_saltational_upper = f_test$conf.int[2],
           OR_saltational_modif = f_test_modif$estimate,
           OR_saltational_modif_lower = f_test_modif$conf.int[1],
           OR_saltational_modif_upper = f_test_modif$conf.int[2],
           n_mut_saltational = n_mut_saltational,
           n_mut_background = n_mut_background,
           n_other_muts_saltational = n_other_muts_saltational,
           n_other_muts_background = n_other_muts_background,
           p_value_salt = f_test$p.value)
  }))
  
  # Get OR of mutations in persistent infections
  df_OR_chronic_modif <- Reduce('bind_rows', lapply(1:nrow(df_chronic_muts_for_OR_comput), FUN = function(i_row){
    curr_gene <- df_chronic_muts_for_OR_comput$gene[i_row]
    curr_mut <- df_chronic_muts_for_OR_comput$mutation[i_row]
    
    tmp_df <- df_chronic_muts %>% 
      mutate(is_mut = (mutation == curr_mut) & (gene == curr_gene)) %>% 
      group_by(is_mut) %>% 
      summarise(n_chronic = sum(n_mut_chronic),
                n_high_quality_sequences = sum(n_mut_high_quality_sequence))
    
    n_mut_chronic <- tmp_df$n_chronic[tmp_df$is_mut]
    n_mut_background <- tmp_df$n_high_quality_sequences[tmp_df$is_mut]
    
    n_other_muts_chronic <- tmp_df$n_chronic[!tmp_df$is_mut]
    n_other_muts_background <- tmp_df$n_high_quality_sequences[!tmp_df$is_mut]
    
    mat_fisher <- matrix(c(n_mut_chronic, n_other_muts_chronic, n_mut_background, n_other_muts_background),
                         ncol = 2, nrow = 2)
    f_test <- fisher.test(x = mat_fisher)
    
    mat_fisher_modif <- matrix(c(n_mut_chronic + 1, n_other_muts_chronic + 1,
                                 n_mut_background + 1, n_other_muts_background + 1),
                               ncol = 2, nrow = 2)
    f_test_modif <- fisher.test(x = mat_fisher_modif)
    
    tibble(gene = curr_gene, 
           mutation = curr_mut,
           n_mut_chronic = n_mut_chronic, 
           n_mut_background = n_mut_background,
           n_other_muts_chronic = n_other_muts_chronic,
           n_other_muts_background = n_other_muts_background,
           OR_chronic = f_test$estimate,
           OR_chronic_lower = f_test$conf.int[1],
           OR_chronic_upper = f_test$conf.int[2],
           OR_chronic_modif = f_test_modif$estimate,
           OR_chronic_modif_lower = f_test_modif$conf.int[1],
           OR_chronic_modif_upper = f_test_modif$conf.int[2],
           p_value_chronic = f_test$p.value)
  }))
  
  df_OR_chronic_saltational <- df_chronic_muts %>% select(gene, mutation) %>% 
    left_join(df_OR_chronic_modif %>%
                rename(n_mut_high_quality_sequence = n_mut_background,
                       n_other_muts_high_quality_sequence = n_other_muts_background),
              by = c('gene', 'mutation')) %>% 
    left_join(df_modif_OR_saltational, by = c('gene', 'mutation')) %>% 
    filter(n_mut_saltational + n_mut_background >= 1, 
           n_mut_chronic + n_mut_high_quality_sequence >= 1)
  
  return(df_OR_chronic_saltational)
}
compute_OR_muts_chronic_infections_saltational_internal <- function(file_large_chronic_mut_dataset,
                                                                    file_aa_muts, 
                                                                    name_saltational_branches){
  # Load chronic mutation counts
  df_chronic_muts <- get_chronic_muts_count(file_large_chronic_mut_dataset)
  
  # Get mutations in viridian tree
  df_aa_muts <- read_tsv(file_aa_muts) %>% 
    filter(is_coded) %>% # Coded AA mut
    filter(aa_new != aa_initial) %>% # Non-synonymous
    mutate(full_mut = paste0(aa_muts_position, aa_new),
           is_saltational = node_id %in% name_saltational_branches) %>% 
    filter(! is_saltational | startsWith(node_id, 'node_')) |> 
    group_by(gene, full_mut, is_saltational) %>% 
    summarise(n_muts = n()) %>% 
    group_by(is_saltational) %>% 
    mutate(n_muts_tot = sum(n_muts)) %>% 
    ungroup()
  
  # Get name of mutations observed in viridian tree
  df_muts_in_persistent_inf <- df_aa_muts %>% 
    group_by(gene, full_mut) %>% 
    summarise(is_within_saltational = sum(is_saltational) > 0)
  
  # Get mutations observed in persistent infections dataset or in saltational branches
  # but removing the ones that are never observed in any of the datasets
  df_chronic_muts <- df_chronic_muts %>% 
    full_join(df_muts_in_persistent_inf, by = c('gene', 'mutation' = 'full_mut')) %>% 
    filter(! is.na(is_within_saltational)) %>%   # Remove mutations not observed in the viridian tree
    filter(! is.na(n_mut_chronic)) # Remove mutations not observed in the EPCI dataset
  
  df_chronic_muts_for_OR_comput <- df_chronic_muts %>% 
    filter(is_within_saltational | n_mut_chronic > 0)
  
  # Get OR of mut in saltational branches
  df_modif_OR_saltational <- Reduce('bind_rows', lapply(1:nrow(df_chronic_muts_for_OR_comput), FUN = function(i_row){
    curr_gene <- df_chronic_muts_for_OR_comput$gene[i_row]
    curr_mut <- df_chronic_muts_for_OR_comput$mutation[i_row]
    
    tmp_df <- df_aa_muts %>% 
      filter((full_mut == curr_mut) & (gene == curr_gene))
    
    n_mut_saltational <- tmp_df$n_muts[tmp_df$is_saltational]
    n_mut_saltational <- ifelse(length(n_mut_saltational) == 0, 0, n_mut_saltational)
    
    n_mut_background <- tmp_df$n_muts[(!tmp_df$is_saltational)]
    n_mut_background <- ifelse(length(n_mut_background) == 0, 0, n_mut_background)
    
    n_other_muts_saltational <- tmp_df$n_muts_tot[tmp_df$is_saltational] - n_mut_saltational
    n_other_muts_background <-  tmp_df$n_muts_tot[! tmp_df$is_saltational] - n_mut_background
    
    n_other_muts_saltational <- ifelse(length(n_other_muts_saltational) == 0, df_aa_muts$n_muts_tot[df_aa_muts$is_saltational][1], n_other_muts_saltational)
    n_other_muts_background <-  ifelse(length(n_other_muts_background) == 0, df_aa_muts$n_muts_tot[!df_aa_muts$is_saltational][1], n_other_muts_background)
    
    mat_fisher <- matrix(c(n_mut_saltational, n_other_muts_saltational, n_mut_background, n_other_muts_background),
                         ncol = 2, nrow = 2)
    f_test <- fisher.test(x = mat_fisher)
    
    mat_fisher_modif <- matrix(c(n_mut_saltational + 1, n_other_muts_saltational + 1, 
                                 n_mut_background + 1, n_other_muts_background + 1),
                               ncol = 2, nrow = 2)
    f_test_modif <- fisher.test(x = mat_fisher_modif)
    
    tibble(gene = curr_gene, 
           mutation = curr_mut,
           OR_saltational = f_test$estimate,
           OR_saltational_lower = f_test$conf.int[1],
           OR_saltational_upper = f_test$conf.int[2],
           OR_saltational_modif = f_test_modif$estimate,
           OR_saltational_modif_lower = f_test_modif$conf.int[1],
           OR_saltational_modif_upper = f_test_modif$conf.int[2],
           n_mut_saltational = n_mut_saltational,
           n_mut_background = n_mut_background,
           n_other_muts_saltational = n_other_muts_saltational,
           n_other_muts_background = n_other_muts_background,
           p_value_salt = f_test$p.value)
  }))
  
  # Get OR of mutations in persistent infections
  df_OR_chronic_modif <- Reduce('bind_rows', lapply(1:nrow(df_chronic_muts_for_OR_comput), FUN = function(i_row){
    curr_gene <- df_chronic_muts_for_OR_comput$gene[i_row]
    curr_mut <- df_chronic_muts_for_OR_comput$mutation[i_row]
    
    tmp_df <- df_chronic_muts %>% 
      mutate(is_mut = (mutation == curr_mut) & (gene == curr_gene)) %>% 
      group_by(is_mut) %>% 
      summarise(n_chronic = sum(n_mut_chronic),
                n_high_quality_sequences = sum(n_mut_high_quality_sequence))
    
    n_mut_chronic <- tmp_df$n_chronic[tmp_df$is_mut]
    n_mut_background <- tmp_df$n_high_quality_sequences[tmp_df$is_mut]
    
    n_other_muts_chronic <- tmp_df$n_chronic[!tmp_df$is_mut]
    n_other_muts_background <- tmp_df$n_high_quality_sequences[!tmp_df$is_mut]
    
    mat_fisher <- matrix(c(n_mut_chronic, n_other_muts_chronic, n_mut_background, n_other_muts_background),
                         ncol = 2, nrow = 2)
    f_test <- fisher.test(x = mat_fisher)
    
    mat_fisher_modif <- matrix(c(n_mut_chronic + 1, n_other_muts_chronic + 1,
                                 n_mut_background + 1, n_other_muts_background + 1),
                               ncol = 2, nrow = 2)
    f_test_modif <- fisher.test(x = mat_fisher_modif)
    
    tibble(gene = curr_gene, 
           mutation = curr_mut,
           n_mut_chronic = n_mut_chronic, 
           n_mut_background = n_mut_background,
           n_other_muts_chronic = n_other_muts_chronic,
           n_other_muts_background = n_other_muts_background,
           OR_chronic = f_test$estimate,
           OR_chronic_lower = f_test$conf.int[1],
           OR_chronic_upper = f_test$conf.int[2],
           OR_chronic_modif = f_test_modif$estimate,
           OR_chronic_modif_lower = f_test_modif$conf.int[1],
           OR_chronic_modif_upper = f_test_modif$conf.int[2],
           p_value_chronic = f_test$p.value)
  }))
  
  df_OR_chronic_saltational <- df_chronic_muts %>% select(gene, mutation) %>% 
    left_join(df_OR_chronic_modif %>%
                rename(n_mut_high_quality_sequence = n_mut_background,
                       n_other_muts_high_quality_sequence = n_other_muts_background),
              by = c('gene', 'mutation')) %>% 
    left_join(df_modif_OR_saltational, by = c('gene', 'mutation')) %>% 
    filter(n_mut_saltational + n_mut_background >= 1, 
           n_mut_chronic + n_mut_high_quality_sequence >= 1)
  
  return(df_OR_chronic_saltational)
}

compute_OR_muts_chronic_infections_saltational_terminal <- function(file_large_chronic_mut_dataset,
                                                                    file_aa_muts, 
                                                                    name_saltational_branches){
  # Load chronic mutation counts
  df_chronic_muts <- get_chronic_muts_count(file_large_chronic_mut_dataset)
  
  # Get mutations in viridian tree
  df_aa_muts <- read_tsv(file_aa_muts) %>% 
    filter(is_coded) %>% # Coded AA mut
    filter(aa_new != aa_initial) %>% # Non-synonymous
    mutate(full_mut = paste0(aa_muts_position, aa_new),
           is_saltational = node_id %in% name_saltational_branches) %>% 
    filter(! is_saltational | !startsWith(node_id, 'node_')) |> 
    group_by(gene, full_mut, is_saltational) %>% 
    summarise(n_muts = n()) %>% 
    group_by(is_saltational) %>% 
    mutate(n_muts_tot = sum(n_muts)) %>% 
    ungroup()
  
  # Get name of mutations observed in viridian tree
  df_muts_in_persistent_inf <- df_aa_muts %>% 
    group_by(gene, full_mut) %>% 
    summarise(is_within_saltational = sum(is_saltational) > 0)
  
  # Get mutations observed in persistent infections dataset or in saltational branches
  # but removing the ones that are never observed in any of the datasets
  df_chronic_muts <- df_chronic_muts %>% 
    full_join(df_muts_in_persistent_inf, by = c('gene', 'mutation' = 'full_mut')) %>% 
    filter(! is.na(is_within_saltational)) %>%   # Remove mutations not observed in the viridian tree
    filter(! is.na(n_mut_chronic)) # Remove mutations not observed in the EPCI dataset
  
  df_chronic_muts_for_OR_comput <- df_chronic_muts %>% 
    filter(is_within_saltational | n_mut_chronic > 0)
  
  # Get OR of mut in saltational branches
  df_modif_OR_saltational <- Reduce('bind_rows', lapply(1:nrow(df_chronic_muts_for_OR_comput), FUN = function(i_row){
    curr_gene <- df_chronic_muts_for_OR_comput$gene[i_row]
    curr_mut <- df_chronic_muts_for_OR_comput$mutation[i_row]
    
    tmp_df <- df_aa_muts %>% 
      filter((full_mut == curr_mut) & (gene == curr_gene))
    
    n_mut_saltational <- tmp_df$n_muts[tmp_df$is_saltational]
    n_mut_saltational <- ifelse(length(n_mut_saltational) == 0, 0, n_mut_saltational)
    
    n_mut_background <- tmp_df$n_muts[(!tmp_df$is_saltational)]
    n_mut_background <- ifelse(length(n_mut_background) == 0, 0, n_mut_background)
    
    n_other_muts_saltational <- tmp_df$n_muts_tot[tmp_df$is_saltational] - n_mut_saltational
    n_other_muts_background <-  tmp_df$n_muts_tot[! tmp_df$is_saltational] - n_mut_background
    
    n_other_muts_saltational <- ifelse(length(n_other_muts_saltational) == 0, df_aa_muts$n_muts_tot[df_aa_muts$is_saltational][1], n_other_muts_saltational)
    n_other_muts_background <-  ifelse(length(n_other_muts_background) == 0, df_aa_muts$n_muts_tot[!df_aa_muts$is_saltational][1], n_other_muts_background)
    
    mat_fisher <- matrix(c(n_mut_saltational, n_other_muts_saltational, n_mut_background, n_other_muts_background),
                         ncol = 2, nrow = 2)
    f_test <- fisher.test(x = mat_fisher)
    
    mat_fisher_modif <- matrix(c(n_mut_saltational + 1, n_other_muts_saltational + 1, 
                                 n_mut_background + 1, n_other_muts_background + 1),
                               ncol = 2, nrow = 2)
    f_test_modif <- fisher.test(x = mat_fisher_modif)
    
    tibble(gene = curr_gene, 
           mutation = curr_mut,
           OR_saltational = f_test$estimate,
           OR_saltational_lower = f_test$conf.int[1],
           OR_saltational_upper = f_test$conf.int[2],
           OR_saltational_modif = f_test_modif$estimate,
           OR_saltational_modif_lower = f_test_modif$conf.int[1],
           OR_saltational_modif_upper = f_test_modif$conf.int[2],
           n_mut_saltational = n_mut_saltational,
           n_mut_background = n_mut_background,
           n_other_muts_saltational = n_other_muts_saltational,
           n_other_muts_background = n_other_muts_background,
           p_value_salt = f_test$p.value)
  }))
  
  # Get OR of mutations in persistent infections
  df_OR_chronic_modif <- Reduce('bind_rows', lapply(1:nrow(df_chronic_muts_for_OR_comput), FUN = function(i_row){
    curr_gene <- df_chronic_muts_for_OR_comput$gene[i_row]
    curr_mut <- df_chronic_muts_for_OR_comput$mutation[i_row]
    
    tmp_df <- df_chronic_muts %>% 
      mutate(is_mut = (mutation == curr_mut) & (gene == curr_gene)) %>% 
      group_by(is_mut) %>% 
      summarise(n_chronic = sum(n_mut_chronic),
                n_high_quality_sequences = sum(n_mut_high_quality_sequence))
    
    n_mut_chronic <- tmp_df$n_chronic[tmp_df$is_mut]
    n_mut_background <- tmp_df$n_high_quality_sequences[tmp_df$is_mut]
    
    n_other_muts_chronic <- tmp_df$n_chronic[!tmp_df$is_mut]
    n_other_muts_background <- tmp_df$n_high_quality_sequences[!tmp_df$is_mut]
    
    mat_fisher <- matrix(c(n_mut_chronic, n_other_muts_chronic, n_mut_background, n_other_muts_background),
                         ncol = 2, nrow = 2)
    f_test <- fisher.test(x = mat_fisher)
    
    mat_fisher_modif <- matrix(c(n_mut_chronic + 1, n_other_muts_chronic + 1,
                                 n_mut_background + 1, n_other_muts_background + 1),
                               ncol = 2, nrow = 2)
    f_test_modif <- fisher.test(x = mat_fisher_modif)
    
    tibble(gene = curr_gene, 
           mutation = curr_mut,
           n_mut_chronic = n_mut_chronic, 
           n_mut_background = n_mut_background,
           n_other_muts_chronic = n_other_muts_chronic,
           n_other_muts_background = n_other_muts_background,
           OR_chronic = f_test$estimate,
           OR_chronic_lower = f_test$conf.int[1],
           OR_chronic_upper = f_test$conf.int[2],
           OR_chronic_modif = f_test_modif$estimate,
           OR_chronic_modif_lower = f_test_modif$conf.int[1],
           OR_chronic_modif_upper = f_test_modif$conf.int[2],
           p_value_chronic = f_test$p.value)
  }))
  
  df_OR_chronic_saltational <- df_chronic_muts %>% select(gene, mutation) %>% 
    left_join(df_OR_chronic_modif %>%
                rename(n_mut_high_quality_sequence = n_mut_background,
                       n_other_muts_high_quality_sequence = n_other_muts_background),
              by = c('gene', 'mutation')) %>% 
    left_join(df_modif_OR_saltational, by = c('gene', 'mutation')) %>% 
    filter(n_mut_saltational + n_mut_background >= 1, 
           n_mut_chronic + n_mut_high_quality_sequence >= 1)
  
  return(df_OR_chronic_saltational)
}


compute_OR_muts_chronic_infections_saltational_by_gene <- function(file_large_chronic_mut_dataset, file_aa_muts, name_saltational_branches, gene_of_interest){
  # Load chronic mutation counts
  df_chronic_muts <- get_chronic_muts_count(file_large_chronic_mut_dataset)
  
  # Get mutations in viridian tree
  df_aa_muts <- read_tsv(file_aa_muts) %>% 
    filter(is_coded) %>% # Coded AA mut
    filter(aa_new != aa_initial, gene == gene_of_interest) %>% # Non-synonymous
    mutate(full_mut = paste0(aa_muts_position, aa_new),
           is_saltational = node_id %in% name_saltational_branches) %>% 
    group_by(gene, full_mut, is_saltational) %>% 
    summarise(n_muts = n()) %>% 
    group_by(is_saltational) %>% 
    mutate(n_muts_tot = sum(n_muts)) %>% 
    ungroup()
  
  # Get name of mutations observed in viridian tree
  df_muts_in_persistent_inf <- df_aa_muts %>% 
    group_by(gene, full_mut) %>% 
    summarise(is_within_saltational = sum(is_saltational) > 0)
  
  # Get mutations observed in persistent infections dataset or in saltational branches
  # but removing the ones that are never observed in any of the datasets
  df_chronic_muts <- df_chronic_muts %>% 
    filter(gene == gene_of_interest) %>% 
    full_join(df_muts_in_persistent_inf, by = c('gene', 'mutation' = 'full_mut')) %>% 
    filter(! is.na(is_within_saltational)) %>%   # Remove mutations not observed in the viridian tree
    filter(! is.na(n_mut_chronic)) # Remove mutations not observed in the EPCI dataset
  
  df_chronic_muts_for_OR_comput <- df_chronic_muts %>% 
    filter(is_within_saltational | n_mut_chronic > 0)
  
  # Get OR of mut in saltational branches
  df_modif_OR_saltational <- Reduce('bind_rows', lapply(1:nrow(df_chronic_muts_for_OR_comput), FUN = function(i_row){
    curr_gene <- df_chronic_muts_for_OR_comput$gene[i_row]
    curr_mut <- df_chronic_muts_for_OR_comput$mutation[i_row]
    
    tmp_df <- df_aa_muts %>% 
      filter((full_mut == curr_mut) & (gene == curr_gene))
    
    n_mut_saltational <- tmp_df$n_muts[tmp_df$is_saltational]
    n_mut_saltational <- ifelse(length(n_mut_saltational) == 0, 0, n_mut_saltational)
    
    n_mut_background <- tmp_df$n_muts[(!tmp_df$is_saltational)]
    n_mut_background <- ifelse(length(n_mut_background) == 0, 0, n_mut_background)
    
    n_other_muts_saltational <- tmp_df$n_muts_tot[tmp_df$is_saltational] - n_mut_saltational
    n_other_muts_background <-  tmp_df$n_muts_tot[! tmp_df$is_saltational] - n_mut_background
    
    n_other_muts_saltational <- ifelse(length(n_other_muts_saltational) == 0, df_aa_muts$n_muts_tot[df_aa_muts$is_saltational][1], n_other_muts_saltational)
    n_other_muts_background <-  ifelse(length(n_other_muts_background) == 0, df_aa_muts$n_muts_tot[!df_aa_muts$is_saltational][1], n_other_muts_background)
    
    mat_fisher <- matrix(c(n_mut_saltational, n_other_muts_saltational, n_mut_background, n_other_muts_background),
                         ncol = 2, nrow = 2)
    f_test <- fisher.test(x = mat_fisher)
    
    mat_fisher_modif <- matrix(c(n_mut_saltational + 1, n_other_muts_saltational + 1, 
                                 n_mut_background + 1, n_other_muts_background + 1),
                               ncol = 2, nrow = 2)
    f_test_modif <- fisher.test(x = mat_fisher_modif)
    
    tibble(gene = curr_gene, 
           mutation = curr_mut,
           OR_saltational = f_test$estimate,
           OR_saltational_lower = f_test$conf.int[1],
           OR_saltational_upper = f_test$conf.int[2],
           OR_saltational_modif = f_test_modif$estimate,
           OR_saltational_modif_lower = f_test_modif$conf.int[1],
           OR_saltational_modif_upper = f_test_modif$conf.int[2],
           n_mut_saltational = n_mut_saltational,
           n_mut_background = n_mut_background,
           n_other_muts_saltational = n_other_muts_saltational,
           n_other_muts_background = n_other_muts_background,
           p_value_salt = f_test$p.value)
  }))
  
  # Get OR of mutations in persistent infections
  df_OR_chronic_modif <- Reduce('bind_rows', lapply(1:nrow(df_chronic_muts_for_OR_comput), FUN = function(i_row){
    curr_gene <- df_chronic_muts_for_OR_comput$gene[i_row]
    curr_mut <- df_chronic_muts_for_OR_comput$mutation[i_row]
    
    tmp_df <- df_chronic_muts %>% 
      mutate(is_mut = (mutation == curr_mut) & (gene == curr_gene)) %>% 
      group_by(is_mut) %>% 
      summarise(n_chronic = sum(n_mut_chronic),
                n_high_quality_sequences = sum(n_mut_high_quality_sequence))
    
    n_mut_chronic <- tmp_df$n_chronic[tmp_df$is_mut]
    n_mut_background <- tmp_df$n_high_quality_sequences[tmp_df$is_mut]
    
    n_other_muts_chronic <- tmp_df$n_chronic[!tmp_df$is_mut]
    n_other_muts_background <- tmp_df$n_high_quality_sequences[!tmp_df$is_mut]
    
    mat_fisher <- matrix(c(n_mut_chronic, n_other_muts_chronic, n_mut_background, n_other_muts_background),
                         ncol = 2, nrow = 2)
    f_test <- fisher.test(x = mat_fisher)
    
    mat_fisher_modif <- matrix(c(n_mut_chronic + 1, n_other_muts_chronic + 1,
                                 n_mut_background + 1, n_other_muts_background + 1),
                               ncol = 2, nrow = 2)
    f_test_modif <- fisher.test(x = mat_fisher_modif)
    
    tibble(gene = curr_gene, 
           mutation = curr_mut,
           n_mut_chronic = n_mut_chronic, 
           n_mut_background = n_mut_background,
           n_other_muts_chronic = n_other_muts_chronic,
           n_other_muts_background = n_other_muts_background,
           OR_chronic = f_test$estimate,
           OR_chronic_lower = f_test$conf.int[1],
           OR_chronic_upper = f_test$conf.int[2],
           OR_chronic_modif = f_test_modif$estimate,
           OR_chronic_modif_lower = f_test_modif$conf.int[1],
           OR_chronic_modif_upper = f_test_modif$conf.int[2],
           p_value_chronic = f_test$p.value)
  }))
  
  df_OR_chronic_saltational <- df_chronic_muts %>% select(gene, mutation) %>% 
    left_join(df_OR_chronic_modif %>%
                rename(n_mut_high_quality_sequence = n_mut_background,
                       n_other_muts_high_quality_sequence = n_other_muts_background),
              by = c('gene', 'mutation')) %>% 
    left_join(df_modif_OR_saltational, by = c('gene', 'mutation')) %>% 
    filter(n_mut_saltational + n_mut_background >= 1, 
           n_mut_chronic + n_mut_high_quality_sequence >= 1)
  
  return(df_OR_chronic_saltational)
}

get_prop_salt_branches_per_location <- function(name_saltational_branches, df_inferred_countries){
  df_prop_salt_per_location <- df_inferred_countries %>% 
    filter(strain != 'node_1') %>% 
    mutate(is_saltational = strain %in% name_saltational_branches) %>% 
    group_by(inferred_location, scale, is_saltational) %>% 
    summarise(n_branches = n()) %>% 
    group_by(inferred_location, scale) %>% 
    summarise(n_saltational_branches = ifelse(sum(is_saltational) == 0, 0, n_branches[is_saltational]),
              n_tot_branches = sum(n_branches)) %>% 
    mutate(prop_saltational_branches = n_saltational_branches/n_tot_branches)
 
  
  df_CI_prop_salt_branches <- binom.confint(
    x = df_prop_salt_per_location$n_saltational_branches,
    n = df_prop_salt_per_location$n_tot_branches,
    methods = "wilson"
  ) %>% 
    select(lower, upper) %>% 
    rename(lower_prop_salt = lower, upper_prop_salt = upper) 

  return(df_prop_salt_per_location %>% bind_cols(df_CI_prop_salt_branches))

}

get_prop_salt_branches_HIV_burden_per_location <- function(name_saltational_branches, 
                                                           df_inferred_countries,
                                                           threshold_n_branches_inclusion, 
                                                           file_path_country_demographics,
                                                           file_path_HIV_burden){
  # Get proportion of saltational branches by country
  df_prop_salt_per_location <- get_prop_salt_branches_per_location(name_saltational_branches, df_inferred_countries) %>% 
    filter(n_tot_branches >= threshold_n_branches_inclusion, inferred_location != 'Multiple descendant locations') %>% 
    mutate(inferred_location = ifelse(inferred_location == 'United Kingdom', 'UK', inferred_location))
  
  # Load country demographic data
  countries_demographics <- read_csv(file_path_country_demographics) %>% 
    select(country, population_size_2022_worldbank)
  countries_ihme_deaths <- read_csv(file_path_HIV_burden) %>% 
    filter(cause == 'HIV/AIDS') %>% 
    group_by(country) %>% summarise(n_deaths = sum(n_deaths)) %>% 
    rename(n_deaths_hiv = n_deaths) %>% 
    mutate(country = ifelse(country == 'United Kingdom', 'UK', country),
           country = ifelse(country == 'United States of America', 'USA', country))
  
  
  df_for_ecological_plot <- df_prop_salt_per_location %>% 
    left_join(countries_ihme_deaths, by = c('inferred_location' = 'country')) %>% 
    left_join(countries_demographics, by = c('inferred_location' = 'country')) %>% 
    group_by(inferred_location, scale) %>% 
    mutate(hiv_deaths_per_1000 = n_deaths_hiv / population_size_2022_worldbank * 1000,
           modif_saltational_branches = (n_saltational_branches + 1) / (n_tot_branches + 1))
  
  return(df_for_ecological_plot)
 
}
# Functions to make plots
plot_OR_func_nb_mut_colour_recurrent_mut <- function(df_OR_AA_mut){
  
  label_125 <- function(x) {
    m <- signif(x / 10^floor(log10(x)), 12) 
    labs <- x
    ifelse(m %in% c(1, 2, 5), labs, "")
  }
  
  plt_OR_func_mut_frequency_coloured_by_profile_type <- df_OR_AA_mut %>% 
    ggplot(aes(x = n_mut_of_interest_adapt_branches, y = OR, 
               colour = type, alpha = is_signif)) +
    geom_point() +
    geom_text_repel(data = df_OR_AA_mut %>% 
                      filter(is_signif, type != 'Others'),
                    aes(label = aa_muts_tot),
                    size = 3.5,
                    max.time = 20, min.segment.length = 0., max.overlaps = Inf) +
    scale_x_continuous(labels = label_125,
                       trans = 'log',
                       expand = expansion(mult = c(0.065, 0.06)),
                       breaks = c(seq(1, 9, 1), seq(1e1, 9e1, 1e1), seq(1e2, 9e2, 1e2)),
                       name = 'Number of saltational branches with mutation')  +
    scale_y_continuous(name = 'OR',
                       labels = label_125,
                       trans = 'log', 
                       expand = expansion(mult = c(0.05, 0.02)),
                       breaks = c(seq(1e-1, 9e-1, 1e-1), seq(1, 9, 1), seq(10, 90, 10),
                                  seq(100, 900, 100), seq(1000, 9000, 1000))
    ) +
    scale_colour_manual(breaks = c('Persistent infections-like', 'Deer-like', 'Mink-like'), 
                        values = c(
                          'Persistent infections-like' = col_persistent,
                          'Deer-like' = col_deers,
                          'Mink-like' = col_minks,
                          'Others' = col_others
                        ),
                        labels = c('Persistent infections', 'Deers', 'Minks'),
                        name = 'Recurrent mutations in') +
    scale_alpha_manual(breaks = c(T, F), values = c(0.9, .01), guide = 'none') +
    theme_classic() +
    theme(axis.text = element_text(size = 11),
          legend.background = element_blank(),
          legend.position = 'inside',
          legend.position.inside = c(0.75, 0.98),
          legend.key.size = unit(0.1, "cm")) +
    coord_cartesian(ylim = c(1, NA)) +
    guides(colour = guide_legend(override.aes = list(label = '')))
  
  return(plt_OR_func_mut_frequency_coloured_by_profile_type)
}

plot_OR_chronic_func_OR_saltational <- function(df_OR_chronic_saltational){
  label_125 <- function(x) {
    x <- ifelse(abs(x - round(x)) < 1e-5, round(x), x)
    m <- as.integer(signif(x / 10^floor(log10(x) + 1e-8), 12))
    ifelse(m %in% as.integer(c(1, 2, 5)), x, "")
  }
  label_1 <- function(x) {
    x <- ifelse(abs(x - round(x)) < 1e-5, round(x), x)
    m <- as.integer(signif(x / 10^floor(log10(x) + 1e-8), 12))
    ifelse(m %in% as.integer(c(1)), x, "")
  }
  
  plt_comp_OR_saltational_OR_chronic <- df_OR_chronic_saltational %>% 
    mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
    arrange(n_mut_high_quality_sequence) %>% 
    ggplot(aes(x = OR_chronic_modif, y = OR_saltational_modif)) +
    geom_point(aes(colour = gene), alpha = 0.8) +
    stat_cor(method = 'pearson', cor.coef.name = 'r') +
    scale_x_continuous(trans = 'log', name = expression(OR['chronic']),
                     labels = label_1,
                     breaks = c(seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 9e-1, 1e-1),
                                seq(1e0, 9e0, 1e0), seq(1e1, 9e1, 1e1), seq(1e2, 9e2, 1e2), 
                                seq(1e3, 9e3, 1e3))) +
    scale_y_continuous(trans = 'log', name = expression(OR['saltational']),
                       labels = label_125,
                       breaks = c(seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 9e-1, 1e-1),
                                  seq(1e0, 9e0, 1e0), seq(1e1, 9e1, 1e1), seq(1e2, 9e2, 1e2))) +
    scale_colour_manual(values = palette_genes, name = 'Gene') +
    theme_classic() +
    coord_fixed() +
    theme(axis.text = element_text(size = 11))
  
  return(plt_comp_OR_saltational_OR_chronic)
}

plot_OR_chronic_func_OR_saltational_cor_by_gene <- function(df_OR_chronic_saltational){
  label_125 <- function(x) {
    x <- ifelse(abs(x - round(x)) < 1e-5, round(x), x)
    m <- as.integer(signif(x / 10^floor(log10(x) + 1e-8), 12))
    ifelse(m %in% as.integer(c(1, 2, 5)), x, "")
  }
  label_1 <- function(x) {
    x <- ifelse(abs(x - round(x)) < 1e-5, round(x), x)
    m <- as.integer(signif(x / 10^floor(log10(x) + 1e-8), 12))
    ifelse(m %in% as.integer(c(1)), x, "")
  }
  
  plt_comp_OR_saltational_OR_chronic <- df_OR_chronic_saltational %>% 
    mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
    arrange(n_mut_high_quality_sequence) %>% 
    ggplot(aes(x = OR_chronic_modif, y = OR_saltational_modif, colour = gene)) +
    geom_point(alpha = 0.8) +
    stat_cor(method = 'pearson') +
    scale_x_continuous(trans = 'log', name = expression(OR['chronic']),
                       labels = label_1,
                       breaks = c(seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 9e-1, 1e-1),
                                  seq(1e0, 9e0, 1e0), seq(1e1, 9e1, 1e1), seq(1e2, 9e2, 1e2), 
                                  seq(1e3, 9e3, 1e3))) +
    scale_y_continuous(trans = 'log', name = expression(OR['saltational']),
                       labels = label_125,
                       breaks = c(seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 9e-1, 1e-1),
                                  seq(1e0, 9e0, 1e0), seq(1e1, 9e1, 1e1), seq(1e2, 9e2, 1e2))) +
    scale_colour_manual(values = palette_genes, name = 'Gene') +
    theme_classic() +
    coord_fixed() +
    theme(axis.text = element_text(size = 11))
  
  return(plt_comp_OR_saltational_OR_chronic)
}

plot_prop_salt_branches_func_HIV_burden <- function(df_prop_salt_branches_HIV_burden_per_country){
  label_125 <- function(x) {
    m <- signif(x / 10^floor(log10(x)), 12) 
    labs <- x
    ifelse(m %in% c(1, 2, 5), labs, "")
  }
  label_125_percent <- function(x) {
    m <- signif(x / 10^floor(log10(x)), 12)   # mantissa
    labs <- scales::label_percent(accuracy = 1e-2)(x)
    ifelse(m %in% c(1, 2, 5), labs, "")
  }
  
  df_country_mapping <- read_csv('../../params/mapping/country_to_continent_mapping.csv') %>% 
    select(country, continent, subregion) %>% distinct() %>% 
    mutate(country = ifelse(country == 'United Kingdom', 'UK', country))
  
  plt_prop_func_HIV_burden_by_country <- df_prop_salt_branches_HIV_burden_per_country %>% 
      left_join(df_country_mapping, by = c('inferred_location' = 'country')) %>% 
      ggplot(aes(x = hiv_deaths_per_1000, 
                 y = modif_saltational_branches)) +
      geom_point(aes(colour = continent)) +
      geom_text_repel(aes(label = inferred_location, colour = continent), max.overlaps = Inf) +
      scale_colour_manual(breaks = c('North America', 'South America', 'Europe', 'Africa', 'Asia', 'Oceania'),
                          values = c(palette_north_america[10], palette_south_america[10], 
                                     palette_europe[10], palette_africa[10], 
                                     palette_asia[10], palette_oceania[10]),
                          guide = 'none') +
      scale_x_continuous(trans = 'log', name = 'Yearly HIV deaths (per 1,000)',
                         breaks = c(seq(1e-4, 9e-4, 1e-4), seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 9e-1, 1e-1),
                                    seq(1e0, 9e0, 1e0)),
                         labels = c(expression(10^{-4}), rep('', 8),
                                    expression(10^{-3}), rep('', 8),
                                    expression(10^{-2}), rep('', 8),
                                    expression(10^{-1}), rep('', 8),
                                    expression(10^{0}), rep('', 8))) + 
      scale_y_continuous(trans = 'log', name = 'Proportion of branches that\nare saltational (adjusted)',
                         breaks = c(seq(1e-4, 9e-4, 1e-4), seq(1e-3, 9e-3, 1e-3)),
                         labels = label_125_percent) +
      stat_cor(method = 'spearman', label.x = -2.5, label.y = -9., cor.coef.name = 'rho', 
               aes(x = hiv_deaths_per_1000, y = prop_saltational_branches + 1e-20)) + # Add pseudo-count to ensure Australia is not dropped from computation
      theme_classic() +
      theme(panel.grid.minor = element_blank(),
            axis.text = element_text(size = 11)) 
 
  return(plt_prop_func_HIV_burden_by_country)   
}
