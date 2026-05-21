get_mut_counts_per_gene_adapt_branches <- function(file_path_aa_muts, name_adaptive_branches){
  df_aa_muts_count_adaptive_branches <- read_tsv(file_path_aa_muts) %>% 
    filter(is_coded, node_id %in% name_adaptive_branches) %>% 
    rename(strain = node_id) %>% 
    group_by(strain) %>% 
    summarise(n_muts = n(), n_non_syn = sum(aa_initial != aa_new),
              n_muts_S = sum(gene == 'S'), n_non_syn_S = sum((aa_initial != aa_new) & (gene == 'S')),
              n_muts_ORF1a = sum(gene == 'ORF1a'), n_non_syn_ORF1a = sum((aa_initial != aa_new) & (gene == 'ORF1a')),
              n_muts_ORF1b = sum(gene == 'ORF1b'), n_non_syn_ORF1b = sum((aa_initial != aa_new) & (gene == 'ORF1b')),
              n_muts_ORF3a = sum(gene == 'ORF3a'), n_non_syn_ORF3a = sum((aa_initial != aa_new) & (gene == 'ORF3a')),
              n_muts_ORF6 = sum(gene == 'ORF6'), n_non_syn_ORF6 = sum((aa_initial != aa_new) & (gene == 'ORF6')),
              n_muts_ORF7a = sum(gene == 'ORF7a'), n_non_syn_ORF7a = sum((aa_initial != aa_new) & (gene == 'ORF7a')),
              n_muts_ORF7b = sum(gene == 'ORF7b'), n_non_syn_ORF7b = sum((aa_initial != aa_new) & (gene == 'ORF7b')),
              n_muts_ORF8 = sum(gene == 'ORF8'), n_non_syn_ORF8 = sum((aa_initial != aa_new) & (gene == 'ORF8')),
              n_muts_ORF9b = sum(gene == 'ORF9b'), n_non_syn_ORF9b = sum((aa_initial != aa_new) & (gene == 'ORF9b')),
              n_muts_M = sum(gene == 'M'), n_non_syn_M = sum((aa_initial != aa_new) & (gene == 'M')),
              n_muts_N = sum(gene == 'N'), n_non_syn_N = sum((aa_initial != aa_new) & (gene == 'N')),
              n_muts_E = sum(gene == 'E'), n_non_syn_E = sum((aa_initial != aa_new) & (gene == 'E')),
              n_syn_S = n_muts_S - n_non_syn_S,
              n_syn_M = n_muts_M - n_non_syn_M,
              n_syn_N = n_muts_N - n_non_syn_N,
              n_syn_E = n_muts_E - n_non_syn_E,
              n_syn_ORF1a = n_muts_ORF1a - n_non_syn_ORF1a,
              n_syn_ORF1b = n_muts_ORF1b - n_non_syn_ORF1b,
              n_syn_ORF3a = n_muts_ORF3a - n_non_syn_ORF3a,
              n_syn_ORF6 = n_muts_ORF6 - n_non_syn_ORF6,
              n_syn_ORF7a = n_muts_ORF7a - n_non_syn_ORF7a,
              n_syn_ORF7b = n_muts_ORF7b - n_non_syn_ORF7b,
              n_syn_ORF8 = n_muts_ORF8 - n_non_syn_ORF8,
              n_syn_ORF9b = n_muts_ORF9b - n_non_syn_ORF9b)
  
  return(df_aa_muts_count_adaptive_branches)
  
}

get_AIC_poisson_reg <-  function(df_saltational_descendants, 
                                 name_cols,
                                 threshold_n_descendants_inclusion = Inf){
  # Remove outliers and select predictor columns
  df_for_reg <- df_saltational_descendants %>% filter(n_terminal_descendants_counts < threshold_n_descendants_inclusion) %>% 
    dplyr::select(n_terminal_descendants_counts, all_of(name_cols))
  
  # Run model
  poisson_model <-  glmmTMB(data = df_for_reg, 
                            n_terminal_descendants_counts ~ 
                              n_non_syn_ORF1a + n_non_syn_ORF1b + n_non_syn_S + 
                              n_non_syn_E + n_non_syn_M +
                              n_non_syn_ORF3a + n_non_syn_ORF6 +
                              n_non_syn_ORF7a + n_non_syn_ORF7b + n_non_syn_ORF8 + n_non_syn_ORF9b + n_non_syn_N +
                              n_syn_ORF1a + n_syn_ORF1b + n_syn_S + n_syn_ORF3a + n_syn_ORF6 +
                              n_syn_ORF7a + n_syn_ORF7b + n_syn_ORF8 + n_syn_ORF9b + n_syn_N +
                              n_syn_E + n_syn_M, 
                            ziformula = ~ 0, family = poisson)
  
  return(AIC(poisson_model))
}

get_AIC_negbin_reg <-  function(df_saltational_descendants, 
                                name_cols,
                                threshold_n_descendants_inclusion = Inf){
  # Remove outliers and select predictor columns
  df_for_reg <- df_saltational_descendants %>% filter(n_terminal_descendants_counts < threshold_n_descendants_inclusion) %>% 
    dplyr::select(n_terminal_descendants_counts, all_of(name_cols))
  
  # Run model
  negbin_model <-  glmmTMB(data = df_for_reg, 
                           n_terminal_descendants_counts ~ 
                             n_non_syn_ORF1a + n_non_syn_ORF1b + n_non_syn_S + 
                             n_non_syn_E + n_non_syn_M +
                             n_non_syn_ORF3a + n_non_syn_ORF6 +
                             n_non_syn_ORF7a + n_non_syn_ORF7b + n_non_syn_ORF8 + n_non_syn_ORF9b + n_non_syn_N +
                             n_syn_ORF1a + n_syn_ORF1b + n_syn_S + n_syn_ORF3a + n_syn_ORF6 +
                             n_syn_ORF7a + n_syn_ORF7b + n_syn_ORF8 + n_syn_ORF9b + n_syn_N +
                             n_syn_E + n_syn_M, 
                           ziformula = ~ 0, family = nbinom2)
  
  return(AIC(negbin_model))
}

likelihood_ratio_test_poisson_negbin <-  function(df_saltational_descendants, 
                                                  name_cols,
                                                  threshold_n_descendants_inclusion = Inf){
  # Remove outliers and select predictor columns
  df_for_reg <- df_saltational_descendants %>% filter(n_terminal_descendants_counts < threshold_n_descendants_inclusion) %>% 
    dplyr::select(n_terminal_descendants_counts, all_of(name_cols))
  
  # Run Poisson model
  poisson_model <-  glmmTMB(data = df_for_reg, 
                            n_terminal_descendants_counts ~ 
                              n_non_syn_ORF1a + n_non_syn_ORF1b + n_non_syn_S + 
                              n_non_syn_E + n_non_syn_M +
                              n_non_syn_ORF3a + n_non_syn_ORF6 +
                              n_non_syn_ORF7a + n_non_syn_ORF7b + n_non_syn_ORF8 + n_non_syn_ORF9b + n_non_syn_N +
                              n_syn_ORF1a + n_syn_ORF1b + n_syn_S + n_syn_ORF3a + n_syn_ORF6 +
                              n_syn_ORF7a + n_syn_ORF7b + n_syn_ORF8 + n_syn_ORF9b + n_syn_N +
                              n_syn_E + n_syn_M, 
                            ziformula = ~ 0, family = poisson)
  
  # Run Negative Binomial model
  negbin_model <-  glmmTMB(data = df_for_reg, 
                           n_terminal_descendants_counts ~ 
                             n_non_syn_ORF1a + n_non_syn_ORF1b + n_non_syn_S + n_non_syn_ORF3a + n_non_syn_ORF6 +
                             n_non_syn_ORF7a + n_non_syn_ORF7b + n_non_syn_ORF8 + n_non_syn_ORF9b + n_non_syn_N +
                             n_syn_ORF1a + n_syn_ORF1b + n_syn_S + n_syn_ORF3a + n_syn_ORF6 +
                             n_syn_ORF7a + n_syn_ORF7b + n_syn_ORF8 + n_syn_ORF9b + n_syn_N, 
                           ziformula = ~ 0, family = nbinom2)
  
  anova_poisson_negbin <- anova(poisson_model, negbin_model)
  return(anova_poisson_negbin)
}


run_negbin_reg <- function(df_saltational_descendants, 
                           name_cols,
                           threshold_n_descendants_inclusion = Inf){
  
  # Remove outliers and select predictor columns
  df_for_reg <- df_saltational_descendants %>% filter(n_terminal_descendants_counts < threshold_n_descendants_inclusion) %>% 
    dplyr::select(n_terminal_descendants_counts, all_of(name_cols))
  
  # Run model
  negbin_model <-  glmmTMB(data = df_for_reg, 
                           n_terminal_descendants_counts ~ 
                             n_non_syn_ORF1a + n_non_syn_ORF1b + n_non_syn_S + 
                             n_non_syn_E + n_non_syn_M +
                             n_non_syn_ORF3a + n_non_syn_ORF6 +
                             n_non_syn_ORF7a + n_non_syn_ORF7b + n_non_syn_ORF8 + n_non_syn_ORF9b + n_non_syn_N +
                             n_syn_ORF1a + n_syn_ORF1b + n_syn_S + n_syn_ORF3a + n_syn_ORF6 +
                             n_syn_ORF7a + n_syn_ORF7b + n_syn_ORF8 + n_syn_ORF9b + n_syn_N +
                             n_syn_E + n_syn_M, 
                           ziformula = ~ 0, family = nbinom2)
  # Get model parameters
  df_param_estim <- confint(negbin_model)[- 1, ] %>% 
    as_tibble() %>% 
    rename(param_estim = Estimate) %>% 
    rename(lower_param = `2.5 %`, upper_param = `97.5 %`) %>% 
    mutate(name_col = rownames(confint(negbin_model)[- 1, ])) %>% 
    group_by(name_col) %>% 
    mutate(gene = strsplit(name_col, '_')[[1]][length(strsplit(name_col, '_')[[1]])],
           mut_type = ifelse(startsWith(name_col, 'n_syn'), 'Synonymous', 'Nonsynonymous')) %>%
    ungroup() %>% 
    dplyr::select(- name_col)
  
  return(df_param_estim) 
}

