library(tidyverse)
library(scales)
library(glmnet)
source('utils_descendants_analysis.R')
source('config_plot.R')

# Load saltational branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.7_BF_20.tsv', header = F) %>%
  unlist() %>% as.character()

# Load descendants information
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(has_descendants = n_terminal_descendants > 1) %>% 
  select(strain, has_descendants)

# Load AA mutation information
df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  rename(strain = node_id) %>% 
  filter(is_coded, strain %in% name_saltational_branches) %>% 
  mutate(is_synonymous = (aa_initial == aa_new),
         aa_muts_direction = paste0(aa_muts_position, aa_new))

get_OR_having_descendants_having_muts <- function(df_aa_muts, 
                                                  df_descendants,
                                                  threshold_mutation_occurence_inclusion = 4){
  
  df_muts_to_include <- df_aa_muts %>% 
    group_by(gene, aa_muts_direction, is_synonymous) %>% 
    summarise(n_occ = n()) %>% 
    filter(n_occ >= threshold_mutation_occurence_inclusion)
  
  df_OR_having_descendants_by_mut <- Reduce('bind_rows', lapply(1:nrow(df_muts_to_include), FUN = function(i_row){
    tmp_df_has_mut <- df_aa_muts %>% 
      filter(gene == df_muts_to_include$gene[i_row],
             aa_muts_direction == df_muts_to_include$aa_muts_direction[i_row],
             is_synonymous == df_muts_to_include$is_synonymous[i_row]) %>% 
      mutate(has_mut = T) %>% 
      select(strain, has_mut)
    
    tmp_df_descendants_muts <- df_descendants %>% 
      left_join(tmp_df_has_mut, by = 'strain') %>% 
      mutate(has_mut = replace_na(has_mut, F))
    
    n_with_mut_with_descendants <- sum((tmp_df_descendants_muts$has_descendants) & (tmp_df_descendants_muts$has_mut))
    n_with_mut_without_descendants <- sum((!tmp_df_descendants_muts$has_descendants) & (tmp_df_descendants_muts$has_mut))
    n_without_mut_with_descendants <- sum((tmp_df_descendants_muts$has_descendants) & (!tmp_df_descendants_muts$has_mut))
    n_without_mut_without_descendants <- sum((!tmp_df_descendants_muts$has_descendants) & (!tmp_df_descendants_muts$has_mut))
    
    mat_fisher <- matrix(c(n_with_mut_with_descendants, n_with_mut_without_descendants, n_without_mut_with_descendants, n_without_mut_without_descendants),
                         nrow = 2, ncol = 2)
    f_test <- fisher.test(x = mat_fisher)
    
    tibble(gene = df_muts_to_include$gene[i_row],
           aa_muts_direction = df_muts_to_include$aa_muts_direction[i_row],
           is_synonymous = df_muts_to_include$is_synonymous[i_row],
           OR = as.numeric(f_test$estimate),
           OR_lower = as.numeric(f_test$conf.int[1]),
           OR_upper = as.numeric(f_test$conf.int[2]),
           p_value = as.numeric(f_test$p.value))
  }))
  
  return(df_OR_having_descendants_by_mut)
}


df_OR_having_descendants_by_mut <- 
  get_OR_having_descendants_having_muts(df_aa_muts, 
                                        df_descendants,
                                        threshold_mutation_occurence_inclusion = 4)

df_OR_having_descendants_by_mut %>% 
  ggplot(aes(x = log10(OR), y = -log10(p_value), colour = gene, shape = is_synonymous)) +
  geom_point() +
  ggrepel::geom_text_repel(data = df_OR_having_descendants_by_mut %>% filter(abs(log10(OR)) > 0.5),
                           aes(label = aa_muts_direction), 
                           max.overlaps = Inf) +
  facet_wrap(. ~ gene) +
  scale_colour_manual(values = palette_genes)



####
run_lasso_log_reg <- function(df_aa_muts, 
                              df_descendants,
                              threshold_mutation_occurence_inclusion = 4){
  
  # Create feature matrix only retaining the most frequent mutations
  df_muts_to_include <- df_aa_muts %>% 
    group_by(gene, aa_muts_direction, is_synonymous) %>% 
    summarise(n_occ = n()) %>% 
    filter(n_occ > threshold_mutation_occurence_inclusion) %>% 
    mutate(feature = paste(gene, aa_muts_direction, is_synonymous, sep = "__"))
  
  df_X <- df_aa_muts %>% 
    inner_join(
      df_muts_to_include %>% select(gene, aa_muts_direction, is_synonymous, feature),
      by = c("gene", "aa_muts_direction", "is_synonymous")
    ) %>% 
    select(strain, feature) %>% 
    distinct() %>% 
    mutate(value = 1L) %>% 
    tidyr::pivot_wider(
      names_from = feature,
      values_from = value,
      values_fill = 0
    )
  
  # Define model dataframe
  df_model <- df_descendants %>% 
    select(strain, has_descendants) %>% 
    left_join(df_X, by = "strain")
  
  mutation_cols <- setdiff(colnames(df_model), c("strain", "has_descendants"))
  
  df_model <- df_model %>% 
    mutate(across(all_of(mutation_cols), ~ replace_na(., 0L)))
  
  # Run LASSO
  X_matrix <- df_model %>% 
    select(all_of(mutation_cols)) %>% 
    as.matrix()
  
  y_outcome <- df_model %>% 
    pull(has_descendants) %>% 
    as.integer()
  
  cross_validation_lasso <- cv.glmnet(x = X_matrix, y = y_outcome,
                                      family = 'binomial', alpha = 1.)
  coef_lasso_lambda_min <- coef(cross_validation_lasso, s = cross_validation_lasso$lambda.min)
  
  cross_validation_elastic_net <- cv.glmnet(x = X_matrix, y = y_outcome,
                                            family = 'binomial', alpha = 0.5)
  coef_elastic_net_lambda_min <- coef(cross_validation_elastic_net, 
                                      s = cross_validation_elastic_net$lambda.min)
  
  cross_validation_ridge <- cv.glmnet(x = X_matrix, y = y_outcome,
                                            family = 'binomial', alpha = 0.)
  coef_ridge_lambda_min <- coef(cross_validation_ridge, 
                                s = cross_validation_ridge$lambda.min)
  
  df_ridge <- tibble(param = rownames(coef_ridge_lambda_min), 
                     coef = coef_ridge_lambda_min[,1]) %>% 
    filter(param != '(Intercept)') %>% 
    tidyr::separate(
      param,
      into = c("gene", "aa_muts_direction", "is_synonymous"),
      sep = "__",
      remove = FALSE
    )
}

df_ridge %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  ggplot(aes(x = gene, y = coef)) +
  geom_jitter() +
  facet_wrap(. ~ is_synonymous)
