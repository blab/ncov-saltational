library(extraDistr)
library(tidyverse)
library(argparse)
library(zipfR)

source('scripts/bayesian-model/utils_model.R')

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--data_viridian_grouped', type = 'character', help = 'TSV file with data used to run model')
  parser$add_argument('--df_prior_param', type = 'character', help = 'TSV file with prior parametrization')
  parser$add_argument('--CV_prior', type = 'double', help = 'CV for prior parametrization')
  parser$add_argument('--genome_region', type = 'character', help = 'Genome region')
  parser$add_argument('--output_param_estim', type = 'character', help = 'TSV file for model run')
  return(parser$parse_args())
}

args <- collect_args()

# Load prior parametrisation
df_prior <- read_tsv(args$df_prior_param)

# Load gene data
gene_data <- read_tsv(args$data_viridian_grouped) %>% 
  mutate(i_row = 1:n())

# Function to run model
run_model <- function(df_prior, data_gene, gene, CV_prior){
  
  # Define prior parameters
  prior_param <- list('alpha_N' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dN[df_prior$gene == gene]), cv = CV_prior)['alpha'], 
                      'beta_N' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dN[df_prior$gene == gene]), cv = CV_prior)['beta'],
                      'alpha_S' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dS[df_prior$gene == gene]), cv = CV_prior)['alpha'], 
                      'beta_S' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dS[df_prior$gene == gene]), cv = CV_prior)['beta'])
 
  df_param_estim <- Reduce('bind_rows', lapply(1:nrow(gene_data), FUN = function(i_row){
    
    data_gene <- list('N' = gene_data$non_synonymous_muts[i_row], 
                      'S' = gene_data$synonymous_muts[i_row], 
                      'LN' = gene_data$nonsyn_sites[i_row],
                      'LS' = gene_data$syn_sites[i_row])
    
    mod <- create_model(data_gene, prior_param)
    
    df_dN <- tibble(posterior_mean = mod$mean_posterior_dN,
                   posterior_mode = mod$mode_posterior_dN,
                   posterior_median = mod$get_quantile_posterior_dN(0.5),
                   posterior_2.5 = mod$get_quantile_posterior_dN(0.025),
                   posterior_97.5 = mod$get_quantile_posterior_dN(0.975)) %>% 
      mutate(param = 'dN',
             i_row = i_row)
    
    df_dS <- tibble(posterior_mean = mod$mean_posterior_dS,
                    posterior_mode = mod$mode_posterior_dS,
                    posterior_median = mod$get_quantile_posterior_dS(0.5),
                    posterior_2.5 = mod$get_quantile_posterior_dS(0.025),
                    posterior_97.5 = mod$get_quantile_posterior_dS(0.975)) %>% 
      mutate(param = 'dS',
             i_row = i_row)
    
    df_omega <- tibble(posterior_mean = mod$mean_posterior_omega,
                       posterior_mode = mod$mode_posterior_omega,
                       posterior_median = mod$get_quantile_posterior_omega(0.5),
                       posterior_2.5 = mod$get_quantile_posterior_omega(0.025),
                       posterior_97.5 = mod$get_quantile_posterior_omega(0.975),
                       posterior_proba_omega_greater_1 = mod$posterior_proba_omega_greater_1,
                       prior_proba_omega_greater_1 = mod$prior_proba_omega_greater_1) %>% 
      mutate(param = 'omega',
             i_row = i_row)
    
    bind_rows(df_dN, df_dS, df_omega) %>% 
      mutate(marginal_LH_N = mod$marginal_likelihood_N,
             marginal_LH_S = mod$marginal_likelihood_S)
  })) %>% 
    left_join(gene_data) %>% 
    mutate(gene = gene, CV_prior = CV_prior)
  
  return(df_param_estim)
}

# Model run
df_param_estim <- run_model(df_prior, gene_data, args$genome_region, args$CV_prior)

# Save model run
write_tsv(df_param_estim, args$output_param_estim)
