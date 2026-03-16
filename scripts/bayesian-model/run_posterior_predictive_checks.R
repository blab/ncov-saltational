library(extraDistr)
library(tidyverse)
library(argparse)
library(zipfR)
library(doParallel)
library(foreach)

source('scripts/bayesian-model/utils_model.R')

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--df_prior_param', type = 'character', help = 'TSV file with prior parametrization')
  parser$add_argument('--genome_region', type = 'character', help = 'Genome region')
  parser$add_argument('--CV_prior', type = 'double', help = 'CV for prior parametrization')
  parser$add_argument('--data_viridian_grouped', type = 'character', help = 'TSV file with data used to run model')
  parser$add_argument('--output_posterior_simul', type = 'character', help = 'TSV file for posterior simul')
  return(parser$parse_args())
}

args <- collect_args()

# Load info for prior parametrisation
df_prior <- read_tsv(args$df_prior_param)

# Load gene data
gene_data <- read_tsv(args$data_viridian_grouped) %>% 
  mutate(i_row = 1:n())

# Function to run posterior simulation
simulate_from_posterior <- function(df_prior, gene_data, gene, CV_prior){
  
  # Define prior parameters
  prior_param <- list('alpha_N' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dN[df_prior$gene == gene]), cv = CV_prior)['alpha'], 
                      'beta_N' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dN[df_prior$gene == gene]), cv = CV_prior)['beta'],
                      'alpha_S' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dS[df_prior$gene == gene]), cv = CV_prior)['alpha'], 
                      'beta_S' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dS[df_prior$gene == gene]), cv = CV_prior)['beta'])
  
  # Load gene data
  n_replicate_simul <- 100
  
  n_cores <- 3
  cl <- makeForkCluster(n_cores)
  registerDoParallel(cl)
  t0 <- Sys.time()
  df_simul_from_posterior <- Reduce('bind_rows', foreach(i_row = 1:nrow(gene_data), .packages = c('tidyverse')) %dopar% {
    data_gene <- list('N' = gene_data$non_synonymous_muts[i_row], 
                      'S' = gene_data$synonymous_muts[i_row], 
                      'LN' = gene_data$nonsyn_sites[i_row],
                      'LS' = gene_data$syn_sites[i_row])
    
    mod <- create_model(data_gene, prior_param)
    
    Reduce('bind_rows', lapply(1:n_replicate_simul, FUN = function(i_simul){
      tibble(id_seq = 1:gene_data$n_branches[i_row],
             S_new = mod$simulate_from_posterior_S(n_rep = gene_data$n_branches[i_row]),
             N_new = mod$simulate_from_posterior_N(n_rep = gene_data$n_branches[i_row])) %>% 
        mutate(LN = gene_data$nonsyn_sites[i_row],
               LS = gene_data$syn_sites[i_row]) %>% 
        group_by(S_new, N_new, LN, LS) %>% 
        summarise(n_branches_new = n()) %>% 
        mutate(i_row = i_row, i_simul = i_simul)
    }))
  }) %>% 
    mutate(gene = gene, CV_prior = CV_prior)
  
  return(df_simul_from_posterior)
}

df_simul_from_posterior <- simulate_from_posterior(df_prior, gene_data, args$genome_region, args$CV_prior)
write_tsv(df_simul_from_posterior, args$output_posterior_simul)
