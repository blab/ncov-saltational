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
  parser$add_argument('--genome_region', type = 'character', help = 'Genome region (should be "Genome")')
  parser$add_argument('--output_param_estim', type = 'character', help = 'TSV file for model run')
  return(parser$parse_args())
}

args <- collect_args()

# Run hierarchical Bayesian model for genome-wide analysis
# Uses gene-level posteriors combined with site-weighted averaging
run_genome_model <- function(df_prior, genome_data_grouped, CV_prior){

  # Define gene vector
  vec_genes <- c('ORF1a', 'ORF1b', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF9b')

  # Create gene-specific prior parameters list
  prior_params_list <- lapply(vec_genes, FUN = function(curr_gene){
    list(
      'alpha_N' = get_gamma_param_from_mean_CV(
        mu = exp(df_prior$log_mean_dN[df_prior$gene == curr_gene]),
        cv = CV_prior)['alpha'],
      'beta_N' = get_gamma_param_from_mean_CV(
        mu = exp(df_prior$log_mean_dN[df_prior$gene == curr_gene]),
        cv = CV_prior)['beta'],
      'alpha_S' = get_gamma_param_from_mean_CV(
        mu = exp(df_prior$log_mean_dS[df_prior$gene == curr_gene]),
        cv = CV_prior)['alpha'],
      'beta_S' = get_gamma_param_from_mean_CV(
        mu = exp(df_prior$log_mean_dS[df_prior$gene == curr_gene]),
        cv = CV_prior)['beta']
    )
  })
  names(prior_params_list) <- vec_genes

  # Helper function to summarize Monte Carlo draws
  summarize_draws <- function(x) {
    c(
      posterior_mean = mean(x),
      posterior_mode = NA_real_,  # No closed form for hierarchical model
      posterior_median = median(x),
      posterior_2.5 = unname(quantile(x, 0.025)),
      posterior_97.5 = unname(quantile(x, 0.975))
    )
  }

  # Run model for each unique branch pattern
  df_param_estim <- Reduce('bind_rows', lapply(1:nrow(genome_data_grouped), FUN = function(i_row){

    branch_row <- genome_data_grouped[i_row, ]

    # Generate posterior and prior draws using hierarchical approach
    # (compute_genome_draws imported from utils_model.R)
    posterior_draws <- compute_genome_draws(
      branch_row, vec_genes, prior_params_list,
      n_rep = 1000, draw_type = 'posterior'
    )

    prior_draws <- compute_genome_draws(
      branch_row, vec_genes, prior_params_list,
      n_rep = 1000, draw_type = 'prior'
    )

    # Aggregate observed quantities across genes for output metadata
    synonymous_muts <- sum(sapply(vec_genes, function(g) branch_row[[paste0("synonymous_muts_", g)]]))
    non_synonymous_muts <- sum(sapply(vec_genes, function(g) branch_row[[paste0("non_synonymous_muts_", g)]]))
    syn_sites <- sum(sapply(vec_genes, function(g) branch_row[[paste0("syn_sites_", g)]]))
    nonsyn_sites <- sum(sapply(vec_genes, function(g) branch_row[[paste0("nonsyn_sites_", g)]]))

    # Summarize posterior draws for each parameter
    summary_posterior_dN <- summarize_draws(posterior_draws$dN_gen)
    summary_posterior_dS <- summarize_draws(posterior_draws$dS_gen)
    summary_posterior_omega <- summarize_draws(posterior_draws$omega_gen)

    # Compute probability that omega > 1 (positive selection)
    posterior_proba_omega_above_1 <- sum(posterior_draws$omega_gen > 1) / length(posterior_draws$omega_gen)
    prior_proba_omega_above_1 <- sum(prior_draws$omega_gen > 1) / length(prior_draws$omega_gen)

    # Format output to match run_bayesian_model_by_gene.R structure
    # Creates 3 rows per branch pattern: dN, dS, omega
    bind_rows(
      # dN row
      tibble(
        posterior_mean = summary_posterior_dN["posterior_mean"],
        posterior_mode = summary_posterior_dN["posterior_mode"],
        posterior_median = summary_posterior_dN["posterior_median"],
        posterior_2.5 = summary_posterior_dN["posterior_2.5"],
        posterior_97.5 = summary_posterior_dN["posterior_97.5"],
        param = "dN",
        i_row = i_row,
        posterior_proba_omega_greater_1 = NA_real_,
        prior_proba_omega_greater_1 = NA_real_,
        marginal_LH_N = NA_real_,
        marginal_LH_S = NA_real_,
        synonymous_muts = synonymous_muts,
        non_synonymous_muts = non_synonymous_muts,
        syn_sites = syn_sites,
        nonsyn_sites = nonsyn_sites,
        n_branches = branch_row[["n_branches"]],
        gene = "Genome",
        CV_prior = CV_prior
      ),
      # dS row
      tibble(
        posterior_mean = summary_posterior_dS["posterior_mean"],
        posterior_mode = summary_posterior_dS["posterior_mode"],
        posterior_median = summary_posterior_dS["posterior_median"],
        posterior_2.5 = summary_posterior_dS["posterior_2.5"],
        posterior_97.5 = summary_posterior_dS["posterior_97.5"],
        param = "dS",
        i_row = i_row,
        posterior_proba_omega_greater_1 = NA_real_,
        prior_proba_omega_greater_1 = NA_real_,
        marginal_LH_N = NA_real_,
        marginal_LH_S = NA_real_,
        synonymous_muts = synonymous_muts,
        non_synonymous_muts = non_synonymous_muts,
        syn_sites = syn_sites,
        nonsyn_sites = nonsyn_sites,
        n_branches = branch_row[["n_branches"]],
        gene = "Genome",
        CV_prior = CV_prior
      ),
      # omega row
      tibble(
        posterior_mean = summary_posterior_omega["posterior_mean"],
        posterior_mode = summary_posterior_omega["posterior_mode"],
        posterior_median = summary_posterior_omega["posterior_median"],
        posterior_2.5 = summary_posterior_omega["posterior_2.5"],
        posterior_97.5 = summary_posterior_omega["posterior_97.5"],
        param = "omega",
        i_row = i_row,
        posterior_proba_omega_greater_1 = posterior_proba_omega_above_1,
        prior_proba_omega_greater_1 = prior_proba_omega_above_1,
        marginal_LH_N = NA_real_,
        marginal_LH_S = NA_real_,
        synonymous_muts = synonymous_muts,
        non_synonymous_muts = non_synonymous_muts,
        syn_sites = syn_sites,
        nonsyn_sites = nonsyn_sites,
        n_branches = branch_row[["n_branches"]],
        gene = "Genome",
        CV_prior = CV_prior
      )
    )
  }))

  return(df_param_estim)
}

# Load prior parametrization
df_prior <- read_tsv(args$df_prior_param)

# Load genome-level grouped data
genome_data_grouped <- read_tsv(args$data_viridian_grouped)

# Run genome-wide model
df_param_estim <- run_genome_model(df_prior, genome_data_grouped, args$CV_prior)

# Save results
write_tsv(df_param_estim, args$output_param_estim)
