library(extraDistr)
library(tidyverse)
library(argparse)
library(zipfR)
library(parallel)

source('scripts/bayesian-model/utils_model.R')

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--df_prior_param', type = 'character', help = 'TSV file with prior parametrization')
  parser$add_argument('--genome_region', type = 'character', help = 'Genome region (should be "Genome")')
  parser$add_argument('--CV_prior', type = 'double', help = 'CV for prior parametrization')
  parser$add_argument('--data_viridian_grouped', type = 'character', help = 'TSV file with grouped genome data')
  parser$add_argument('--output_posterior_simul', type = 'character', help = 'TSV file for posterior simulations')
  parser$add_argument('--n_cores', type = 'integer', default = 1, help = 'Number of cores for parallel processing')
  return(parser$parse_args())
}

args <- collect_args()

# Function to run posterior predictive checks for genome-wide model
run_genome_posterior_predictive_checks <- function(df_prior, genome_data_grouped, CV_prior, n_cores = 1){

  # Define gene vector
  vec_genes <- c('ORF1a', 'ORF1b', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF9b')

  cat(sprintf("Running genome-wide posterior predictive checks with %d cores\n", n_cores))
  cat(sprintf("Processing %d unique branch patterns\n", nrow(genome_data_grouped)))

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

  # Number of replicate simulations per pattern
  n_replicate_simul <- 100

  # Process each unique branch pattern
  process_branch_ppc <- function(i_row) {
    branch_row <- genome_data_grouped[i_row, ]

    # For each replicate simulation
    replicate_results <- lapply(1:n_replicate_simul, function(i_simul) {
      # Initialize genome-wide totals
      N_genome <- 0
      S_genome <- 0
      total_LN <- 0
      total_LS <- 0

      # Sum across all genes
      for (g in vec_genes) {
        # Get gene-specific data
        data_g <- list(
          N  = branch_row[[paste0("non_synonymous_muts_", g)]],
          S  = branch_row[[paste0("synonymous_muts_", g)]],
          LN = branch_row[[paste0("nonsyn_sites_", g)]],
          LS = branch_row[[paste0("syn_sites_", g)]]
        )

        # Create gene-specific model
        model_g <- create_model(data_g, prior_params_list[[g]])

        # Sample ONE count from negative binomial posterior for this gene
        N_new_g <- model_g$simulate_from_posterior_N(n_rep = 1)
        S_new_g <- model_g$simulate_from_posterior_S(n_rep = 1)

        # Add to genome totals
        N_genome <- N_genome + N_new_g
        S_genome <- S_genome + S_new_g
        total_LN <- total_LN + data_g$LN
        total_LS <- total_LS + data_g$LS
      }

      # Return simulated genome-wide counts
      tibble(
        S_new = S_genome,
        N_new = N_genome,
        LN = total_LN,
        LS = total_LS,
        i_simul = i_simul
      )
    })

    # Combine replicates and group by unique (S_new, N_new) combinations
    bind_rows(replicate_results) %>%
      group_by(S_new, N_new, LN, LS) %>%
      summarise(n_branches_new = n(), .groups = 'drop') %>%
      mutate(i_row = i_row)
  }

  # Run simulations (parallel or serial)
  if (n_cores > 1) {
    cat(sprintf("Using parallel processing with %d cores\n", n_cores))
    results_list <- mclapply(1:nrow(genome_data_grouped),
                             process_branch_ppc,
                             mc.cores = n_cores)
    df_simul_from_posterior <- bind_rows(results_list)
  } else {
    cat("Using serial processing (1 core)\n")
    df_simul_from_posterior <- bind_rows(
      lapply(1:nrow(genome_data_grouped), process_branch_ppc)
    )
  }

  # Add metadata (match gene-level output format)
  df_simul_from_posterior <- df_simul_from_posterior %>%
    mutate(gene = "Genome", CV_prior = CV_prior)

  return(df_simul_from_posterior)
}

# Load prior parametrization
df_prior <- read_tsv(args$df_prior_param, show_col_types = FALSE)

# Load genome-level grouped data
genome_data_grouped <- read_tsv(args$data_viridian_grouped, show_col_types = FALSE)

# Run posterior predictive checks
df_simul_from_posterior <- run_genome_posterior_predictive_checks(
  df_prior,
  genome_data_grouped,
  args$CV_prior,
  args$n_cores
)

# Save results
write_tsv(df_simul_from_posterior, args$output_posterior_simul)

cat(sprintf("Saved %d simulated patterns to %s\n",
            nrow(df_simul_from_posterior),
            args$output_posterior_simul))
