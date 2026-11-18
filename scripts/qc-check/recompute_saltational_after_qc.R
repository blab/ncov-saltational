library(extraDistr)
library(tidyverse)
library(argparse)
library(zipfR)
library(parallel)
library(data.table)

source('scripts/bayesian-model/utils_model.R')

# Saltational status of each branch once the mutations flagged by quality control are
# removed. Only the affected branches are recomputed: every other branch, and every
# internal branch, keeps the status published by the pipeline.
#
# The gene level Bayes Factor is closed form. The genome wide one is approximated by
# Monte Carlo exactly as run_bayesian_model_genome.R does, with the same per row seeding,
# so that recomputing a branch whose counts are unchanged reproduces the published value.

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--df_counts', type = 'character', help = 'TSV file with the mutation counts per branch and gene')
  parser$add_argument('--df_prior_param', type = 'character', help = 'TSV file with prior parametrization')
  parser$add_argument('--name_saltational_branches', type = 'character', help = 'Text file with the name of the saltational branches')
  parser$add_argument('--CV_prior', type = 'double', help = 'CV for prior parametrization')
  parser$add_argument('--BF_threshold', type = 'double', help = 'Bayes Factor threshold used to define adaptive branches')
  parser$add_argument('--min_nonsyn', type = 'integer', default = 5, help = 'Minimum number of non-synonymous mutations for a saltational branch')
  parser$add_argument('--n_draws_genome', type = 'integer', default = 10000, help = 'Number of Monte Carlo draws for the genome wide Bayes Factor')
  parser$add_argument('--seed', type = 'integer', default = 1, help = 'Random seed; branch i is drawn with seed + i, as in the genome wide model')
  parser$add_argument('--n_cores', type = 'integer', default = 1, help = 'Number of cores for parallel processing')
  parser$add_argument('--output_name_saltational_branches', type = 'character', help = 'Text file with the saltational branches after filtering')
  parser$add_argument('--output_summary', type = 'character', help = 'TSV file comparing the branches before and after filtering')
  return(parser$parse_args())
}

args <- collect_args()

vec_genes <- c('ORF1a', 'ORF1b', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF9b')

# Load input data
df_counts <- fread(args$df_counts)
df_prior <- read_tsv(args$df_prior_param, show_col_types = FALSE)
vec_name_saltational_branches <- read.table(args$name_saltational_branches) %>% unlist() %>% as.character()

# Prior parameters for each gene
prior_params_list <- lapply(vec_genes, FUN = function(curr_gene){
  param_N <- get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dN[df_prior$gene == curr_gene]), cv = args$CV_prior)
  param_S <- get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dS[df_prior$gene == curr_gene]), cv = args$CV_prior)
  list('alpha_N' = unname(param_N['alpha']), 'beta_N' = unname(param_N['beta']),
       'alpha_S' = unname(param_S['alpha']), 'beta_S' = unname(param_S['beta']))
})
names(prior_params_list) <- vec_genes

# Gene level Bayes Factor for omega > 1, in closed form
get_BF_gene <- function(curr_gene, N, S, LN, LS){
  pr <- prior_params_list[[curr_gene]]
  prior_proba <- 1. - pbetapr(pr$beta_N / pr$beta_S, pr$alpha_N, pr$alpha_S)
  posterior_proba <- 1. - pbetapr((pr$beta_N + LN) / (pr$beta_S + LS), pr$alpha_N + N, pr$alpha_S + S)
  (posterior_proba / (1. - posterior_proba)) / (prior_proba / (1. - prior_proba))
}

# Only the branches that lost at least one mutation need recomputing
vec_branches_affected <- df_counts[n_nonsyn_flagged > 0 | n_syn_flagged > 0, unique(strain)]
cat(sprintf('Branches affected by the filtering: %d of %d saltational branches\n',
            length(vec_branches_affected), length(vec_name_saltational_branches)))

df_affected <- df_counts[strain %in% vec_branches_affected]

# Gene level: is the branch still adaptive in at least one gene?
df_affected[, BF_gene := get_BF_gene(gene, non_synonymous_muts_qc_filtered, synonymous_muts_qc_filtered,
                                     nonsyn_sites, syn_sites), by = gene]
df_adaptive_gene <- df_affected[, .(is_adaptive_any_gene = any(BF_gene > args$BF_threshold),
                                    genes_adaptive = paste(sort(gene[BF_gene > args$BF_threshold]), collapse = '|'),
                                    n_nonsyn_total = sum(non_synonymous_muts_qc_filtered)),
                                by = strain]

# Genome wide: Monte Carlo over the gene specific posteriors, as in the genome wide model.
# The branch is laid out in the wide format compute_genome_draws expects.
get_BF_genome <- function(curr_strain, i_branch){
  d <- df_affected[strain == curr_strain]
  branch_row <- as.list(setNames(
    c(d$non_synonymous_muts_qc_filtered, d$synonymous_muts_qc_filtered, d$nonsyn_sites, d$syn_sites),
    c(paste0('non_synonymous_muts_', d$gene), paste0('synonymous_muts_', d$gene),
      paste0('nonsyn_sites_', d$gene), paste0('syn_sites_', d$gene))))

  # Same per branch seeding as run_bayesian_model_genome.R, so that an unchanged branch
  # reproduces the published Bayes Factor
  set.seed(args$seed + i_branch)
  posterior_draws <- compute_genome_draws(branch_row, vec_genes, prior_params_list,
                                          n_rep = args$n_draws_genome, draw_type = 'posterior')
  prior_draws <- compute_genome_draws(branch_row, vec_genes, prior_params_list,
                                      n_rep = args$n_draws_genome, draw_type = 'prior')
  posterior_proba <- mean(posterior_draws$omega_gen > 1)
  prior_proba <- mean(prior_draws$omega_gen > 1)
  (posterior_proba / (1. - posterior_proba)) / (prior_proba / (1. - prior_proba))
}

cat(sprintf('Computing the genome wide Bayes Factor for %d branches with %d draws on %d cores\n',
            length(vec_branches_affected), args$n_draws_genome, args$n_cores))
vec_BF_genome <- unlist(mclapply(seq_along(vec_branches_affected), FUN = function(i_branch){
  get_BF_genome(vec_branches_affected[i_branch], i_branch)
}, mc.cores = args$n_cores))

df_recomputed <- df_adaptive_gene %>%
  left_join(tibble(strain = vec_branches_affected, BF_genome = vec_BF_genome), by = 'strain') %>%
  mutate(is_adaptive_genome = BF_genome > args$BF_threshold,
         # The saltational criteria, applied to the filtered counts
         is_saltational_qc_filtered = (is_adaptive_any_gene | is_adaptive_genome) &
           n_nonsyn_total >= args$min_nonsyn)

# Branches that were not affected keep their status
vec_saltational_qc_filtered <- c(
  setdiff(vec_name_saltational_branches, vec_branches_affected),
  df_recomputed$strain[df_recomputed$is_saltational_qc_filtered])

df_summary <- tibble(
  metric = c('saltational_branches_published',
             'saltational_branches_qc_filtered',
             'branches_lost',
             'branches_affected_by_filtering',
             'saltational_terminal_published',
             'saltational_terminal_qc_filtered',
             'saltational_internal_published',
             'saltational_internal_qc_filtered'),
  value = c(length(vec_name_saltational_branches),
            length(vec_saltational_qc_filtered),
            length(vec_name_saltational_branches) - length(vec_saltational_qc_filtered),
            length(vec_branches_affected),
            sum(! grepl('^node_', vec_name_saltational_branches)),
            sum(! grepl('^node_', vec_saltational_qc_filtered)),
            sum(grepl('^node_', vec_name_saltational_branches)),
            sum(grepl('^node_', vec_saltational_qc_filtered)))) %>%
  mutate(CV_prior = args$CV_prior, BF_threshold = args$BF_threshold)

cat('\n')
print(as.data.frame(df_summary %>% select(metric, value)), row.names = FALSE)

# Save results
write.table(sort(vec_saltational_qc_filtered), args$output_name_saltational_branches,
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write_tsv(df_summary, args$output_summary)
