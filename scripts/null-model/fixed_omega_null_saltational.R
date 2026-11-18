library(extraDistr)
library(tidyverse)
library(argparse)
library(zipfR)
library(parallel)
library(data.table)

source('scripts/bayesian-model/utils_model.R')

# Expected number of adaptive and saltational branches under a fixed dN/dS null.
#
# Each branch keeps its observed number of coding mutations IN EACH GENE. Only the
# synonymous / non-synonymous split is redrawn, using a dN/dS fixed to that gene's
# tree-wide average:
#
#   n_nonsyn[branch, gene] ~ Binomial(n_mutations[branch, gene], p_gene)
#   p_gene = omega_gene * L_N[gene] / (omega_gene * L_N[gene] + L_S[gene])
#
# omega_gene is the observed ratio of tree-wide totals, read from the prior
# parametrization file; nothing is drawn from the Gamma prior. This is the Poisson
# likelihood of the model conditioned on the number of mutations on the branch.
#
# The null therefore asserts: given how many mutations a branch has in each gene, no
# branch is more adaptive than that gene's tree-wide average. The adaptive and
# saltational criteria are then applied exactly as they are to the observed counts.

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--dir_input_model', type = 'character', help = 'Directory with data_viridian_{gene}.tsv files')
  parser$add_argument('--df_prior_param', type = 'character', help = 'TSV file with prior parametrization')
  parser$add_argument('--CV_prior', type = 'double', help = 'CV for prior parametrization')
  parser$add_argument('--BF_threshold', type = 'double', help = 'Bayes Factor threshold used to define adaptive branches')
  parser$add_argument('--min_nonsyn', type = 'integer', default = 5, help = 'Minimum number of non-synonymous mutations for a saltational branch')
  parser$add_argument('--n_draws_genome', type = 'integer', default = 1000, help = 'Number of Monte Carlo draws per branch for the genome-wide Bayes Factor (same as the genome-wide model)')
  parser$add_argument('--n_replicates', type = 'integer', default = 100, help = 'Number of simulated replicates')
  parser$add_argument('--seed', type = 'integer', default = 1, help = 'Random seed; replicate r uses seed + r')
  parser$add_argument('--n_cores', type = 'integer', default = 1, help = 'Number of cores for parallel processing')
  parser$add_argument('--output_replicates', type = 'character', help = 'TSV file with one row per replicate and metric')
  parser$add_argument('--output_summary', type = 'character', help = 'TSV file comparing observed and simulated counts')
  parser$add_argument('--dir_output_branch_counts', type = 'character', help = 'Directory for the per-branch counts of each replicate (one gzipped TSV per replicate)')
  return(parser$parse_args())
}

args <- collect_args()

dir.create(args$dir_output_branch_counts, recursive = TRUE, showWarnings = FALSE)

vec_genes <- c('ORF1a', 'ORF1b', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF9b')
n_genes <- length(vec_genes)
n_draws_genome <- args$n_draws_genome
n_draws_prior_genome <- 2e6
vec_spike_thresholds <- 4:10
vec_quantiles <- c(0.025, 0.25, 0.5, 0.75, 0.975)   # Same as the boxplots of mutation counts per branch type

# Load prior parametrisation
df_prior <- read_tsv(args$df_prior_param, show_col_types = FALSE)

# Load observed counts per branch and gene
mat_N <- NULL
mat_S <- NULL
vec_LN <- c()
vec_LS <- c()
for(curr_gene in vec_genes){
  df_gene <- fread(paste0(args$dir_input_model, '/data_viridian_', curr_gene, '.tsv'),
                   select = c('strain', 'syn_sites', 'nonsyn_sites', 'synonymous_muts', 'non_synonymous_muts'))
  if(is.null(mat_N)){
    vec_strains <- df_gene$strain
    mat_N <- matrix(0L, nrow = nrow(df_gene), ncol = n_genes, dimnames = list(NULL, vec_genes))
    mat_S <- mat_N
  }
  stopifnot(identical(df_gene$strain, vec_strains))
  mat_N[, curr_gene] <- df_gene$non_synonymous_muts
  mat_S[, curr_gene] <- df_gene$synonymous_muts
  # Site counts vary very little across branches: use the tree-wide mean for each gene
  vec_LN[curr_gene] <- mean(df_gene$nonsyn_sites)
  vec_LS[curr_gene] <- mean(df_gene$syn_sites)
}
rm(df_gene)
n_branches <- nrow(mat_N)

# Fixed dN/dS per gene: the observed ratio of tree-wide totals. Converted to the
# probability that a mutation on a branch is non-synonymous, given the site counts.
vec_omega <- sapply(vec_genes, FUN = function(curr_gene){
  exp(df_prior$log_mean_dN[df_prior$gene == curr_gene]) / exp(df_prior$log_mean_dS[df_prior$gene == curr_gene])
})
vec_p_nonsyn <- vec_omega * vec_LN / (vec_omega * vec_LN + vec_LS)
cat('Fixed dN/dS and implied P(non-synonymous) per gene:\n')
print(data.frame(gene = vec_genes, omega = round(vec_omega, 4), p_nonsyn = round(vec_p_nonsyn, 4)), row.names = FALSE)

# Prior parameters for each gene
prior_params_list <- lapply(vec_genes, FUN = function(curr_gene){
  param_N <- get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dN[df_prior$gene == curr_gene]), cv = args$CV_prior)
  param_S <- get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dS[df_prior$gene == curr_gene]), cv = args$CV_prior)
  list('alpha_N' = unname(param_N['alpha']), 'beta_N' = unname(param_N['beta']),
       'alpha_S' = unname(param_S['alpha']), 'beta_S' = unname(param_S['beta']))
})
names(prior_params_list) <- vec_genes

# Gene-level rule: for each gene, is the Bayes Factor for omega > 1 above the threshold
# for a branch with N non-synonymous and S synonymous mutations? (closed form)
max_count <- min(max(rowSums(mat_N) + rowSums(mat_S)), 300)
list_is_adaptive <- lapply(vec_genes, FUN = function(curr_gene){
  pr <- prior_params_list[[curr_gene]]
  prior_proba <- 1. - pbetapr(pr$beta_N / pr$beta_S, pr$alpha_N, pr$alpha_S)
  outer(0:max_count, 0:max_count, FUN = function(N, S){
    posterior_proba <- 1. - pbetapr((pr$beta_N + vec_LN[curr_gene]) / (pr$beta_S + vec_LS[curr_gene]), pr$alpha_N + N, pr$alpha_S + S)
    BF <- (posterior_proba / (1. - posterior_proba)) / (prior_proba / (1. - prior_proba))
    BF > args$BF_threshold
  })
})
names(list_is_adaptive) <- vec_genes

# Genome-level rule: site-weighted average of gene-specific draws (Monte Carlo)
weight_N <- vec_LN / sum(vec_LN)
weight_S <- vec_LS / sum(vec_LS)

set.seed(args$seed)
prior_proba_genome <- mean(
  rowSums(sapply(vec_genes, FUN = function(g) weight_N[g] * rgamma(n_draws_prior_genome, prior_params_list[[g]]$alpha_N, prior_params_list[[g]]$beta_N))) /
    rowSums(sapply(vec_genes, FUN = function(g) weight_S[g] * rgamma(n_draws_prior_genome, prior_params_list[[g]]$alpha_S, prior_params_list[[g]]$beta_S))) > 1
)

is_adaptive_genome <- function(sub_N, sub_S){
  vec_is_adaptive <- logical(nrow(sub_N))
  if(nrow(sub_N) == 0) return(vec_is_adaptive)
  # Chunks sized to keep the matrices of draws small
  chunk_size <- max(1, floor(2e6 / n_draws_genome))
  list_chunks <- split(seq_len(nrow(sub_N)), ceiling(seq_len(nrow(sub_N)) / chunk_size))
  for(idx in list_chunks){
    n <- length(idx)
    dN_genome <- matrix(0., n, n_draws_genome)
    dS_genome <- matrix(0., n, n_draws_genome)
    for(g in vec_genes){
      pr <- prior_params_list[[g]]
      dN_genome <- dN_genome + weight_N[g] * matrix(rgamma(n * n_draws_genome, pr$alpha_N + sub_N[idx, g], pr$beta_N + vec_LN[g]), n, n_draws_genome)
      dS_genome <- dS_genome + weight_S[g] * matrix(rgamma(n * n_draws_genome, pr$alpha_S + sub_S[idx, g], pr$beta_S + vec_LS[g]), n, n_draws_genome)
    }
    posterior_proba <- rowMeans(dN_genome / dS_genome > 1)
    BF <- (posterior_proba / (1. - posterior_proba)) / (prior_proba_genome / (1. - prior_proba_genome))
    vec_is_adaptive[idx] <- BF > args$BF_threshold
  }
  return(vec_is_adaptive)
}

# Apply the adaptive and saltational criteria to a set of counts
get_metrics <- function(curr_N, curr_S){
  mat_adaptive <- sapply(vec_genes, FUN = function(g){
    list_is_adaptive[[g]][cbind(pmin(curr_N[, g], max_count) + 1L, pmin(curr_S[, g], max_count) + 1L)]
  })
  # Only branches with enough non-synonymous mutations can be saltational
  idx_candidates <- which(rowSums(curr_N) >= args$min_nonsyn)
  adaptive_gene <- mat_adaptive[idx_candidates, , drop = FALSE]
  adaptive_any_gene <- rowSums(adaptive_gene) > 0
  adaptive_genome <- is_adaptive_genome(curr_N[idx_candidates, , drop = FALSE], curr_S[idx_candidates, , drop = FALSE])

  # Mutation profile of saltational and other branches: number of mutations per branch
  # (quantiles) and number of mutations per gene
  is_saltational <- logical(nrow(curr_N))
  is_saltational[idx_candidates[adaptive_genome | adaptive_any_gene]] <- TRUE
  vec_nonsyn <- rowSums(curr_N)
  vec_syn <- rowSums(curr_S)
  get_profile <- function(idx, prefix){
    c(setNames(quantile(vec_nonsyn[idx] + vec_syn[idx], vec_quantiles, names = FALSE), paste0(prefix, '_all_muts_q', 100 * vec_quantiles)),
      setNames(quantile(vec_nonsyn[idx], vec_quantiles, names = FALSE), paste0(prefix, '_nonsyn_muts_q', 100 * vec_quantiles)),
      setNames(quantile(vec_syn[idx], vec_quantiles, names = FALSE), paste0(prefix, '_syn_muts_q', 100 * vec_quantiles)),
      setNames(colSums(curr_N[idx, , drop = FALSE]) + colSums(curr_S[idx, , drop = FALSE]), paste0(prefix, '_n_muts_', vec_genes)),
      setNames(colSums(curr_N[idx, , drop = FALSE]), paste0(prefix, '_n_nonsyn_', vec_genes)))
  }
  
  vec_metrics <- c(setNames(colSums(mat_adaptive), paste0('adaptive_', vec_genes)),
    setNames(colSums(adaptive_gene), paste0('saltational_via_', vec_genes)),
    'saltational_via_genome_only' = sum(adaptive_genome & !adaptive_any_gene),
    'saltational_total' = sum(adaptive_genome | adaptive_any_gene),
    'branches_with_min_nonsyn' = length(idx_candidates),
    setNames(sapply(vec_spike_thresholds, FUN = function(k) sum(curr_N[, 'S'] >= k)), paste0('branches_spike_nonsyn_ge_', vec_spike_thresholds)),
    get_profile(which(is_saltational), 'profile_saltational'),
    get_profile(which(!is_saltational), 'profile_other'))

  # Saltational status per branch, and the route by which it qualified. A branch can be
  # adaptive in several genes at once, so vec_saltational_via lists them separated by '|'.
  vec_saltational_via <- rep(NA_character_, nrow(curr_N))
  if(length(idx_candidates) > 0){
    vec_via_gene <- apply(adaptive_gene, 1, FUN = function(x) paste(vec_genes[x], collapse = '|'))
    vec_via_gene[adaptive_genome] <- ifelse(vec_via_gene[adaptive_genome] == '', 'Genome',
                                            paste0(vec_via_gene[adaptive_genome], '|Genome'))
    vec_saltational_via[idx_candidates] <- ifelse(vec_via_gene == '', NA_character_, vec_via_gene)
  }

  list('metrics' = vec_metrics,
       'is_saltational' = is_saltational,
       'saltational_via' = vec_saltational_via)
}

# Total number of mutations per branch and gene. Fixed by the null: only the
# synonymous / non-synonymous split is redrawn.
mat_total <- mat_N + mat_S

draw_counts_fixed_omega <- function(){
  new_N <- mat_N
  for(curr_gene in vec_genes){
    new_N[, curr_gene] <- rbinom(n_branches, mat_total[, curr_gene], vec_p_nonsyn[curr_gene])
  }
  new_S <- mat_total - new_N

  # The draw keeps the number of mutations of each branch in each gene
  stopifnot(all(new_N + new_S == mat_total), all(new_N >= 0L), all(new_S >= 0L))

  list('N' = new_N, 'S' = new_S)
}

# Counts drawn for each branch and gene, with the saltational status they imply.
# One gzipped file per replicate: 5.3M branches is too large to hold them all in memory.
save_branch_counts <- function(curr_N, curr_S, curr_metrics, file_path){
  df_branches <- as.data.table(curr_N)
  setnames(df_branches, paste0('n_nonsyn_', vec_genes))
  df_branches[, (paste0('n_syn_', vec_genes)) := as.data.table(curr_S)]
  df_branches[, `:=`(strain = vec_strains,
                     n_nonsyn_total = rowSums(curr_N),
                     n_syn_total = rowSums(curr_S),
                     is_saltational = curr_metrics$is_saltational,
                     saltational_via = curr_metrics$saltational_via)]
  setcolorder(df_branches, c('strain', 'n_nonsyn_total', 'n_syn_total', 'is_saltational', 'saltational_via'))
  fwrite(df_branches, file_path, sep = '\t', compress = 'gzip')
}

# Observed counts
cat(sprintf('Applying criteria to observed counts on %d branches\n', n_branches))
set.seed(args$seed)
metrics_observed <- get_metrics(mat_N, mat_S)
vec_observed <- metrics_observed$metrics
save_branch_counts(mat_N, mat_S, metrics_observed, file.path(args$dir_output_branch_counts, 'branch_counts_observed.tsv.gz'))

# Simulated counts
cat(sprintf('Running %d replicates with %d cores\n', args$n_replicates, args$n_cores))
list_replicates <- mclapply(seq_len(args$n_replicates), FUN = function(i_replicate){
  # Seed per replicate so that results do not depend on the number of cores
  set.seed(args$seed + i_replicate)
  simulated <- draw_counts_fixed_omega()
  curr_metrics <- get_metrics(simulated$N, simulated$S)
  save_branch_counts(simulated$N, simulated$S, curr_metrics,
                     file.path(args$dir_output_branch_counts, sprintf('branch_counts_replicate_%d.tsv.gz', i_replicate)))
  curr_metrics$metrics
}, mc.cores = args$n_cores)
stopifnot(all(sapply(list_replicates, is.numeric)))

df_replicates <- Reduce('bind_rows', lapply(seq_along(list_replicates), FUN = function(i_replicate){
  tibble(replicate = i_replicate, metric = names(list_replicates[[i_replicate]]), value = unname(list_replicates[[i_replicate]]))
})) %>%
  mutate(CV_prior = args$CV_prior, BF_threshold = args$BF_threshold)

df_summary <- df_replicates %>%
  group_by(metric) %>%
  summarise(null_mean = mean(value),
            null_2.5 = quantile(value, 0.025),
            null_97.5 = quantile(value, 0.975),
            n_replicates = n()) %>%
  left_join(tibble(metric = names(vec_observed), observed = unname(vec_observed)), by = 'metric') %>%
  mutate(ratio_observed_over_null = observed / null_mean,
         CV_prior = args$CV_prior, BF_threshold = args$BF_threshold) %>%
  arrange(match(metric, names(vec_observed))) %>%
  select(metric, observed, null_mean, null_2.5, null_97.5, ratio_observed_over_null, n_replicates, CV_prior, BF_threshold)

# Save results
write_tsv(df_replicates, args$output_replicates)
write_tsv(df_summary, args$output_summary)
