library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--df_param_estim', type = 'character', help = 'TSV file with model run')
  parser$add_argument('--df_BF', type = 'character', help = 'TSV file with BF values')
  return(parser$parse_args())
}

args <- collect_args()

# Load model run and compute BF
df_param_estim_with_BF_omega_greater_1 <- read_tsv(args$df_param_estim) %>% 
  filter(param == 'omega') %>% 
  mutate(log_posterior_odds_omega_greater_1 = log(posterior_proba_omega_greater_1) - log(1. - posterior_proba_omega_greater_1),
         log_prior_odds_omega_greater_1 = log(prior_proba_omega_greater_1) - log(1. - prior_proba_omega_greater_1),
         log_BF_omega_greater_1 = log_posterior_odds_omega_greater_1 - log_prior_odds_omega_greater_1,
         BF_omega_greater_1 = exp(log_BF_omega_greater_1))

write_tsv(df_param_estim_with_BF_omega_greater_1, args$df_BF)
