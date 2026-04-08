library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--name_adaptive_branches', type = 'character', help = 'TSV file containing chronumental inferred dates')
  parser$add_argument('--df_aa_muts_coded', type = 'character', help = 'TSV file')
  parser$add_argument('--percentile_threshold_nonsyn', type = 'character', help = 'Name of nodes with adapt branches')
  parser$add_argument('--output_name_saltational_branches', type = 'character', help = 'TSV file for output metadata')
  return(parser$parse_args())
}

args <- collect_args()

# Load name of adaptive branches
vec_adaptive_branches <- read.table(args$name_adaptive_branches) %>% unlist() %>% as.character()

# Load AA muts information
df_aa_muts <- read_tsv(args$df_aa_muts_coded) %>% 
  filter(is_coded, aa_initial != aa_new) %>% 
  group_by(node_id) %>% summarise(n_non_syn = n())

# Define threshold to get saltational branches
percentile_non_syn_branches <- quantile(df_aa_muts$n_non_syn, probs = as.numeric(args$percentile_threshold_nonsyn))

# Define saltational branches
name_saltational_branches <- df_aa_muts %>%
  filter(node_id %in% vec_adaptive_branches, n_non_syn >= percentile_non_syn_branches) %>% 
  select(node_id) %>% unlist() %>% as.character()

# Save results
write.table(name_saltational_branches, 
            args$output_name_saltational_branches, 
            col.names = F, row.names = F, quote = F)
