library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--name_adaptive_branches', type = 'character', help = 'TSV file containing chronumental inferred dates')
  parser$add_argument('--df_aa_muts_coded', type = 'character', help = 'TSV file')
  parser$add_argument('--percentile_threshold_nonsyn', type = 'character', help = 'Name of nodes with adapt branches')
  parser$add_argument('--df_aa_muts_coded_reference', type = 'character', default = NULL,
                      help = 'Optional TSV file with AA muts used to compute the non-syn percentile threshold (default: --df_aa_muts_coded)')
  parser$add_argument('--output_name_saltational_branches', type = 'character', help = 'TSV file for output metadata')
  return(parser$parse_args())
}

args <- collect_args()

# Load name of adaptive branches
if(file.info(args$name_adaptive_branches)$size > 0){
  
  vec_adaptive_branches <- read.table(args$name_adaptive_branches) %>% unlist() %>% as.character()
  
  # Load AA muts information
  df_aa_muts <- read_tsv(args$df_aa_muts_coded) %>% 
    filter(is_coded, aa_initial != aa_new) %>% 
    group_by(node_id) %>% summarise(n_non_syn = n())
  
  # Load AA muts used to define the threshold (e.g. observed tree when calling saltational branches on a null tree)
  if(is.null(args$df_aa_muts_coded_reference)){
    df_aa_muts_reference <- df_aa_muts
  } else{
    df_aa_muts_reference <- read_tsv(args$df_aa_muts_coded_reference) %>% 
      filter(is_coded, aa_initial != aa_new) %>% 
      group_by(node_id) %>% summarise(n_non_syn = n())
  }
  
  # Define threshold to get saltational branches
  percentile_non_syn_branches <- quantile(df_aa_muts_reference$n_non_syn, probs = as.numeric(args$percentile_threshold_nonsyn))
  print(paste0('Threshold on number of non-synonymous mutations: ', percentile_non_syn_branches))
  
  # Define saltational branches
  name_saltational_branches <- df_aa_muts %>%
    filter(node_id %in% vec_adaptive_branches, n_non_syn >= percentile_non_syn_branches) %>% 
    select(node_id) %>% unlist() %>% as.character()
  
  # Save results
  write.table(name_saltational_branches, 
              args$output_name_saltational_branches, 
              col.names = F, row.names = F, quote = F)
  
} else{
  write.table(character(0), args$output_name_saltational_branches, 
              col.names = F, row.names = F, quote = F)
}
