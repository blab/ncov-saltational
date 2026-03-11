library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--chronumental', type = 'character', help = 'TSV file containing chronumental inferred dates')
  parser$add_argument('--metadata', type = 'character', help = 'TSV file')
  parser$add_argument('--names_adapt_branches', type = 'character', help = 'Name of nodes with adapt branches')
  parser$add_argument('--output_metadata', type = 'character', help = 'TSV file for output metadata')
  return(parser$parse_args())
}

args <- collect_args()


# Load name of adaptive branches
vec_adaptive_branches <- read.table(args$names_adapt_branches) %>% unlist() %>% as.character()

# Load Taxonium metadata
metadata_tips <- read_tsv(args$metadata) %>% 
  filter(strain %in% vec_adaptive_branches)

vec_name_tips_lacking_date <- metadata_tips$strain[metadata_tips$date == '.']

metadata_tips <- metadata_tips %>% filter(date != '.')

# Load chronumental data and only keep internal nodes or tips with lacking dates
chronumental <- read_tsv(args$chronumental) %>% 
  filter(strain %in% vec_adaptive_branches) %>% 
  filter(startsWith(strain, 'node_') | strain %in% vec_name_tips_lacking_date)

# Aggregate date information
metadata_date_output <- 
  bind_rows(
    chronumental %>%
      mutate(predicted_date = as.character(as.Date(predicted_date))) %>% 
      rename(date = predicted_date),
    metadata_tips %>% select(strain, date)
  ) %>% 
  mutate(branch_type = 'Adaptive')

# Save results
write_tsv(metadata_date_output, args$output_metadata)