library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--metadata', type = 'character', help = 'TSV file')
  parser$add_argument('--names_adapt_branches', type = 'character', help = 'Name of nodes with adapt branches')
  parser$add_argument('--output_metadata', type = 'character', help = 'TSV file for output metadata')
  parser$add_argument('--vec_tips_non_adaptive_with_dates', type = 'character', help = 'Text file for strain names')
  return(parser$parse_args())
}

args <- collect_args()

# Load name of adaptive branches
vec_adaptive_branches <- read.table(args$names_adapt_branches) %>% unlist() %>% as.character()

# Load metadata
metadata_non_adaptive <- read_tsv(args$metadata) %>% 
  filter(date != '.', strain != 'MN908947', ! strain %in% vec_adaptive_branches) %>% 
  select(strain, date) %>% 
  mutate(branch_type = 'Non adaptive')

# Save results
write_tsv(metadata_non_adaptive, args$output_metadata)
write.table(metadata_non_adaptive$strain, args$vec_tips_non_adaptive_with_dates,
            row.names = FALSE, col.names = FALSE, quote = FALSE)