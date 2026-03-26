library(data.table)
library(argparse)
library(tidyverse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--metadata', type = 'character', help = 'TSV file containing metadata for which to clean incorrect dates')
  parser$add_argument('--output_metadata', type = 'character', help = 'TSV file with final metadata')
  return(parser$parse_args())
}

args <- collect_args()

# Load metadata 
sra_metadata <- read_tsv(args$metadata, col_select = c('strain', 'date'))

# Clean metadata
output_metadata <- sra_metadata %>% 
  mutate(date = ifelse(as.Date(date) < as.Date('2019-12-26'), 'XXXX-XX-XX', date))

# Save results
write_tsv(output_metadata, args$output_metadata)