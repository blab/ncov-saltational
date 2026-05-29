library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--sra_metadata_dates', type = 'character', help = 'TSV file containing SRA metadata')
  parser$add_argument('--seed', type = 'character', help = 'Seed for downsampling')
  parser$add_argument('--name_subsamples', type = 'character', help = 'Name of tips to be included in the downsampled tree')
  return(parser$parse_args())
}

args <- collect_args()

# Load metadata with dates
sra_metadata_dates <- read_tsv(args$sra_metadata_dates) %>% 
  mutate(year = substr(date, 1, 4))

# Number of branches from 2024 (year with the least amount of tips)
n_branches_2024 <- sra_metadata_dates %>% 
  group_by(year) %>% summarise(n_branches = n()) %>% 
  filter(year == '2024') %>% 
  select(n_branches) %>% unlist() %>% as.integer()

# Downsample to same number of tips per year
vec_years <- c('2020', '2021', '2022', '2023', '2024')

set.seed(args$seed)
df_subsample_same_number_per_year <- Reduce('bind_rows', lapply(vec_years, FUN = function(curr_year){
  sra_metadata_dates %>% 
    filter(year == curr_year) %>% 
    sample_n(size = n_branches_2024, replace = F)
}))

# Save results
write.table(df_subsample_same_number_per_year$strain, file = args$name_subsamples, quote = F, row.names = F)