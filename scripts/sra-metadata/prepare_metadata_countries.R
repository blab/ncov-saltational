library(argparse)
library(tidyverse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--sra_metadata', type = 'character', help = 'TSV file containing metadata downloaded from the SRA')
  parser$add_argument('--taxonium_metadata', type = 'character', help = 'TSV file containing metadata parsed from the viridian Taxonium tree')
  parser$add_argument('--output_metadata', type = 'character', help = 'TSV file with final metadata')
  return(parser$parse_args())
}

args <- collect_args()

# Get SRA metadata
sra_metadata <- read_tsv(args$sra_metadata,
                         col_select = c('strain', 'geo_loc_name_country_calc', 'geo_loc_name_sam', 'jattr'))

# Get Taxonium metadata
taxonium_metadata <- read_tsv(args$taxonium_metadata)

# Get country from SRA if available and from Taxonium otherwise
sra_metadata_with_country <- sra_metadata %>% 
  filter((geo_loc_name_country_calc != 'uncalculated'), ! is.na(geo_loc_name_country_calc)) %>% 
  rename(country_usher = geo_loc_name_country_calc) %>% 
  select(strain, country_usher)

sra_metadata_without_country <- sra_metadata %>% 
  filter((geo_loc_name_country_calc == 'uncalculated') | is.na(geo_loc_name_country_calc)) %>% 
  left_join(taxonium_metadata, by = 'strain') %>% 
  rename(country_usher = country) %>% 
  select(strain, country_usher)

metadata_full <- bind_rows(sra_metadata_with_country, sra_metadata_without_country)

# Save results
write_tsv(metadata_full, args$output_metadata)