#!/usr/bin/env Rscript
# Aggregate genome-wide gene-level data to simple format
# Converts wide format (gene-specific columns) to simple format (totals)

library(tidyverse)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--data_viridian_wide', type = 'character', help = 'TSV file with wide-format genome data')
parser$add_argument('--output', type = 'character', help = 'Output file with simple aggregated format')
args <- parser$parse_args()

# Define gene vector
vec_genes <- c('ORF1a', 'ORF1b', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF9b')

# Load wide-format data
df_wide <- read_tsv(args$data_viridian_wide, show_col_types = FALSE)

# Aggregate across genes
df_simple <- df_wide %>%
  mutate(
    # Sum mutations across all genes
    synonymous_muts = rowSums(select(., starts_with("synonymous_muts_")), na.rm = TRUE),
    non_synonymous_muts = rowSums(select(., starts_with("non_synonymous_muts_")), na.rm = TRUE),
    # Sum sites across all genes
    syn_sites = rowSums(select(., starts_with("syn_sites_")), na.rm = TRUE),
    nonsyn_sites = rowSums(select(., starts_with("nonsyn_sites_")), na.rm = TRUE)
  ) %>%
  # Keep only strain and aggregated columns
  select(strain, syn_sites, nonsyn_sites, synonymous_muts, non_synonymous_muts)

# Write output
write_tsv(df_simple, args$output)
