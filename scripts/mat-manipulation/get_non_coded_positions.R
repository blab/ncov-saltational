library(tidyverse)

args <- commandArgs(trailingOnly = T)
df_aa_muts_path <- as.character(args[1])
df_stop_codons_path <- as.character(args[2])
df_aa_muts_with_non_coded_path <- as.character(args[3])

# Read AA change file
df_aa_muts <- read_tsv(df_aa_muts_path) %>% 
  mutate(aa_muts_position = as.numeric(aa_muts_position))
# Read stop codon file
df_stop_codons <- read_tsv(df_stop_codons_path)

# Get earliest stop codons by gene by branch
df_min_stop_codons_position_by_gene <- df_stop_codons %>% 
  group_by(node_id, gene) %>% 
  summarise(min_stop_codon_position = min(position))

# Add coded status for each aa mutation
df_aa_muts_with_non_codes <- df_aa_muts %>% 
  left_join(df_min_stop_codons_position_by_gene, by = c('node_id', 'gene')) %>% 
  mutate(is_coded = case_when(is.na(min_stop_codon_position) ~ T,
                              aa_muts_position > min_stop_codon_position ~ F,
                              T ~ T)) %>% 
  select(- min_stop_codon_position)

# Save results
write_tsv(df_aa_muts_with_non_codes, df_aa_muts_with_non_coded_path)
