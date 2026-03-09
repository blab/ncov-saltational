library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--nodes_translations', type = 'character', help = 'TSV file containing translations info from matUtils')
  parser$add_argument('--output_nt_mutations', type = 'character', help = 'Output file for nt mutations information')
  parser$add_argument('--output_aa_mutations', type = 'character', help = 'Output file for aa mutations information')
  return(parser$parse_args())
}

args <- collect_args()

# Read node translations
nodes_translations <- read_tsv(args$nodes_translations)

# Convert to format easier to manipulate
df_nt_by_positions <- nodes_translations %>% 
  select(node_id, nt_mutations) %>% 
  group_by(node_id) %>% 
  separate_longer_delim(cols = 'nt_mutations', delim = ';') %>% 
  separate_longer_delim(cols = 'nt_mutations', delim = ',') %>% 
  group_by(node_id, nt_mutations) %>% 
  mutate(nt_mutations_position = as.integer(substr(nt_mutations, 2, nchar(nt_mutations) - 1)),
         nt_initial = substr(nt_mutations, 1, 1),
         nt_new = substr(nt_mutations, nchar(nt_mutations), nchar(nt_mutations))) %>% 
  ungroup()

df_aa_by_positions <- nodes_translations %>% 
  select(node_id, aa_mutations) %>% 
  group_by(node_id) %>% 
  separate_longer_delim(cols = 'aa_mutations', delim = ';') %>% 
  separate_wider_delim(cols = 'aa_mutations', delim = ':', names = c('gene', 'aa_muts')) %>% 
  group_by(node_id, gene, aa_muts) %>% 
  mutate(aa_muts_position = as.integer(substr(aa_muts, 2, nchar(aa_muts) - 1)),
         aa_initial = substr(aa_muts, 1, 1),
         aa_new = substr(aa_muts, nchar(aa_muts), nchar(aa_muts))) %>% 
  ungroup()

# Save results
write_tsv(df_nt_by_positions, args$output_nt_mutations)
write_tsv(df_aa_by_positions, args$output_aa_mutations)
