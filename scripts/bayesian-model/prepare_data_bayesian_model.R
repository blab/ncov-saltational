library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--df_sites', type = 'character', help = 'TSV file containing site information')
  parser$add_argument('--df_aa_muts_coded', type = 'character', help = 'TSV file containing AA muts information')
  parser$add_argument('--genome_region', type = 'character', help = 'Gene region')
  parser$add_argument('--output_data_viridian', type = 'character', help = 'Output file ')
  parser$add_argument('--output_data_viridian_grouped', type = 'character', help = 'Output file ')
  return(parser$parse_args())
}

args <- collect_args()

# Load number of sites
df_site_counts <- read_tsv(args$df_sites)

# Get number of coded AA mutations by gene by branch
df_aa_mut_all <- read_tsv(args$df_aa_muts_coded) %>% 
  mutate(is_synonymous = (aa_initial == aa_new)) %>% 
  filter(is_coded) %>% 
  rename(strain = node_id) %>% 
  group_by(strain, gene) %>% 
  summarise(synonymous_muts = sum(is_synonymous),
            non_synonymous_muts = sum(! is_synonymous))

if(args$genome_region != 'Genome'){
  # Get genome region dataframe for model run and save it
  df_data_viridian_region <- df_site_counts %>% 
    filter(gene == args$genome_region) %>% 
    rename(strain = node_id) %>% 
    left_join(df_aa_mut_all %>% filter(gene == args$genome_region), by = c('strain', 'gene')) %>% 
    mutate(synonymous_muts = replace_na(synonymous_muts, 0),
           non_synonymous_muts = replace_na(non_synonymous_muts, 0))
  
  write_tsv(df_data_viridian_region, args$output_data_viridian)
  
  # Save scenarios with only the unique combinations for tge reguih
  df_data_viridian_region_by_group <- df_data_viridian_region %>% 
    group_by(synonymous_muts, non_synonymous_muts, syn_sites, nonsyn_sites) %>% 
    summarise(n_branches = n()) %>% arrange(- n_branches) %>% 
    ungroup() %>% 
    mutate(i_row = 1:n())
  
  write_tsv(df_data_viridian_region_by_group, args$output_data_viridian_grouped)
  
} else {
  
  # Genome wide dataframe
  df_whole_genome <- df_site_counts %>% 
    rename(strain = node_id) %>% 
    group_by(strain) %>% 
    summarise(syn_sites = sum(syn_sites), nonsyn_sites = sum(nonsyn_sites)) %>% 
    left_join(df_aa_mut_all %>% 
                group_by(strain) %>% 
                summarise(synonymous_muts = sum(synonymous_muts),
                          non_synonymous_muts = sum(non_synonymous_muts)), 
              by = c('strain')) %>% 
    mutate(synonymous_muts = replace_na(synonymous_muts, 0),
           non_synonymous_muts = replace_na(non_synonymous_muts, 0))
  
  write_tsv(df_whole_genome, args$output_data_viridian)
  
  df_whole_genome_by_group <- df_whole_genome %>% 
    group_by(synonymous_muts, non_synonymous_muts, syn_sites, nonsyn_sites) %>% 
    summarise(n_branches = n()) %>% arrange(- n_branches) %>% 
    ungroup() %>% 
    mutate(i_row = 1:n())
  
  write_tsv(df_whole_genome_by_group, args$output_data_viridian_grouped)
}