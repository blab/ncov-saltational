library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  
  parser$add_argument('--df_sites', type = 'character', help = 'TSV file containing site information')
  parser$add_argument('--df_aa_muts_coded', type = 'character', help = 'TSV file containing AA muts information')
  parser$add_argument('--name_recombinants_to_remove', type = 'character', help = 'Text file with recombinant branches to remove')
  parser$add_argument('--df_parent_child', type = 'character', help = 'TSV file with child parent relationships')
  parser$add_argument('--df_prior_param', type = 'character', help = 'Output file with prior parametrization')
  return(parser$parse_args())
}

args <- collect_args()

# Load name of strains to remove from analysis
name_recombinants <- read.table(args$name_recombinants_to_remove) %>% unlist() %>% as.character()

# Load number of sites
df_site_counts <- read_tsv(args$df_sites) %>% 
  filter(! node_id %in% name_recombinants)

# Load child parent relationship to get all the nodes
df_parent_child <- read_tsv(args$df_parent_child) %>% 
  filter(strain != 'node_1', ! strain %in% name_recombinants) %>% # Remove root sequence (Not a branch) and recombinants
  select(strain, branch_length)

# Get number of coded AA mutations by gene by branch
df_aa_mut_all <- read_tsv(args$df_aa_muts_coded) %>% 
  mutate(is_synonymous = (aa_initial == aa_new)) %>% 
  filter(is_coded) %>% 
  rename(strain = node_id) %>% 
  filter(! strain %in% name_recombinants) %>% 
  group_by(strain, gene) %>% 
  summarise(n_mut_synonymous = sum(is_synonymous),
            n_mut_non_synonymous = sum(! is_synonymous))

# Add 0 mutation branches to the dataframe
vec_genes <- c('ORF1a', 'ORF1b', 'S', 'E', 'M', 'N', 'ORF3a', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'ORF9b')
df_aa_mut_by_gene_by_branch <- Reduce('bind_rows', lapply(vec_genes, FUN = function(curr_gene){
  df_parent_child %>% 
    left_join(df_aa_mut_all %>% filter(gene == curr_gene)) %>% 
    mutate(n_mut_synonymous = replace_na(n_mut_synonymous, 0),
           n_mut_non_synonymous = replace_na(n_mut_non_synonymous, 0),
           gene = replace_na(gene, curr_gene))
}))

# Get mean number of mutations per branch (whole genome)
df_mean_nb_mut_per_branch_whole_genome <- df_aa_mut_by_gene_by_branch %>% 
  summarise(mean_n_mut_synonymous = mean(n_mut_synonymous),
            mean_n_mut_non_synonymous = mean(n_mut_non_synonymous)) %>% 
  mutate(gene = 'Genome')

# Get mean number of mutations per gene
df_mean_nb_mut_per_branch_per_gene <- df_aa_mut_by_gene_by_branch %>% 
  group_by(gene) %>% 
  summarise(mean_n_mut_synonymous = mean(n_mut_synonymous),
            mean_n_mut_non_synonymous = mean(n_mut_non_synonymous))


# Get mean number of sites per branch (whole genome)
df_mean_nb_sites_per_branch_genome <- df_site_counts %>% 
  group_by(node_id) %>% 
  summarise(syn_sites = sum(syn_sites), 
            nonsyn_sites = sum(nonsyn_sites)) %>% 
  ungroup() %>% 
  summarise(mean_syn_sites = mean(syn_sites),
            mean_nonsyn_sites = mean(nonsyn_sites)) %>% 
  mutate(gene = 'Genome')

# Get mean number of sites per branch per gene
df_mean_nb_sites_per_branch_per_gene <- df_site_counts %>% 
  group_by(gene) %>% 
  summarise(mean_syn_sites = mean(syn_sites),
            mean_nonsyn_sites = mean(nonsyn_sites))

# Get mean dS and dN per genome region for prior parametrization
df_prior_param <- bind_rows(df_mean_nb_mut_per_branch_per_gene, df_mean_nb_mut_per_branch_whole_genome) %>% 
  left_join(
    bind_rows(df_mean_nb_sites_per_branch_per_gene, df_mean_nb_sites_per_branch_genome),
    by = 'gene'
  ) %>% 
  mutate(log_mean_dS = log(mean_n_mut_synonymous) - log(mean_syn_sites),
         log_mean_dN = log(mean_n_mut_non_synonymous) - log(mean_nonsyn_sites))

write_tsv(df_prior_param, args$df_prior_param)