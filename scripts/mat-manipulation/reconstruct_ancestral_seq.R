library(tidyverse)
library(argparse)
library(seqinr)
library(parallel)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--df_nuc_mut', type = 'character', help = 'TSV file containing nucleotide mutation information')
  parser$add_argument('--ref_fasta', type = 'character', help = 'FASTA file with reference sequence information')
  parser$add_argument('--nodes_to_reconstruct', type = 'character', help = 'TSV file with names of ancestral nodes to be reconstructed')
  parser$add_argument('--df_child_parent', type = 'character', help = 'TSV file with child parent relationships')
  parser$add_argument('--output_ancestral_fasta', type = 'character', help = 'FASTA file for output ancestral sequence')
  parser$add_argument('--n_cores', type = 'integer', default = 1, help = 'Number of cores to use for parallelization')
  return(parser$parse_args())
}

args <- collect_args()

# Load nucleotide mutation information
df_nuc_mut <- read_tsv(args$df_nuc_mut) %>% 
  group_by(strain, nuc_mutation) %>% 
  mutate(nuc_init = substr(nuc_mutation, 1, 1),
         nuc_new = substr(nuc_mutation, nchar(nuc_mutation), nchar(nuc_mutation)),
         nuc_position = as.integer(substr(nuc_mutation, 2, nchar(nuc_mutation) - 1))) %>% 
  ungroup()

# Load FASTA file
ref_alignment <- read.fasta(args$ref_fasta)$MN908947

# Load nodes to reconstruct 
vec_nodes_to_reconstruct <- read.table(args$nodes_to_reconstruct) %>% unlist() %>% as.character()

# Load child parent relationship
df_child_parent <- read_tsv(args$df_child_parent)

get_path_to_root <- function(node_name, root_name = 'node_1', df_child_parent){
  
  # Initialization
  curr_node <- node_name
  name_parent <- df_child_parent$parent_name[df_child_parent$strain == curr_node]
  vec_path <- c(curr_node)
  
  while(name_parent != root_name){
    curr_node <- name_parent
    name_parent <- df_child_parent$parent_name[df_child_parent$strain == curr_node]
    vec_path <- c(vec_path, curr_node)
  }
  
  return(vec_path)
}

get_mutations_from_root <- function(node_name, root_name = 'node_1', df_child_parent, df_nuc_mut){
  
  # Get path information
  path_to_root <- get_path_to_root(node_name, root_name, df_child_parent)
  df_position_within_path <- tibble(position_id = 1:length(path_to_root),
                                    strain = path_to_root)
  
  df_mut_from_root <- df_nuc_mut %>%
    filter(strain %in% path_to_root) %>%
    left_join(df_position_within_path, by = "strain") %>%
    group_by(nuc_position) %>%
    filter(position_id == min(position_id)) %>%  # Only retain last mutated nuc
    ungroup()

  return(df_mut_from_root)
}

get_alignment_node <- function(node_name, root_name, df_child_parent, df_nuc_mut, ref_alignment){
  df_mut_from_root <- get_mutations_from_root(node_name, root_name, df_child_parent, df_nuc_mut)
  
  node_alignment <- as.character(ref_alignment)
  for(i_mut in 1:nrow(df_mut_from_root)){
    node_alignment[df_mut_from_root$nuc_position[i_mut]] <- tolower(df_mut_from_root$nuc_new[i_mut])
  }
  
  node_alignment <- as.SeqFastadna(node_alignment, name = node_name, Annot = node_name)
  
  return(node_alignment)
}

# Determine number of cores to use
n_cores <- min(args$n_cores, parallel::detectCores())
message(sprintf("Using %d cores for parallelization", n_cores))

# Reconstruct ancestral sequences in parallel
list_ancestral_sequences <- mclapply(vec_nodes_to_reconstruct,
                                      FUN = function(node_name){
  toupper(get_alignment_node(node_name, root_name = 'node_1', df_child_parent, df_nuc_mut, ref_alignment))
}, mc.cores = n_cores)

# Save output
write.fasta(list_ancestral_sequences, file.out = args$output_ancestral_fasta, names = vec_nodes_to_reconstruct)
