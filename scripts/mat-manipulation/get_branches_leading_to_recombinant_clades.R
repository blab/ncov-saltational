library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--df_child_parent', type = 'character', help = 'TSV file with child parent relationships')
  parser$add_argument('--df_clades', type = 'character', help = 'TSV file with clade associated with each node')
  parser$add_argument('--name_branches_leading_to_recombinant', type = 'character', help = 'Text file for branche names')
  return(parser$parse_args())
}

args <- collect_args()

# Load clade information
df_clades <- read_tsv(args$df_clades)

# Load child parent relationships
df_child_parent <- read_tsv(args$df_child_parent)

vec_name_recombinants <- df_clades %>% 
  left_join(df_child_parent %>% select(strain, parent_name, branch_length)) %>% 
  left_join(df_clades %>% rename(clade_parent = clade), by = c('parent_name' = 'strain')) %>% 
  filter(clade != clade_parent) %>% 
  mutate(is_recombinant = grepl('^X[A-Z0-9]*$', clade)) %>% 
  filter(is_recombinant) %>% 
  select(strain) %>% unlist() %>% as.character()

# Save results
write.table(vec_name_recombinants, args$name_branches_leading_to_recombinant, 
            col.names = F, row.names = F, quote = F)