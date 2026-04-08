df_nuc_mut <- read_tsv('../../results/mutations/df_nuc_mut.tsv') %>% 
  group_by(strain, nuc_mutation) %>% 
  mutate(nuc_position = as.integer(substr(nuc_mutation, start = 2, stop = nchar(nuc_mutation) - 1)))
  mutate(is_non_coding = ((nuc_position < 266) | ))
  
coding_positions <- read_csv('../param/proteins_nucleotide_coordinates.csv') %>% 
  filter(gene != 'ORF9b') %>% 
  filter(startsWith(gene, 'ORF') | gene %in% c('S', 'E', 'M', 'N'))

vec_positions <- min(df_nuc_mut$nuc_position):max(df_nuc_mut$nuc_position)
vec_is_coding <- sapply(vec_positions, FUN = function(pos){
  is_coding <- sum((coding_positions$start <= pos) & (coding_positions$stop >= pos)) > 0
})
vec_non_coding_positions <- vec_positions[! vec_is_coding]

df_mut_non_coding <- df_nuc_mut %>% 
  filter(nuc_position %in% vec_non_coding_positions) %>% 
  mutate(is_adaptive = strain %in% vec_name_node_adapt) %>% 
  group_by(is_adaptive, nuc_position) %>% summarise(n_mut = n())

df_mut_non_coding %>% filter(n_mut > 3) %>% 
  group_by(nuc_position) %>% filter(n() > 1) %>% 
  arrange(nuc_position) %>% View()


df_nuc_mut %>% 
  filter(nuc_position %in% c(28271, 28272, 28273)) %>% 
  mutate(is_adaptive = strain %in% vec_name_node_adapt) %>% 
  filter(is_adaptive) %>% View()

df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  filter(is_coded)

df_aa_muts %>% 
  filter(node_id == 'ERR4877118')
  
  filter(node_id %in% vec_name_node_adapt, gene == 'ORF7a') %>% View()
