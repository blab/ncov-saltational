library(data.table)
library(argparse)
library(tidyverse)
collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--metadata', type = 'character', help = 'TSV file containing metadata downloaded from the SRA')
  parser$add_argument('--missing_metadata', type = 'character', help = 'CSV file containing the missing metadata manually curated')
  parser$add_argument('--strains_viridian', type = 'character', help = 'Text file with all viridian strains')
  parser$add_argument('--output_strains_to_keep', type = 'character', help = 'Text file with name of tips to include')
  parser$add_argument('--output_subset_metadata', type = 'character', help = 'TSV file with subset metadata')
  return(parser$parse_args())
}

args <- collect_args()

# Load metadata 
sra_metadata <- fread(args$metadata)
print('Loaded metadata')

# Load missing metadata
missing_metadata <- fread(args$missing_metadata)
print('Loaded missing metadata')

# Load viridian strains
viridian_strains <- fread(args$strains_viridian, header = FALSE)[[1]]
print('Loaded viridian strains')

# Check for missing IDs - Stop if some are missing
combined_acc <- unique(c(sra_metadata$acc, missing_metadata$acc))
missing_ids <- setdiff(viridian_strains, combined_acc)
if(length(missing_ids) > 0){
  stop(paste('ERROR:', length(missing_ids), 'SRA ids are missing. First few:', 
             paste(head(missing_ids, 5), collapse=", ")))
}
print('Finished checking for missing IDs')

# Merge SRA metadata with missing metadata
sra_metadata <- sra_metadata %>% 
  bind_rows(missing_metadata %>% mutate(collection_date_sam = as.character(collection_date_sam)))
sra_metadata <- unique(sra_metadata, by = "acc")

print('Finished merging metadata')

# Only keep one runs for biosamples with multiple runs (except for ones with NAs)
setorder(sra_metadata, biosample, na.last = TRUE)
sra_metadata[, keep := is.na(biosample) | .I == .I[1], by = biosample]
sra_metadata <- sra_metadata[keep == TRUE]
sra_metadata[, keep := NULL]

print('Finished keeping on run per biosample')

# Tidy date information
sra_metadata <- sra_metadata %>% 
  mutate(date = case_when(collection_date_sam == '2020' ~ '2020-XX-XX',
                          collection_date_sam == '2021' ~ '2021-XX-XX',
                          collection_date_sam == '2022' ~ '2022-XX-XX',
                          collection_date_sam == '2023' ~ '2023-XX-XX',
                          collection_date_sam == '2024' ~ '2024-XX-XX',
                          collection_date_sam == '2021-12' ~ '2021-12-XX',
                          collection_date_sam == '2022-04' ~ '2022-04-XX',
                          collection_date_sam == 'missing' ~ 'XXXX-XX-XX',
                          collection_date_sam == 'Missing' ~ 'XXXX-XX-XX',
                          collection_date_sam == 'not applicable' ~ 'XXXX-XX-XX',
                          collection_date_sam == 'not provided' ~ 'XXXX-XX-XX',
                          is.na(collection_date_sam) ~ 'XXXX-XX-XX',
                          grepl("^\\[", collection_date_sam) ~ sub("^\\[", "", collection_date_sam),
                          T ~ collection_date_sam),
         date = ifelse(grepl("\\]$", date), sub("\\]$", "", date), date))

print('Finished tidying dates')

# Isolation sources to remove
vec_isolation_sources_to_remove <- c('Human nasal airway epithelial cell cultures', 
                                     'Human bronchial airway epithelial cell cultures',
                                     'supernatant from infected VeroE6-TMPRSS2 cells',
                                     'wastewater', 'Human precision-cut lung slice cultures',
                                     'Cell Culture Propagation', 
                                     'Environmental Surface', 
                                     'Propagated Cell Line', 
                                     'Reverse genetics derived infectious clone',
                                     'Vacuum Bag Dust, Building01', 
                                     'Vacuum Bag Dust, Building03', 
                                     'Vacuum Bag Dust, Building04', 
                                     'Vacuum Bag Dust, Building05',
                                     'Vacuum Bag Dust, Building11',
                                     'Vacuum Bag Dust, Building14',
                                     'Vacuum Bag Dust, Building16',
                                     'Vacuum Bag Dust, Building18',
                                     'Vacuum Bag Dust, Building19',
                                     'Vacuum Bag Dust, Building21', 
                                     'Vacuum Bag Dust, Building22',
                                     'Vacuum Bag Dust, Building24',
                                     'Vacuum Bag Dust, Building26',
                                     'Vacuum Bag Dust, Building27',
                                     'Vacuum Bag Dust, Building30',
                                     'Vacuum Bag Dust, Building31',
                                     'Vacuum Bag Dust, Building32',
                                     'Vacuum Bag Dust, Building33',
                                     'Vacuum Bag Dust, Building36',
                                     'Vacuum Bag Dust, Building35-11',
                                     'Cell culture', 
                                     'Cell free ', 
                                     'Cercopithecus aethiops',
                                     'Placenta tissue', 'Wastewater treatment plant',
                                     'cytobrush', 'feces', 'Faces', 'urine',
                                     'Pooled nasal swab', 'lung', 'BSL3',
                                     'bronchoalveolar lavage', 'Bronchial washings', 
                                     'Obtained from BEI Resources (NR-52281)', 
                                     'Obtained from BEI Resources (NR-55329)', 
                                     'nasal swab from Felis catus',
                                     'Cell free RNA', 'amniotic fluid')

metadata_tidy <- sra_metadata %>% filter(! isolation_source %in% vec_isolation_sources_to_remove)

print('Finished filtering isolation sources')

# Save results
# Rename 'acc' to 'strain' for compatibility with chronumental
metadata_tidy <- metadata_tidy %>% rename(strain = acc)
write_tsv(metadata_tidy, args$output_subset_metadata)
write.table(metadata_tidy$strain, args$output_strains_to_keep, col.names = F, row.names = F, quote = F)
