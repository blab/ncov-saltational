library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--df_child_parent', type = 'character', help = 'TSV file with child parent relationships')
  parser$add_argument('--name_saltational_branches', type = 'character', help = 'Name of saltational branches')
  parser$add_argument('--inferred_country', type = 'character', help = 'Inferred country information')
  parser$add_argument('--inferred_subregion', type = 'character', help = 'Inferred subregion information')
  parser$add_argument('--inferred_continent', type = 'character', help = 'Inferred continent information')
  parser$add_argument('--df_aa_muts', type = 'character', help = 'TSV file with number of AA muts')
  parser$add_argument('--sra_metadata', type = 'character', help = 'Metadata downloaded from the SRA')
  parser$add_argument('--output_metadata', type = 'character', help = '')
  
  return(parser$parse_args())
}

args <- collect_args()

# Load name of saltational events
name_saltational_branches <- read.table(args$name_saltational_branches) %>% unlist() %>% as.character()

# Load child parent relationships to get the list of all tips
df_child_parent <- read_tsv(args$df_child_parent)

# Load inferred locations information
inferred_country <- read_tsv(args$inferred_country)
inferred_subregion <- read_tsv(args$inferred_subregion)
inferred_continent <- read_tsv(args$inferred_continent)

# Load AA muts information
df_aa_muts <- read_tsv(args$df_aa_muts) %>% 
  filter(is_coded) %>% rename(strain = node_id) %>% 
  group_by(strain) %>% 
  summarise(n_synonymous = sum(aa_initial == aa_new), n_non_synonymous = sum(aa_initial != aa_new))

# Load SRA metadata and tidy age, sex and subdivision info
sra_metadata <- read_tsv(args$sra_metadata)

sra_metadata_subdivisions <-
  sra_metadata %>% filter(! is.na(geo_loc_name_sam)) %>%
  select(strain, geo_loc_name_sam) %>%
  mutate(geo_clean = str_remove_all(geo_loc_name_sam, "^\\[|\\]$"),
         geo_clean = ifelse(geo_clean %in% c('unknown', 'B.7 cidade', 'missing', 'Nasopharyngeal swab and Oropharyngeal swab', 
                                             'Nasopharyngeal swab'),
                            NA, geo_clean),
         subdivision_full = str_trim(str_extract(geo_clean, "(?<=:)\\s*.*")),
         subdivision_1 = str_trim(str_extract(subdivision_full, "^[^,]+")),
         subdivision_2 = str_trim(str_extract(subdivision_full, "(?<=,)\\s*.*")),
         subdivision_1 = case_when(geo_clean == 'USA: CA, Humboldt County' ~ 'California',
                                   T ~ subdivision_1)) %>% 
  select(strain, subdivision_full, subdivision_1, subdivision_2)

sra_metadata_age_sex <- sra_metadata %>% 
  mutate(sex = case_when(host_sex == 'male' ~ 'Male',
                         host_sex == 'female' ~ 'Female',
                         T ~ NA)) %>% 
  select(strain, host_age, sex) %>%
  mutate(age = case_when(host_age %in% c('not collected', '#VALUE!', 'Not Provided', 'Missing', 'Restricted Access',
                                         'Restricted Access Restricted access', 'Unidentified', 'missing', 'U',
                                         '-4 years', '-5 years', '-6 years', '-7 years') ~ NA,
                         grepl("^\\s*\\d+\\s*y\\s*$", host_age) ~ sub("^\\s*(\\d+)\\s*y\\s*$", "\\1", host_age),
                         grepl("^\\s*\\d+\\s*Y\\s*$", host_age) ~ sub("^\\s*(\\d+)\\s*Y\\s*$", "\\1", host_age),
                         grepl("^\\s*\\d+\\s*years\\s*$", host_age) ~ sub("^\\s*(\\d+)\\s*years\\s*$", "\\1", host_age),
                         grepl("^\\s*\\d+\\s*months\\s*$", host_age) ~ as.character(floor(as.numeric(sub("^\\s*(\\d+)\\s*months\\s*$", "\\1", host_age))/12)),
                         grepl("^\\s*\\d+\\s*M\\s*$", host_age) ~ as.character(floor(as.numeric(sub("^\\s*(\\d+)\\s*M\\s*$", "\\1", host_age))/12)),
                         grepl("^\\s*\\d+\\s*month\\s*$", host_age) ~ as.character(floor(as.numeric(sub("^\\s*(\\d+)\\s*month\\s*$", "\\1", host_age))/12)),
                         grepl("^\\s*\\d+\\s*Months\\s*$", host_age) ~ as.character(floor(as.numeric(sub("^\\s*(\\d+)\\s*Months\\s*$", "\\1", host_age))/12)),
                         grepl("^\\s*\\d+\\s*mo\\s*$", host_age) ~ as.character(floor(as.numeric(sub("^\\s*(\\d+)\\s*mo\\s*$", "\\1", host_age))/12)),
                         grepl("^\\s*\\d+\\s*d\\s*$", host_age) ~ as.character(floor(as.numeric(sub("^\\s*(\\d+)\\s*d\\s*$", "\\1", host_age))/365)),
                         grepl("^\\s*\\d+\\s*days\\s*$", host_age) ~ as.character(floor(as.numeric(sub("^\\s*(\\d+)\\s*days\\s*$", "\\1", host_age))/365)),
                         host_age == '<1' ~ '0',
                         host_age == '5 Y 4 M' ~ '5',
                         host_age == '2 Y 1 M' ~ '2',
                         host_age == '15.5y' ~ '15',
                         host_age == '25 d\\u00c3\\u00adas' ~ '25',
                         host_age == '12 d\\u00c3\\u00adas' ~ '12',
                         host_age == '80 - 89' ~ '80-89',
                         host_age == '30 - 49' ~ '30-39',
                         host_age == '40 - 49' ~ '40-49',
                         host_age == '20 - 29' ~ '20-29',
                         host_age == '50 - 59' ~ '50-59',
                         host_age == '60 - 69' ~ '60-69',
                         host_age == '0 - 9' ~ '0-9',
                         host_age == '30 - 39' ~ '30-3',
                         host_age == '10 - 19' ~ '10-19',
                         host_age == '>89' ~ '90+',
                         T ~ host_age),
         age = ifelse((! is.na(as.numeric(age))) & (as.numeric(age) > 116), NA, age),
         age = ifelse((! is.na(as.numeric(age))) & (as.numeric(age) < 0), NA, age)) %>% 
  select(strain, age, sex)


# Merge metadata
metadata_for_taxonium <- df_child_parent %>% 
  select(- node_id, - parent_name) %>% 
  left_join(inferred_country %>% select(strain, inferred_location) %>% rename(country = inferred_location), by = 'strain') %>% 
  left_join(inferred_subregion %>% select(strain, inferred_location) %>% rename(subregion = inferred_location), by = 'strain') %>% 
  left_join(inferred_continent %>% select(strain, inferred_location) %>% rename(continent = inferred_location), by = 'strain') %>% 
  mutate(is_saltational = strain %in% name_saltational_branches) %>% 
  left_join(df_aa_muts, by = 'strain') %>% 
  mutate(n_synonymous = replace_na(n_synonymous, 0), n_non_synonymous = replace_na(n_non_synonymous, 0)) %>% 
  left_join(sra_metadata_age_sex, by = 'strain') %>% 
  left_join(sra_metadata_subdivisions, by = 'strain') %>% 
  filter((! startsWith(strain, 'node_')) | is_saltational) %>% # Only keep saltational internal nodes and remove the other ones
  mutate(strain = ifelse(startsWith(strain, 'node_'), paste0('adaptive_', strain), strain)) # Rename saltational branches to match their name in the tree

# Write results
write_tsv(metadata_for_taxonium, args$output_metadata)