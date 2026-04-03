library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--proximity_samples', type = 'character', help = 'Text file containing list of all proximity samples')
  parser$add_argument('--chronumental', type = 'character', help = 'TSV file containing chronumental inferred dates')
  parser$add_argument('--metadata', type = 'character', help = 'TSV file containing metadata for tips/strains')
  parser$add_argument('--output_metadata', type = 'character', help = 'TSV file for output metadata')
  return(parser$parse_args())
}

args <- collect_args()

# Load proximity samples list
vec_proximity_samples <- read.table(args$proximity_samples) %>%
  unlist() %>%
  as.character()

# Load chronumental inferred dates for internal nodes
chronumental_raw <- read_tsv(args$chronumental) %>%
  mutate(predicted_date = as.character(as.Date(predicted_date))) %>%
  rename(date = predicted_date) %>%
  select(strain, date)

# For adaptive_node_* entries, we need to map to original node_* names in chronumental
# Create a lookup table for adaptive_node_XXXXX -> node_XXXXX
adaptive_nodes <- vec_proximity_samples[startsWith(vec_proximity_samples, 'adaptive_node_')]
adaptive_node_mapping <- tibble(
  adaptive_name = adaptive_nodes,
  original_name = sub('^adaptive_node_', 'node_', adaptive_nodes)
)

# Get dates for adaptive nodes by looking up their original node_* names
chronumental_adaptive <- chronumental_raw %>%
  filter(strain %in% adaptive_node_mapping$original_name) %>%
  left_join(adaptive_node_mapping, by = c('strain' = 'original_name')) %>%
  select(strain = adaptive_name, date)

# Get dates for regular internal nodes (non-adaptive)
chronumental_regular <- chronumental_raw %>%
  filter(strain %in% vec_proximity_samples) %>%
  filter(!startsWith(strain, 'adaptive_node_'))

# Load metadata for tips (strains from Viridian)
metadata_tips <- read_tsv(args$metadata, col_select = c('strain', 'date')) %>%
  filter(strain %in% vec_proximity_samples) %>%
  filter(date != '.') %>%
  select(strain, date)

# Combine: adaptive nodes from chronumental (mapped), regular internal nodes, and tips
metadata_output <- bind_rows(
  chronumental_adaptive,
  chronumental_regular,
  metadata_tips %>% filter(!startsWith(strain, 'node_'))  # Exclude node_* from metadata (use chronumental instead)
) %>%
  distinct(strain, .keep_all = TRUE) %>%
  # Replace dates after 2024-07-01 with XXXX-XX-XX
  mutate(date = ifelse(as.Date(date) > as.Date("2024-07-01"), "XXXX-XX-XX", date))

# Verify we have dates for all samples
n_samples_total <- length(vec_proximity_samples)
n_samples_with_dates <- nrow(metadata_output)

cat(sprintf("Total proximity samples: %d\n", n_samples_total))
cat(sprintf("Samples with dates: %d\n", n_samples_with_dates))
cat(sprintf("Missing dates: %d\n", n_samples_total - n_samples_with_dates))

# Save results
write_tsv(metadata_output, args$output_metadata)
