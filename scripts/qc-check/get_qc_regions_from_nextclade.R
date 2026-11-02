library(tidyverse)
library(argparse)

# Regions of each terminal branch where an amino acid mutation call may be an artefact
# rather than a genuine mutation, from the Nextclade quality control output.
#
# Five sources are parsed into a common format (one row per branch, gene and codon
# interval): frameshifts, premature stop codons, clustered SNPs, mixed (ambiguous)
# sites and missing data. Frameshifts and stop codons are reported by Nextclade in
# codon positions within a gene; the other three are genomic nucleotide positions and
# are converted here.
#
# Every region Nextclade reports is kept, whatever quality control status it was given:
# the columns are only populated when there is something to report, so this is the
# conservative choice.

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--nextclade', type = 'character', help = 'TSV file with the Nextclade output for the adaptive branches')
  parser$add_argument('--gene_coordinates', type = 'character', help = 'CSV file with the amino acid coordinates of each gene')
  parser$add_argument('--output_qc_regions', type = 'character', help = 'TSV file with the regions flagged by quality control')
  return(parser$parse_args())
}

# Read arguments
args <- collect_args()

# Load input data
nextclade <- read_tsv(args$nextclade, show_col_types = FALSE)

# The file starts with a byte order mark, so the first column is renamed by position
df_gene_coordinates <- read_csv(args$gene_coordinates, show_col_types = FALSE) %>%
  rename(gene = 1)
gene_coords_aa <- df_gene_coordinates %>% select(gene, gene_length_aa = gene_length)
gene_coords_nuc <- df_gene_coordinates %>% select(gene, gene_start = start_CDS, gene_stop = end_CDS)

# Split a column of comma separated entries into one row per entry
separate_qc_entries <- function(df_nextclade, col_name){
  df_nextclade %>%
    select(seqName, entry = all_of(col_name)) %>%
    filter(! is.na(entry), entry != '') %>%
    separate_longer_delim(entry, delim = ',')
}

# Positions are reported either as a single value or as a range
get_start_from_range <- function(x) as.numeric(str_extract(x, '^[0-9]+'))
get_end_from_range <- function(x) as.numeric(str_extract(x, '(?<=-)[0-9]+$'))

# Nucleotide positions are converted to codon positions within each gene. A region is
# assigned to every gene it overlaps, so that regions in overlapping reading frames
# (N / ORF9b, ORF7a / ORF7b) are recorded under both genes.
convert_nucleotide_region_to_codons <- function(df_regions, gene_coords_nuc){
  df_regions %>%
    left_join(gene_coords_nuc, join_by(overlaps(nuc_start, nuc_end, gene_start, gene_stop))) %>%
    filter(! is.na(gene)) %>%
    mutate(aa_start = floor((pmax(nuc_start, gene_start) - gene_start) / 3) + 1,
           aa_end = floor((pmin(nuc_end, gene_stop) - gene_start) / 3) + 1)
}

# Frameshifts: codon positions within a gene, as GENE:start-end.
# An empty range means the gene product is lost entirely, so the whole gene is flagged.
df_frameshifts <- separate_qc_entries(nextclade, 'frameShifts') %>%
  separate_wider_delim(entry, delim = ':', names = c('gene', 'range')) %>%
  mutate(is_empty_range = range == 'empty range',
         aa_start = get_start_from_range(range),
         aa_end = get_end_from_range(range),
         # A frameshift at a single codon is reported without a range
         aa_end = if_else(is.na(aa_end) & ! is.na(aa_start), aa_start, aa_end)) %>%
  left_join(gene_coords_aa, by = 'gene') %>%
  mutate(aa_start = if_else(is_empty_range, 1, aa_start),
         aa_end = if_else(is_empty_range, gene_length_aa, aa_end)) %>%
  transmute(seqName, gene, aa_start, aa_end, qc_reason = 'frameshift')

# Premature stop codons: everything downstream of the earliest stop is not translated
df_stop_codons <- separate_qc_entries(nextclade, 'qc.stopCodons.stopCodons') %>%
  separate_wider_delim(entry, delim = ':', names = c('gene', 'position')) %>%
  mutate(position = as.numeric(position)) %>%
  slice_min(position, by = c(seqName, gene), with_ties = FALSE) %>%
  left_join(gene_coords_aa, by = 'gene') %>%
  transmute(seqName, gene, aa_start = position + 1, aa_end = gene_length_aa,
            qc_reason = 'after_stop_codon')

# Clustered SNPs: nucleotide positions, as start-end:n_snps
df_snp_clusters <- separate_qc_entries(nextclade, 'qc.snpClusters.clusteredSNPs') %>%
  separate_wider_delim(entry, delim = ':', names = c('range', 'n_snps')) %>%
  mutate(nuc_start = get_start_from_range(range),
         nuc_end = get_end_from_range(range),
         nuc_end = if_else(is.na(nuc_end), nuc_start, nuc_end)) %>%
  convert_nucleotide_region_to_codons(gene_coords_nuc) %>%
  transmute(seqName, gene, aa_start, aa_end, qc_reason = 'snp_cluster')

# Mixed sites: nucleotide positions of ambiguous bases, as CODE:position
df_mixed_sites <- separate_qc_entries(nextclade, 'nonACGTNs') %>%
  separate_wider_delim(entry, delim = ':', names = c('ambiguity_code', 'range')) %>%
  mutate(nuc_start = get_start_from_range(range),
         nuc_end = get_end_from_range(range),
         nuc_end = if_else(is.na(nuc_end), nuc_start, nuc_end)) %>%
  convert_nucleotide_region_to_codons(gene_coords_nuc) %>%
  transmute(seqName, gene, aa_start, aa_end, qc_reason = 'mixed_site')

# Missing data: nucleotide ranges, without a gene name
df_missing_data <- separate_qc_entries(nextclade, 'missing') %>%
  mutate(nuc_start = get_start_from_range(entry),
         nuc_end = get_end_from_range(entry),
         nuc_end = if_else(is.na(nuc_end), nuc_start, nuc_end)) %>%
  convert_nucleotide_region_to_codons(gene_coords_nuc) %>%
  transmute(seqName, gene, aa_start, aa_end, qc_reason = 'missing_data')

df_qc_regions <- bind_rows(df_frameshifts, df_stop_codons, df_snp_clusters,
                           df_mixed_sites, df_missing_data) %>%
  arrange(seqName, gene, aa_start)

cat('Regions flagged by quality control:\n')
print(as.data.frame(df_qc_regions %>% count(qc_reason, name = 'n_regions')), row.names = FALSE)
cat(sprintf('\n%d regions on %d branches\n', nrow(df_qc_regions), n_distinct(df_qc_regions$seqName)))

# Save results
write_tsv(df_qc_regions, args$output_qc_regions)
