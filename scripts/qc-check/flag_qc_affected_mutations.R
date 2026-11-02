library(tidyverse)
library(argparse)
library(data.table)

# Amino acid mutations of the saltational branches that fall in a region flagged by
# quality control, and the mutation counts that remain once they are removed.

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--qc_regions', type = 'character', help = 'TSV file with the regions flagged by quality control')
  parser$add_argument('--df_aa_muts', type = 'character', help = 'TSV file with AA muts info')
  parser$add_argument('--name_saltational_branches', type = 'character', help = 'Text file with the name of the saltational branches')
  parser$add_argument('--dir_input_model', type = 'character', help = 'Directory with the data_viridian_{gene}.tsv files')
  parser$add_argument('--output_flagged_mutations', type = 'character', help = 'TSV file with the flagged mutations')
  parser$add_argument('--output_counts', type = 'character', help = 'TSV file with the mutation counts per branch and gene')
  return(parser$parse_args())
}

args <- collect_args()

vec_genes <- c('ORF1a', 'ORF1b', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF9b')

# Load input data
df_qc_regions <- fread(args$qc_regions)
vec_name_saltational_branches <- read.table(args$name_saltational_branches) %>% unlist() %>% as.character()

# Synonymous mutations are kept: a flagged region can contain either kind, and removing
# only the non-synonymous ones would raise dN/dS, biasing the filter in one direction.
df_aa_muts <- fread(args$df_aa_muts)[is_coded == TRUE &
                                       node_id %in% vec_name_saltational_branches]
df_aa_muts[, is_synonymous := aa_initial == aa_new]

# Flag the mutations falling in a flagged region, matching on branch, gene and codon
dt_regions <- df_qc_regions[seqName %in% vec_name_saltational_branches]
setkey(dt_regions, seqName, gene, aa_start, aa_end)

dt_mutations <- df_aa_muts[, .(seqName = node_id, gene, aa_muts, aa_muts_position, is_synonymous,
                               aa_start = aa_muts_position, aa_end = aa_muts_position)]

df_overlaps <- foverlaps(dt_mutations, dt_regions, type = 'within', nomatch = NULL)

# A mutation can fall in several flagged regions, so the reasons are collapsed into a
# single row per mutation. Without this the same mutation would be counted several times.
df_flagged_mutations <- df_overlaps[, .(qc_reasons = paste(sort(unique(qc_reason)), collapse = '|'),
                                        n_qc_reasons = uniqueN(qc_reason)),
                                    by = .(node_id = seqName, gene, aa_muts, aa_muts_position, is_synonymous)]

# Mutation counts per branch and gene, before and after removing the flagged mutations.
# The site counts come from the model input files and are left unchanged.
df_counts_observed <- rbindlist(lapply(vec_genes, FUN = function(curr_gene){
  fread(file.path(args$dir_input_model, paste0('data_viridian_', curr_gene, '.tsv')),
        select = c('strain', 'syn_sites', 'nonsyn_sites', 'synonymous_muts', 'non_synonymous_muts'))[
          strain %in% vec_name_saltational_branches][, gene := curr_gene]
}))

df_n_flagged <- df_flagged_mutations[, .(n_nonsyn_flagged = sum(! is_synonymous),
                                         n_syn_flagged = sum(is_synonymous)),
                                     by = .(strain = node_id, gene)]

df_counts <- merge(df_counts_observed, df_n_flagged, by = c('strain', 'gene'), all.x = TRUE)
df_counts[is.na(n_nonsyn_flagged), n_nonsyn_flagged := 0]
df_counts[is.na(n_syn_flagged), n_syn_flagged := 0]
df_counts[, non_synonymous_muts_qc_filtered := non_synonymous_muts - n_nonsyn_flagged]
df_counts[, synonymous_muts_qc_filtered := synonymous_muts - n_syn_flagged]
df_counts[, is_terminal := ! grepl('^node_', strain)]

# A branch cannot lose more mutations than it carries
stopifnot(all(df_counts$non_synonymous_muts_qc_filtered >= 0),
          all(df_counts$synonymous_muts_qc_filtered >= 0))
# Internal branches are reconstructed sequences and have no Nextclade quality control
stopifnot(all(df_counts[is_terminal == FALSE, n_nonsyn_flagged] == 0),
          all(df_counts[is_terminal == FALSE, n_syn_flagged] == 0))

cat(sprintf('Mutations on saltational branches: %d non-synonymous, %d synonymous\n',
            sum(! df_aa_muts$is_synonymous), sum(df_aa_muts$is_synonymous)))
cat(sprintf('Flagged by quality control: %d non-synonymous (%.1f%%), %d synonymous (%.1f%%)\n',
            sum(! df_flagged_mutations$is_synonymous),
            100 * sum(! df_flagged_mutations$is_synonymous) / sum(! df_aa_muts$is_synonymous),
            sum(df_flagged_mutations$is_synonymous),
            100 * sum(df_flagged_mutations$is_synonymous) / sum(df_aa_muts$is_synonymous)))
cat('\nFlagged mutations by reason:\n')
print(as.data.frame(df_overlaps[, .(n_mutations = uniqueN(paste(seqName, gene, aa_muts))), by = qc_reason][order(-n_mutations)]),
      row.names = FALSE)

# Save results
fwrite(df_flagged_mutations, args$output_flagged_mutations, sep = '\t')
fwrite(df_counts, args$output_counts, sep = '\t')
