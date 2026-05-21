library(tidyverse)
library(ggpubr)
library(patchwork)

source('utils_plot_mutation_profile.R')
source('config_plot.R')

# Load gene coordinates
gene_nuc_coord <- read_csv('../param/proteins_nucleotide_coordinates.csv')
gene_aa_coord <- read_csv('../param/proteins_aa_coordinates.csv')

# Load OR by AA position
df_OR_AA_pos <- read_tsv('../../results/mutation_enrichment/df_AA_pos_enrichment_saltational_CV_1.8_BF_20.tsv')

#############
# Figure 3A #
# Manhattan plot of positions more/less likely to be mutated in saltational branches
#############
manhattan_plot <- get_manhattan_plot(df_OR_AA_pos, threshold_inclusion = 2, 
                                     gene_nuc_coord, gene_aa_coord,
                                     add_gene_regions = T, add_position_outisde_spike = F)


#############
# Figure 3B #
# Number of enriched positions by gene
#############
barplot_nb_enriched_positions_by_gene <- get_barplot_nb_enriched_positions(df_OR_AA_pos, 2)

#############
# Figure 3C #
# Number of enriched positions by Spike subunits
#############
barplot_nb_enriched_positions_S_by_subunit <- get_barplot_nb_enriched_positions_within_S(df_OR_AA_pos, 2)

#############
# Figure 3D #
# Number of enriched positions by NSPs
#############
barplot_nb_enriched_positions_ORF1ab <- get_barplot_nb_enriched_positions_within_ORF1ab(df_OR_AA_pos, 2)

#############
scatter_plot_OR_func_frequency <- plot_OR_spike_func_mut_frequency(df_OR_AA_pos, 2)

full_panel_mutation_enrichment <- manhattan_plot / 
  (barplot_nb_enriched_positions_by_gene + 
     barplot_nb_enriched_positions_S_by_subunit + 
     barplot_nb_enriched_positions_ORF1ab + 
     plot_layout(widths = c(1.6, 1., 0.9))) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = 'bold', size = 20))

pdf('../../figures/figure_mutation_profile.pdf', height = 8., width = 9.5)
plot(full_panel_mutation_enrichment)
dev.off()
