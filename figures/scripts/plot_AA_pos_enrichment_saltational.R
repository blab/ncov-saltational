### This script reproduces figures depicting the 
### AA positions that are more / less likely to 
### have nonsynonymous mutation in saltational 
### than non saltattional branches

library(tidyverse)
library(ggrepel)
source('config_plot.R')


#############################
### LOADING MODEL RESULTS ###
#############################
# Load OR of NS mutation in saltational branches by AA position
df_OR_AA_pos <- read_tsv('../../results/mutation_enrichment/df_AA_pos_enrichment_saltational_CV_1.8_BF_20.tsv')

# Load name of saltational branches
name_saltational_branches <- 
  read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
             header = F) %>% unlist() %>% as.character()

# Load gene coordinates used for visualizations
gene_nuc_coord <- read_csv('../param/proteins_nucleotide_coordinates.csv')
gene_aa_coord <- read_csv('../param/proteins_aa_coordinates.csv')

#############################
### VISUALIZE RESULTS ###
#############################

# Figure 3A
# Manhattan plot with adjusted p-value (Benjamin Hochberg) by AA position
# threshold_inclusion is the required number of times a position has to have NS mutations
# on saltational branches to be included in the analysis
get_manhattan_plot <- function(df_OR_AA_pos, threshold_inclusion = 2, 
                               gene_nuc_coord, gene_aa_coord,
                               add_gene_regions = T, add_position_outisde_spike = T){
  # Add nucleotide coordinates to df_OR_AA_pos
  df_OR_with_coords <- df_OR_AA_pos %>% 
    filter(threshold_for_inclusion == threshold_inclusion) %>% 
    left_join(gene_aa_coord, by = 'gene')
  
  max_y_axis <- max(-log10(df_OR_with_coords$p_value))
  
  baseline_plot <- df_OR_with_coords %>% 
    ggplot() +
    geom_point(aes(x = start_CDS + (position - 1)*3, y = - log10(p_value_BH),
                   colour = gene, alpha = (p_value_BH < 0.05), shape = (OR > 1.))) +
    
    scale_x_continuous(name = 'AA position along the genome', limits = c(0, NA),
                       expand = expansion(mult = c(0., 0.01)),
                       breaks = seq(0, 29000, 2000)) +
    scale_y_continuous(name = expression(-log[10](p[value])),
                       expand = expansion(mult = c(0.01, 0.01)),
                       breaks = seq(0, max_y_axis, 10)) +
    scale_alpha_manual(breaks = c(T, F), values = c(0.8, 0.1),
                       guide = 'none') +
    scale_shape_manual(breaks = c(T, F), values = c(16, 15), name = 'Non-synonymous mutation', 
                       labels = c('More frequent in adaptive branches', 'Less frequent in adaptive branches')) +
    scale_colour_manual(breaks = names(palette_genes), values = palette_genes, name = 'Gene') +
    theme_classic() +
    theme(legend.position = 'bottom',
          legend.direction = 'horizontal',
          axis.text = element_text(size = 11),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11)) +
    guides(colour = guide_legend(ncol = 4, override.aes = list(label = '')),
           shape = guide_legend(ncol = 1))
  
  
  if(add_gene_regions){
    palette_regions_binary <- rep(c('black', 'gray80'), 20)
    
    # NSP coordinates
    coord_NSP <- gene_nuc_coord %>% filter(gene %in% paste0('NSP', 1:20)) %>% 
      group_by(gene) %>% 
      group_by(label_NSP = substr(gene, start = 4, stop = nchar(gene))) %>% 
      mutate(gene = factor(gene, levels = paste0('NSP', 1:16)))
    
    # S1 / S2 coordinates
    coord_S1S2 <- gene_nuc_coord %>% filter(gene %in% c('S1', 'S2'))
    
    # RBD / NTD coordinates
    coord_RBD_NTD <- tibble(gene = c('RBD', 'NTD'),
                            start = c(21563 + (319-1) * 3, 21563 + (13-1) * 3),
                            stop = c(21563 + (542-1) * 3 - 1, 21563 + (302-1) * 3 - 1))
    
    height_boxes <- 3
    
    baseline_plot <- baseline_plot +
      geom_rect(data = coord_NSP,
                aes(xmin = start, xmax = stop, ymin = max_y_axis + height_boxes, ymax = max_y_axis + 2 * height_boxes, fill = gene),
                alpha = 0.2) +
      geom_text(data = coord_NSP,
                aes(x = (start + stop)/2, y = max_y_axis + height_boxes * 3/ 2, label = label_NSP),
                colour = colorspace::darken(palette_regions_binary[1:15], .3)) +
      geom_rect(data = coord_S1S2,
                aes(xmin = start, xmax = stop, ymin = max_y_axis + height_boxes, ymax = max_y_axis + 2 * height_boxes),
                alpha = 0.2,
                fill = palette_spike[c('Other S1', 'S2')]) +
      geom_text(data = coord_S1S2,
                aes(x = (start + stop)/2, y = max_y_axis +height_boxes * 3/ 2, label = gene),
                colour = palette_spike[c('Other S1', 'S2')]) +
      geom_rect(data = coord_RBD_NTD,
                aes(xmin = start, xmax = stop, ymin = max_y_axis, ymax = max_y_axis + height_boxes),
                alpha = 0.2, fill = palette_spike[c('NTD', 'RBD')]) +
      geom_text(data = coord_RBD_NTD,
                aes(x = (start + stop)/2, y = max_y_axis + height_boxes / 2, label = gene),
                colour =  colorspace::darken(palette_spike[c('NTD', 'RBD')], 0.1)) +
      scale_fill_manual(values = palette_regions_binary, guide = 'none') 
  }
  if(add_position_outisde_spike){
    baseline_plot <- baseline_plot +
      geom_text_repel(data = df_OR_with_coords %>% filter(gene != 'S', p_value_BH < 0.05),
                      aes(x = start_CDS + (position - 1)*3,
                          y = - log10(p_value_BH), label = position,
                          colour = gene),
                      size = 3., max.overlaps = Inf)
    
  }
  
  return(baseline_plot)
}


manhattan_plot <- get_manhattan_plot(df_OR_AA_pos, threshold_inclusion = 2, 
                                     gene_nuc_coord, gene_aa_coord,
                                     add_gene_regions = T, add_position_outisde_spike = F)



pdf('../figure_mutation_enrichment/manhattan_plot_more_frequent_AA_pos_saltational.pdf',
    height = 4., width = 9.5)
plot(manhattan_plot)
dev.off()
png('../figure_mutation_enrichment/manhattan_plot_more_frequent_AA_pos_saltational.png',
    height = 4., width = 12., res = 350, units = 'in')
plot(manhattan_plot)
dev.off()

# Figure 3B - Barplot of nb of enriched positions by gene
get_barplot_nb_enriched_positions <- function(df_OR_AA_pos, threshold_inclusion, ymax = 60){
  barplot_nb_enriched_positions_by_gene <- tibble(gene = names(palette_genes)[names(palette_genes) != 'Genome']) %>% 
    left_join(
      df_OR_AA_pos %>% 
        filter(threshold_for_inclusion == threshold_inclusion) %>% filter(p_value_BH < 0.05, OR > 1.) %>% 
        group_by(gene) %>% summarise(n_pos = n()),
      by = 'gene'
    ) %>% 
    mutate(gene = factor(gene, levels = names(palette_genes)[names(palette_genes) != 'Genome']),
           n_pos = replace_na(n_pos, 0)) %>% 
    ggplot(aes(x = gene, y = n_pos, fill = gene)) +
    geom_bar(stat = 'identity') +
    geom_text(aes(label = n_pos, colour = gene), vjust = -0.5) +
    scale_x_discrete(name = 'Gene') +
    scale_y_continuous(name = 'Number of enriched positions', 
                       expand = expansion(mult = c(0., 0.1)),
                       limits = c(0., ymax),
                       breaks = seq(0., ymax, 20)) +
    scale_fill_manual(breaks = names(palette_genes), values = palette_genes, name = 'Gene', 
                      guide = 'none') +
    scale_colour_manual(breaks = names(palette_genes), values = palette_genes, name = 'Gene', 
                        guide = 'none') +
    theme_classic() +
    theme(axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1., vjust = 1.))
  
  return(barplot_nb_enriched_positions_by_gene)
}

barplot_nb_enriched_positions_by_gene <- get_barplot_nb_enriched_positions(df_OR_AA_pos, 2)

pdf('../figure_mutation_enrichment/barplot_nb_enriched_positions_by_gene_saltational.pdf',
    height = 2.8, width = 4.)
plot(barplot_nb_enriched_positions_by_gene)
dev.off()
png('../figure_mutation_enrichment/barplot_nb_enriched_positions_by_gene_saltational.png',
    height = 2.8, width = 4., res = 350, units = 'in')
plot(barplot_nb_enriched_positions_by_gene)
dev.off()

# Figure 3C - Barplot of nb of enriched positions within Spike
get_barplot_nb_enriched_positions_within_S <- function(df_OR_AA_pos, threshold_inclusion, ymax = 60){
  
  df_nb_positions_per_S_region <- df_OR_AA_pos %>% 
    filter(gene == 'S', threshold_for_inclusion == threshold_inclusion) %>% filter(p_value_BH < 0.05, OR > 1.) %>% 
    mutate(region = case_when(position >= 13 & position <= 302 ~ 'NTD',
                              position >= 319 & position <= 541 ~ 'RBD',
                              position >= 13 & position <= 685 ~ 'Other S1',
                              position >= 686 & position <= 1273 ~ 'S2',
                              T ~ 'Other S')) %>% 
    group_by(region) %>% summarise(n_pos = n()) %>% 
    mutate(region = factor(region, levels = c( 'Other S1', 'RBD', 'NTD', 'S2', 'Other S')),
           region_plot = case_when(region %in% c('NTD', 'RBD', 'Other S1') ~ 'S1',
                                   region == 'S2' ~ 'S2',
                                   T ~ 'Other'),
           region_plot = factor(region_plot, levels = c('S1', 'S2', 'Other')))
  
  barplot_nb_enriched_positions_by_gene <- df_nb_positions_per_S_region %>% 
    ggplot(aes(x = region_plot, y = n_pos, fill = region, group = region)) +
    geom_bar(stat = 'identity') +
    scale_x_discrete(name = 'Subunit') +
    scale_y_continuous(name = 'Number of enriched positions', 
                       expand = expansion(mult = c(0., 0.1)),
                       limits = c(0., ymax),
                       breaks = seq(0., ymax, 20)) +
    scale_fill_manual(breaks = c('NTD', 'RBD'), values = c(palette_spike['NTD'], palette_spike['RBD']), name = 'Domain', 
                      na.value = palette_genes['S']) +
    theme_classic() +
    theme(axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1., vjust = 1.))
  
  return(barplot_nb_enriched_positions_by_gene)
}

barplot_nb_enriched_positions_S_by_subunit <- get_barplot_nb_enriched_positions_within_S(df_OR_AA_pos, 2)

pdf('../figure_mutation_enrichment/barplot_nb_enriched_positions_S_by_subunit_saltational.pdf',
    height = 2.8, width = 2.8)
plot(barplot_nb_enriched_positions_S_by_subunit)
dev.off()
png('../figure_mutation_enrichment/barplot_nb_enriched_positions_S_by_subunit_saltational.png',
    height = 2.8, width = 2.8, res = 350, units = 'in')
plot(barplot_nb_enriched_positions_S_by_subunit)
dev.off()

# Figure 3D - Barplot of nb of enriched positions within ORF1ab
get_barplot_nb_enriched_positions_within_ORF1ab <- function(df_OR_AA_pos, threshold_inclusion, ymax = 60){
  
  df_nb_positions_per_NSPs <- df_OR_AA_pos %>% 
    filter(gene %in% c('ORF1a', 'ORF1b'), threshold_for_inclusion == threshold_inclusion) %>% 
    filter(p_value_BH < 0.05, OR > 1.) %>% 
    mutate(position = as.integer(position)) %>% 
    mutate(NSP = case_when((gene == 'ORF1a') & (position < 181) ~ 'NSP1',
                           (gene == 'ORF1a') & (position < 819) ~ 'NSP2',
                           (gene == 'ORF1a') & (position < 2764) ~ 'NSP3',
                           (gene == 'ORF1a') & (position < 3264) ~ 'NSP4',
                           (gene == 'ORF1a') & (position < 3570) ~ 'NSP5',
                           (gene == 'ORF1a') & (position < 3860) ~ 'NSP6',
                           (gene == 'ORF1a') & (position < 3943) ~ 'NSP7',
                           (gene == 'ORF1a') & (position < 4141) ~ 'NSP8',
                           (gene == 'ORF1a') & (position < 4254) ~ 'NSP9',
                           (gene == 'ORF1a') & (position < 4393) ~ 'NSP10',
                           (gene == 'ORF1a') ~ 'NSP12',
                           (gene == 'ORF1b') & (position < 924) ~ 'NSP12',
                           (gene == 'ORF1b') & (position < 1525) ~ 'NSP13',
                           (gene == 'ORF1b') & (position < 2052) ~ 'NSP14',
                           (gene == 'ORF1b') & (position < 2398) ~ 'NSP15',
                           (gene == 'ORF1b') & (position <= 2696) ~ 'NSP16',
                           T ~ 'Other'
    )) %>% 
    group_by(NSP, gene) %>% summarise(n_pos = n()) %>% 
    mutate(NSP = factor(NSP, levels = c(paste0('NSP', 16:1), 'Other')))
  
  barplot_nb_enriched_positions_by_gene <- df_nb_positions_per_NSPs %>% 
    ggplot(aes(x = gene, y = n_pos, fill = NSP, group = NSP)) +
    geom_bar(stat = 'identity') +
    scale_x_discrete(name = 'Gene') +
    scale_y_continuous(name = 'Number of enriched positions', 
                       expand = expansion(mult = c(0., 0.1)),
                       breaks = seq(0., ymax, 20),
                       limits = c(0., ymax)) +
    scale_fill_manual(values = palette_NSPs, limits = rev) +
    theme_classic() +
    theme(axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1., vjust = 1.))
  
  return(barplot_nb_enriched_positions_by_gene)
}

barplot_nb_enriched_positions_ORF1ab <- get_barplot_nb_enriched_positions_within_ORF1ab(df_OR_AA_pos, 2)

pdf('../figure_mutation_enrichment/barplot_nb_enriched_positions_ORF1ab_by_NSP_saltational.pdf',
    height = 2.8, width = 3.4)
plot(barplot_nb_enriched_positions_ORF1ab + guides(fill = guide_legend(ncol = 2)))
dev.off()
png('../figure_mutation_enrichment/barplot_nb_enriched_positions_ORF1ab_by_NSP_saltational.png',
    height = 2.8, width = 3.4, units = 'in', res = 350)
plot(barplot_nb_enriched_positions_ORF1ab + guides(fill = guide_legend(ncol = 2)))
dev.off()

# Table S1 - Get list of enriched positions
get_summary_enriched_positions <- function(df_OR_AA_pos, threshold_inclusion, df_aa_muts){
  df_OR_AA_pos %>% 
    filter(threshold_for_inclusion == threshold_inclusion) %>% 
    left_join(df_aa_muts, by = c('position' = 'aa_muts_position', 'gene')) %>% 
    filter(p_value_BH < 0.05, OR > 1.) %>% 
    arrange(gene) %>% 
    mutate(OR_char = paste0(round(OR, 1), ' (', round(OR_lower, 1), '-', round(OR_upper, 1), ')')) %>% 
    select(gene, position, OR, OR_char, n_mut_of_interest_adapt_branches) %>% 
    arrange(- n_mut_of_interest_adapt_branches) %>% return()
}

get_summary_enriched_positions(df_OR_AA_pos, 2, df_aa_muts) %>% 
  write_csv('../figure_mutation_enrichment/summary_enriched_positions.csv')

# Figure not included in the manuscript - OR as a function of mutation frequency
## Load AA nonsynonymoys mutations information and get proportion of adaptive branches with each mutation
df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  filter(is_coded) %>% # Coded AA mut
  filter(aa_new != aa_initial) %>% # Non-synonymous
  rename(strain = node_id) %>% 
  filter(strain %in% name_saltational_branches) %>% 
  mutate(has_mut = T) %>% 
  select(strain, gene, aa_muts_position, has_mut) %>% 
  group_by(gene, aa_muts_position) %>% 
  summarise(n_adapt_branches_with_mut = n()) %>% 
  mutate(prop_adapt_branches_with_mut = n_adapt_branches_with_mut / length(name_saltational_branches))

## Display results
get_plot_OR_func_prop_adapt_branches_labels <- function(df_OR_AA_pos, threshold_inclusion, df_aa_muts, 
                                                        vec_genes, 
                                                        prop_threshold = 0., OR_threshold = 1., 
                                                        xmin = NA, xmax = NA, ymin = NA, ymax = NA){
  
  OR_value_inf <- max(df_OR_AA_pos$OR[is.finite(df_OR_AA_pos$OR)])*2.
  
  label_125 <- function(x) {
    m <- signif(x / 10^floor(log10(x)), 12)   # mantissa
    labs <- scales::label_percent(accuracy = 0.1)(x)
    ifelse(m %in% c(1, 2, 5), labs, "")
  }
  
  plot_OR_func_prop_branches_with_mut <- df_OR_AA_pos %>% 
    filter(threshold_for_inclusion == threshold_inclusion) %>% 
    filter(gene %in% vec_genes) %>% 
    left_join(df_aa_muts, by = c('position' = 'aa_muts_position', 'gene')) %>% 
    mutate(OR = ifelse(is.infinite(OR), OR_value_inf, OR)) %>% 
    ggplot(aes(x = prop_adapt_branches_with_mut, 
               y = OR, colour = gene)) +
    geom_point(aes(alpha = p_value_BH < 0.05)) +
    geom_text_repel(data = df_OR_AA_pos %>% 
                      filter(threshold_for_inclusion == threshold_inclusion) %>% 
                      left_join(df_aa_muts, by = c('position' = 'aa_muts_position', 'gene')) %>% 
                      mutate(OR = ifelse(is.infinite(OR), OR_value_inf, OR)) %>% 
                      filter(gene %in% vec_genes) %>% 
                      filter(prop_adapt_branches_with_mut > prop_threshold, 
                             OR > OR_threshold, 
                             p_value_BH < 0.05),
                    aes(label = position),
                    size = 3.5,
                    max.time = 20, min.segment.length = 0., max.overlaps = Inf) +
    scale_x_continuous(labels = label_125,
                       trans = 'log',
                       expand = expansion(mult = c(0.065, 0.06)),
                       breaks = c(seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 9e-1, 1e-1)),
                       name = 'Proportion of adaptive branches with mutation at position')  +
    scale_y_continuous(name = 'OR',
                       trans = 'log', 
                       expand = expansion(mult = c(0.02, 0.02)),
                       breaks = c(0.1, 0.2, 0.5, 1.0, 2.0, 5, 10, 20, 50, 100, 200, OR_value_inf),
                       labels = c(0.1, 0.2, 0.5, 1.0, 2.0, 5, 10, 20, 50, 100, 200, 'Inf')
    ) +
    scale_colour_manual(breaks = names(palette_genes), values = palette_genes, name = 'Gene') +
    scale_alpha_manual(breaks = c(T, F), values = c(0.9, .05), guide = 'none') +
    theme_classic() +
    coord_cartesian(ylim = c(ymin, ymax), xlim = c(xmin, xmax)) +
    theme(axis.title = element_text(size = 11), axis.text = element_text(size = 11),
          strip.background = element_blank(), strip.text = element_blank()) +
    guides(colour = guide_legend(override.aes = list(label = '')))
  
  return(plot_OR_func_prop_branches_with_mut)
}


scatter_plot_OR_func_get_plot_OR_func_prop_branches_with_mut <- 
  get_plot_OR_func_prop_adapt_branches_labels(df_OR_AA_pos, 2, df_aa_muts,
                                              names(palette_genes), 
                                              0.01, 5,
                                              NA, NA, 1., NA)

scatter_plot_OR_func_get_plot_OR_func_prop_branches_with_mut_no_zoom <- 
  get_plot_OR_func_prop_adapt_branches_labels(df_OR_AA_pos, 2, df_aa_muts,
                                              names(palette_genes), 
                                              1., 1000,
                                              NA, NA, NA, NA)


pdf('../figure_mutation_enrichment/scatter_plot_OR_func_mut_prevalence_pos_saltational.pdf',
    height = 4., width = 6.)
plot(scatter_plot_OR_func_get_plot_OR_func_prop_branches_with_mut)
dev.off()
png('../figure_mutation_enrichment/scatter_plot_OR_func_mut_prevalence_pos_saltational.png',
    height = 4., width = 6., res = 350, units = 'in')
plot(scatter_plot_OR_func_get_plot_OR_func_prop_branches_with_mut)
dev.off()


  
  