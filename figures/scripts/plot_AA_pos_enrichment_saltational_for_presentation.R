### This script reproduces figures depicting the 
### AA positions that are more / less likely to 
### have nonsynonymous mutation in saltational 
### than non saltational branches

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

# Load sites associated with antibody escape
df_antigenic_sites <- read_csv('../../data/antigenic_sites_escape_antibodies.csv')

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
                aes(xmin = start, xmax = stop, ymin = max_y_axis - height_boxes/2*as.numeric(gene == 'RBD'), ymax = max_y_axis + height_boxes - height_boxes / 2 *as.numeric(gene == 'RBD')),
                alpha = 0.2, fill = palette_spike[c('RBD', 'NTD')]) +
      geom_text(data = coord_RBD_NTD,
                aes(x = (start + stop)/2, 
                    y = max_y_axis + height_boxes / 2 - height_boxes / 2 *as.numeric(gene == 'RBD'), label = gene),
                colour =  colorspace::darken(palette_spike[c('RBD', 'NTD')], 0.1)) +
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

get_manhattan_plot_highlight_mutation_of_interest <- function(df_OR_AA_pos, threshold_inclusion = 2, 
                                                              gene_nuc_coord, gene_aa_coord,
                                                              add_gene_regions = T, add_position_outisde_spike = T,
                                                              spike_position_of_interest = c(444, 484, 501, 681)){
  # Add nucleotide coordinates to df_OR_AA_pos
  df_OR_with_coords <- df_OR_AA_pos %>% 
    filter(threshold_for_inclusion == threshold_inclusion) %>% 
    left_join(gene_aa_coord, by = 'gene')
  
  df_OR_with_coords_for_labels <- df_OR_with_coords %>% 
    filter(gene == 'S', position %in% spike_position_of_interest) %>% 
    mutate(label = paste0('S:', position))
  
  max_y_axis <- max(-log10(df_OR_with_coords$p_value))
  
  baseline_plot <- df_OR_with_coords %>% 
    ggplot() +
    geom_point(aes(x = start_CDS + (position - 1)*3, y = - log10(p_value_BH),
                   colour = gene, alpha = (p_value_BH < 0.05), shape = (OR > 1.))) +
    geom_text_repel(data = df_OR_with_coords_for_labels, min.segment.length = 0.,
                    aes(start_CDS + (position - 1)*3, y = - log10(p_value_BH), label = label)) +
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
  
  baseline_plot
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
                aes(xmin = start, xmax = stop, ymin = max_y_axis - height_boxes/2*as.numeric(gene == 'RBD'), ymax = max_y_axis + height_boxes - height_boxes / 2 *as.numeric(gene == 'RBD')),
                alpha = 0.2, fill = palette_spike[c('RBD', 'NTD')]) +
      geom_text(data = coord_RBD_NTD,
                aes(x = (start + stop)/2, 
                    y = max_y_axis + height_boxes / 2 - height_boxes / 2 *as.numeric(gene == 'RBD'), label = gene),
                colour =  colorspace::darken(palette_spike[c('RBD', 'NTD')], 0.1)) +
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

get_manhattan_plot_antigenic_site_spike <- function(df_OR_AA_pos, threshold_inclusion = 2, 
                                                    gene_nuc_coord, gene_aa_coord,
                                                    add_gene_regions = T, add_position_outisde_spike = T){
  # Add nucleotide coordinates to df_OR_AA_pos
  df_OR_with_coords <- df_OR_AA_pos %>% 
    filter(threshold_for_inclusion == threshold_inclusion) %>% 
    left_join(gene_aa_coord, by = 'gene') %>% 
    left_join(df_antigenic_sites %>% mutate(gene = 'S', is_antigenic = T), by = c('gene', 'position')) %>% 
    mutate(is_antigenic = replace_na(is_antigenic, F))
  
  max_y_axis <- max(-log10(df_OR_with_coords$p_value))
  
  baseline_plot <- df_OR_with_coords %>% 
    ggplot() +
    geom_point(aes(x = start_CDS + (position - 1)*3, y = - log10(p_value_BH),
                   colour = is_antigenic, alpha = (p_value_BH < 0.05), shape = (OR > 1.))) +
    
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
    scale_colour_manual(breaks = c(T,F), values = c('firebrick2', 'darkgrey')) +
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
                aes(x = (start + stop)/2, y = max_y_axis + height_boxes * 3/ 2, label = gene),
                colour = palette_spike[c('Other S1', 'S2')]) +
      geom_rect(data = coord_RBD_NTD,
                aes(xmin = start, xmax = stop, ymin = max_y_axis - height_boxes/2*as.numeric(gene == 'RBD'), ymax = max_y_axis + height_boxes - height_boxes / 2 *as.numeric(gene == 'RBD')),
                alpha = 0.2, fill = palette_spike[c('RBD', 'NTD')]) +
      geom_text(data = coord_RBD_NTD,
                aes(x = (start + stop)/2, 
                    y = max_y_axis + height_boxes / 2 - height_boxes / 2 *as.numeric(gene == 'RBD'), label = gene),
                colour =  colorspace::darken(palette_spike[c('RBD', 'NTD')], 0.1)) +
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
                                     add_gene_regions = T, add_position_outisde_spike = F) +
  theme(legend.position = 'none')

manhattan_plot_antigenic <- get_manhattan_plot_antigenic_site_spike(df_OR_AA_pos, threshold_inclusion = 2, 
                                                                    gene_nuc_coord, gene_aa_coord,
                                                                    add_gene_regions = T, add_position_outisde_spike = F) +
  theme(legend.position = 'none')




png('../figure_mutation_enrichment/manhattan_plot_more_frequent_AA_pos_saltational_presentation.png',
    height = 3., width = 9., res = 350, units = 'in')
plot(manhattan_plot)
dev.off()
png('../figure_mutation_enrichment/manhattan_plot_more_frequent_AA_pos_saltational_presentation_antigenic.png',
    height = 3., width = 9., res = 350, units = 'in')
plot(manhattan_plot_antigenic)
dev.off()

