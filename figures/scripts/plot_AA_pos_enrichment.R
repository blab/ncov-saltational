library(tidyverse)
source('config_plot.R')

# Load gene coordinates
gene_nuc_coord <- read_csv('../param/proteins_nucleotide_coordinates.csv')
gene_aa_coord <- read_csv('../param/proteins_aa_coordinates.csv')

# Load OR by AA position
df_OR_AA_pos <- read_tsv('../../results/mutation_enrichment/df_AA_pos_enrichment_CV_1.7_BF_31.6227766.tsv')

# Manhattan plot
get_manhattan_plot <- function(df_OR_AA_pos, threshold_inclusion = 1, 
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
    
    height_boxes <- 7
    
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


manhattan_plot <- get_manhattan_plot(df_OR_AA_pos, threshold_inclusion = 1, 
                                     gene_nuc_coord, gene_aa_coord,
                                     add_gene_regions = T, add_position_outisde_spike = F)

pdf('../figure_mutation_enrichment/manhattan_plot_more_frequent_AA_pos.pdf',
    height = 4., width = 9.5)
plot(manhattan_plot)
dev.off()
png('../figure_mutation_enrichment/manhattan_plot_more_frequent_AA_pos.png',
    height = 4., width = 12., res = 350, units = 'in')
plot(manhattan_plot)
dev.off()

# Barplot of nb of enriched positions by gene
get_barplot_nb_enriched_positions <- function(df_OR_AA_pos, threshold_inclusion){
  barplot_nb_enriched_positions_by_gene <- tibble(gene = names(palette_genes)[names(palette_genes) != 'Genome']) %>% 
    left_join(
      df_OR_AA_pos %>% 
        filter(threshold_for_inclusion == threshold_inclusion) %>% filter(p_value_BH < 0.05) %>% 
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
                       expand = expansion(mult = c(0., 0.1))) +
    scale_fill_manual(breaks = names(palette_genes), values = palette_genes, name = 'Gene', 
                      guide = 'none') +
    scale_colour_manual(breaks = names(palette_genes), values = palette_genes, name = 'Gene', 
                        guide = 'none') +
    theme_classic() +
    theme(axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1., vjust = 1.))
  
  return(barplot_nb_enriched_positions_by_gene)
}

barplot_nb_enriched_positions_by_gene <- get_barplot_nb_enriched_positions(df_OR_AA_pos, 1)

pdf('../figure_mutation_enrichment/barplot_nb_enriched_positions_by_gene.pdf',
    height = 2.8, width = 4.)
plot(barplot_nb_enriched_positions_by_gene)
dev.off()
png('../figure_mutation_enrichment/barplot_nb_enriched_positions_by_gene.png',
    height = 2.8, width = 4., res = 350, units = 'in')
plot(barplot_nb_enriched_positions_by_gene)
dev.off()

## Make barplot with nb of enriched positions within S by region
## Could make the same plot with prop of time it's mutated.