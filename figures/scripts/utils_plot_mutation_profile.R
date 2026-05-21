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
          axis.text.x = element_text(angle = 45, hjust = 1., vjust = 1.)) +
    theme(legend.position = 'inside',
          legend.background = element_blank(),
          legend.position.inside = c(0.8, 0.8))
  
  return(barplot_nb_enriched_positions_by_gene)
}

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
          axis.text.x = element_text(angle = 45, hjust = 1., vjust = 1.)) +
    guides(fill = guide_legend(ncol = 2)) +
    theme(legend.position = 'inside',
          legend.background = element_blank(),
          legend.position.inside = c(0.75, 0.75))
  
  return(barplot_nb_enriched_positions_by_gene)
}

plot_OR_spike_func_mut_frequency <- function(df_OR_AA_pos, threshold_inclusion){
  df_plot_enriched_pos_Spike <- df_OR_AA_pos %>% filter(threshold_for_inclusion == threshold_inclusion, gene == 'S') %>% 
    mutate(region = case_when(position >= 13 & position <= 302 ~ 'NTD',
                              position >= 319 & position <= 541 ~ 'RBD',
                              position >= 13 & position <= 685 ~ 'Other S1',
                              position >= 686 & position <= 1273 ~ 'S2',
                              T ~ 'Other S')) %>% 
    mutate(region = factor(region, levels = c('Other S1', 'RBD', 'NTD', 'S2', 'Other S')),
           region_plot = case_when(region %in% c('NTD', 'RBD', 'Other S1') ~ 'S1',
                                   region == 'S2' ~ 'S2',
                                   T ~ 'Other'),
           region_plot = factor(region_plot, levels = c('S1', 'S2', 'Other')))
  
  scatter_plot <- df_plot_enriched_pos_Spike %>% 
    mutate(region = factor(region, levels = c('NTD', 'RBD', 'Other S1', 'S2', 'Other S'))) %>% 
    ggplot(aes(x = n_mut_of_interest_adapt_branches, y = OR, colour = region)) +
    geom_point(aes(alpha = p_value_BH < 0.05)) +
    geom_text_repel(data = df_plot_enriched_pos_Spike %>% 
                      filter(p_value_BH < 0.05, OR > 1.),
                    aes(label = position),
                    max.overlaps = Inf, min.segment.length = 0.0, max.time = 20) +
    scale_x_continuous(trans = 'log',
                       name = 'Number of saltational branches\nwith mutation at position',
                       breaks = c(seq(1, 9, 1), seq(10, 90, 10)),
                       labels = c(1, 2, '', '', 5, rep('', 4), 10, 20, '', '', 50, rep('', 4))) +
    scale_y_continuous(trans = 'log', 
                       name = expression(OR['mutation in saltational branches']),
                       breaks = c(seq(1, 9, 1), seq(10, 90, 10)),
                       labels = c(1, 2, '', '', 5, rep('', 4), 10, 20, '', '', 50, rep('', 4))) +
    scale_colour_manual(name = 'Spike region', values = palette_spike) +
    scale_alpha_manual(breaks = c(T, F), values = c(0.8, 0.1),
                       guide = 'none') +
    theme_classic() +
    coord_cartesian(ylim = c(1, NA)) +
    theme(axis.text = element_text(size = 11)) +
    guides(colour = guide_legend(override.aes = list(label = '')))
  
  return(scatter_plot)
}
