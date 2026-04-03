library(tidyverse)
source('config_plot.R')

vec_genes <- names(palette_genes)

df_adapt_branches_per_gene <- Reduce('bind_rows', lapply(vec_genes, FUN = function(curr_gene){
  read_tsv(paste0('../../results/adaptive_branches/name_adaptive_branches_', curr_gene, '_CV_1.7_BF_31.6227766.tsv'),
           col_names = F) %>% 
    summarise(n_adaptive_branches = n()) %>% 
    mutate(gene = curr_gene)
}))
  
barplot_nb_adaptive_branches <- df_adapt_branches_per_gene %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  ggplot(aes(x = gene, y = n_adaptive_branches, fill = gene)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = n_adaptive_branches, color = gene), vjust = -.5) +
  scale_x_discrete(name = '') +
  scale_y_continuous(name = 'Number of adaptive branches',
                     limits = c(0., NA), expand = expansion(mult = c(0., 0.12))) +
  scale_fill_manual(name = 'Gene', breaks = names(palette_genes),
                    values = palette_genes, guide = 'none') +
  scale_colour_manual(name = 'Gene', breaks = names(palette_genes),
                      values = palette_genes, guide = 'none') +
  theme_classic() +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1., hjust = 1.))

pdf('../figure_adaptive_branches/barplot_nb_adaptive_branches_per_gene.pdf',
    height = 2.8, width = 4.)
plot(barplot_nb_adaptive_branches)
dev.off()
png('../figure_adaptive_branches/barplot_nb_adaptive_branches_per_gene.png',
    height = 2.8, width = 4., res = 350, units = 'in')
plot(barplot_nb_adaptive_branches)
dev.off()
