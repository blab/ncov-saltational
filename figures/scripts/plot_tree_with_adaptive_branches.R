library(ggtree)
library(treeio)
library(tidyverse)
source('config_plot.R')

tree <- read.tree('../../results/tree_building/proximity_tree_CV_1.7_BF_31.6227766.nwk')

# Load name of adaptive nodes
name_adaptive_nodes <- read.table('../../results/adaptive_branches/concat_name_adaptive_branches_CV_1.7_BF_31.6227766.tsv', 
                                  header = F) %>% unlist() %>% as.character()

# Load AA muts info
aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  filter(node_id %in% name_adaptive_nodes) %>% 
  rename(strain = node_id) %>% group_by(strain) %>% 
  summarise(n_muts = n(), n_non_syn = sum(aa_initial != aa_new), n_syn = sum(aa_initial == aa_new))

# Get used to plot tree
metadata <- tibble(
  label = c(tree$tip.label, tree$node.label),  # All nodes (tips + internal)
  branch_type = case_when(
    label %in% name_adaptive_nodes ~ "Adaptive",
    TRUE ~ "Background"
  ),
  is_internal = startsWith(label, "node_") | startsWith(label, "adaptive_node_")
) %>% 
  left_join(aa_muts, by = c('label' = 'strain')) %>% 
  mutate(n_muts_bins = cut(n_muts, c(0, 5, 10, 15, 20, 30, 40, 80)),
         n_non_syn_bins = cut(n_non_syn, c(0, 5, 10, 15, 20, 30, 40, 80)),
         n_non_syn_bins = cut(n_non_syn, c(seq(0, 25, 5), 70)),
         n_syn_bins = cut(n_syn, c(0, 5, 10, 15, 20, 30, 40, 80)))

# Display tree
col_background <- 'gray70'
plot_divergence_tree_colour_non_syn <- function(){
  plt_tree <- ggtree(tree, size = 0.3, key_glyph = 'point') %<+% metadata +
    aes(colour = n_non_syn_bins) +
    geom_tippoint(aes(subset = isTip & !is.na(n_non_syn_bins)), size = 0.8, alpha = 0.8) + 
    geom_nodepoint(aes(subset = !isTip & !is.na(n_non_syn_bins)), size = 0.8, alpha = 0.8) + 
    scale_color_manual(name = 'Number of\nAA changes',
                       values = colorRampPalette(colors =  c(brewer.pal(9, 'Reds')[-1], 'gray10'))(6), 
                       na.value = col_background,
                       breaks = function(x) x[!is.na(x)]) +
    theme_tree2() +
    geom_treescale(width = 10, color = 'black', offset = 50, x = 10, y = 9000) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    guides(color = guide_legend(override.aes = list(shape = 15, size = 4)))
  return(plt_tree)
}

tree_colour_non_syn <- plot_divergence_tree_colour_non_syn() +
  theme(legend.position = 'inside', legend.position.inside = c(0.8, 0.3))

pdf('../figure_adaptive_branches/divergence_tree_with_branches.pdf', height = 6, width = 3.5)
plot(tree_colour_non_syn)
dev.off()
png('../figure_adaptive_branches/divergence_tree_with_branches.png', height = 6, width = 3.5,
    res = 350, units = 'in')
plot(tree_colour_non_syn)
dev.off()
