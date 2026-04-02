library(ggtree)
library(treeio)
library(tidyverse)
source('config_plot.R')

# Load newick tree
#tree <- read.tree('../../results/tree_building/proximity_time_tree_CV_1.7_BF_31.6227766.nwk') 
tree <- read.tree('../../results/tree_building/proximity_tree_CV_1.7_BF_31.6227766.nwk')

# Load name of adaptive nodes
name_adaptive_nodes <- read.table('../../results/adaptive_branches/concat_name_adaptive_branches_CV_1.7_BF_31.6227766.tsv', 
                                  header = F) %>% unlist() %>% as.character()

# Load AA muts info
aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  filter(node_id %in% name_adaptive_nodes) %>% 
  rename(strain = node_id) %>% 
  group_by(strain) %>% 
  summarise(n_muts = n(), 
            n_non_syn = sum(aa_initial != aa_new),
            n_syn = sum(aa_initial == aa_new))


# Generate metadata file
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
         n_syn_bins = cut(n_syn, c(0, 5, 10, 15, 20, 30, 40, 80)),)
  

# Plot tree
col_background <- 'darkgrey'
plot_divergence_tree_colour_non_syn <- function(){
  plt_tree <- ggtree(tree, size = 0.3, key_glyph = 'point') %<+% metadata +
    aes(colour = n_non_syn_bins) +
    geom_tippoint(aes(subset = isTip & !is.na(n_non_syn_bins)), size = 0.8, alpha = 0.8) + 
    geom_nodepoint(aes(subset = !isTip & !is.na(n_non_syn_bins)), size = 0.8, alpha = 0.8) + 
    scale_color_manual(name = 'Number of\nnon-synonymous\nmutations',
                       values = brewer.pal(8, 'Reds')[-1], 
                       na.value = col_background,
                       breaks = function(x) x[!is.na(x)]) +
    theme_tree2() +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    guides(color = guide_legend(override.aes = list(shape = 15, size = 4)))
  return(plt_tree)
}
plot_divergence_tree_colour_syn <- function(){
  plt_tree <- ggtree(tree, size = 0.3, key_glyph = 'point') %<+% metadata +
    aes(colour = n_syn_bins) +
    geom_tippoint(aes(subset = isTip & !is.na(n_syn_bins)), size = 0.8, alpha = 0.8) + 
    geom_nodepoint(aes(subset = !isTip & !is.na(n_syn_bins)), size = 0.8, alpha = 0.8) + 
    scale_color_manual(name = 'Number of\nsynonymous\nmutations',
                       values = brewer.pal(8, 'Reds')[-1], 
                       na.value = col_background,
                       breaks = function(x) x[!is.na(x)]) +
    theme_tree2() +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    guides(color = guide_legend(override.aes = list(shape = 15, size = 4)))
  return(plt_tree)
}
plot_divergence_tree_colour_all <- function(){
  plt_tree <- ggtree(tree, size = 0.3, key_glyph = 'point') %<+% metadata +
    aes(colour = n_muts_bins) +
    geom_tippoint(aes(subset = isTip & !is.na(n_muts_bins)), size = 0.8, alpha = 0.8) + 
    geom_nodepoint(aes(subset = !isTip & !is.na(n_muts_bins)), size = 0.8, alpha = 0.8) + 
    scale_color_manual(name = 'Number of\nmutations',
                       values = brewer.pal(8, 'Reds')[-1], 
                       na.value = col_background,
                       breaks = function(x) x[!is.na(x)]) +
    theme_tree2() +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    guides(color = guide_legend(override.aes = list(shape = 15, size = 4)))
  return(plt_tree)
}

tree_colour_syn <- plot_divergence_tree_colour_syn()
tree_colour_non_syn <- plot_divergence_tree_colour_non_syn()
tree_colour_all <- plot_divergence_tree_colour_all()

plot(tree_colour_all + theme(legend.position = 'inside', legend.position.inside = c(0.8, 0.2)))
plot(tree_colour_non_syn + theme(legend.position = 'inside', legend.position.inside = c(0.8, 0.2)))
plot(tree_colour_syn + theme(legend.position = 'inside', legend.position.inside = c(0.8, 0.2)))

pdf('../spatio_temporal_pattern/divergence_tree_colour_nb_non_syn_mut.pdf', height = 6, width = 3.5)
plot(tree_colour_non_syn + theme(legend.position = 'inside', legend.position.inside = c(0.8, 0.3)))
dev.off()
png('../spatio_temporal_pattern/divergence_tree_colour_nb_non_syn_mut.png', height = 6, width = 3.5,
    res = 350, units = 'in')
plot(tree_colour_non_syn + theme(legend.position = 'inside', 
                                 legend.background = element_blank(),
                                 legend.position.inside = c(0.8, 0.3)))
dev.off()

ggpubr::ggarrange(tree_colour_all + theme(legend.position = 'inside', legend.position.inside = c(0.8, 0.2)),
                  tree_colour_non_syn + theme(legend.position = 'inside', legend.position.inside = c(0.8, 0.2)),
                  tree_colour_syn + theme(legend.position = 'inside', legend.position.inside = c(0.8, 0.2)),
                  common.legend = T, legend = 'bottom',
                  ncol = 3, nrow = 1)
