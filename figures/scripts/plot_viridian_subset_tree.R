library(ggtree)
library(treeio)

viridian_nwk <- read.newick('../../results/viridian/small_tree_for_plot.nwk')

plt_newick <- ggtree(viridian_nwk, color = 'gray22', size = 0.2)
plot(plt_newick)

pdf('../figure_model_illustration/viridian_tree.pdf',
    height = 6, width = 4.)
plot(plt_newick)
dev.off()
