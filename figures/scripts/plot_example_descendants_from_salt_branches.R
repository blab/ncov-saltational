library(ggtree)
library(treeio)
library(ape)
library(tidyverse)
library(ggtext)
source('config_plot.R')

# Define name of clusters to be extracted
name_nodes_to_plot <- c('node_176', 'node_2803', 'node_8768', 'node_9379', 
                        'node_38879', 'node_39826', 'node_812479', 
                        'node_810966', 'node_841663', 'node_861367',
                        'node_675670', 'node_7543', 'node_858747',
                        'node_870259', 'node_9009', 'node_369923')

# Get parent information
df_parents_of_nodes_to_plot <- read_tsv('../../results/bte/df_parent_child_relationship.tsv') %>% 
  filter(strain %in% name_nodes_to_plot)

# Load chronic mutation information
df_chronic_muts <- read_csv('../../data/frequent_mutations/frequent_mutations_chronic.csv') %>% 
  filter(source != 'Ghafari et al. Lancet ID. Occurred 4+ times') %>% 
  select(gene, aa_mut_with_origin) %>% mutate(chronic_like = T) %>% distinct()

# Load aa muts information
underline_text <- function(x) {
  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  paste0(chars, "\u0332", collapse = "")
}

df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  rename(strain = node_id) %>% 
  filter(is_coded, strain %in% name_nodes_to_plot) %>% 
  left_join(df_chronic_muts,
            by = c('gene', 'aa_muts' = 'aa_mut_with_origin')) %>% 
  mutate(gene = factor(gene, levels = rev(names(palette_genes))),
         chronic_like = replace_na(chronic_like, F),
         aa_muts_colored = if_else(
           chronic_like,
           paste0("**", aa_muts, "**"),
           aa_muts)
         ) %>% 
  arrange(gene) %>% 
  group_by(strain, gene) %>% 
  summarise(muts = paste(aa_muts, collapse = ','),
            muts_colored = paste(aa_muts_colored, collapse = ',')) %>% 
  mutate(muts_with_gene = paste0(gene, ': ', muts),
         muts_with_gene_colored = paste0(gene, ': ', muts_colored)) %>% 
  group_by(strain) %>% 
  summarise(muts_for_tree = paste(muts_with_gene, collapse = "\n"),
            muts_for_tree_colored = paste(muts_with_gene_colored, collapse = "<br>"),
            .groups = "drop")

# Extract subtrees from these parents
extract_newick_subtree <- function(name_salt_branch, df_parents_of_nodes_to_plot,
                                   matUtils_path = '/Users/ctrankie/miniforge3/bin/matUtils'){
  # Get parent name
  name_parent <- df_parents_of_nodes_to_plot$parent_name[df_parents_of_nodes_to_plot$strain == name_salt_branch]
  
  # Extract newick tree
  curr_command <- paste0(matUtils_path, 
    ' extract -i ../../results/viridian/viridian_tree_subset.pb.gz --get-internal-descendents ',
    name_parent, ' --write-tree ../cluster_spread_from_saltational_events/subtree_', name_salt_branch, '.nwk'
    )
  system(command = curr_command)
  
  # Load newick tree
  tree <- read.newick(paste0('../cluster_spread_from_saltational_events/subtree_', name_salt_branch, '.nwk'))
  phy_tree <- as.phylo(tree)
  sub_tree <- extract.clade(phy_tree, node = name_salt_branch, root.edge = T)
  return(sub_tree)
}

list_subtrees <- lapply(name_nodes_to_plot, FUN = function(name_salt_branch){
  extract_newick_subtree(name_salt_branch, df_parents_of_nodes_to_plot)
})


get_tree_metadata_from_node <- function(name_salt_branch){
  sub_tree <- extract_newick_subtree(name_salt_branch, df_parents_of_nodes_to_plot)
  
  name_tips <- sub_tree$tip.label
  dates_tips <- read_tsv('../../results/sra_metadata/sra_metadata_just_clean_dates.tsv', col_select = c('strain', 'date')) %>% 
    filter(strain %in% name_tips)
  metadata_tips <- read_tsv('../../results/saltational_branches/metadata_taxonium_CV_1.8_BF_20.tsv') %>% 
    filter(strain %in% name_tips) %>% 
    left_join(dates_tips, by = 'strain') %>% 
    mutate(label_age_sex = case_when(! is.na(age) & ! is.na(sex) ~ paste0(age, '', substr(sex, 1, 1)),
                                     ! is.na(age) ~ paste0(age, ''),
                                     ! is.na(sex) ~  substr(sex, 1, 1),
                                     T ~ '')) %>% 
    mutate(label_age_sex_date = case_when((label_age_sex != '') &  (! is.na(date)) ~ paste0(format(date, '%d/%m/%Y'), ' - ', label_age_sex, ''),
                                          (label_age_sex == '') &  (! is.na(date)) ~ paste0(format(date, '%d/%m/%Y')),
                                          label_age_sex != '' ~ paste0(label_age_sex),
                                          T ~ '')) %>% 
    mutate(age_group = cut(age, breaks = c(seq(0, 80, 10), 120), include.lowest = T, right = F))
  
  tree_for_plot <- ggtree(sub_tree)
  df_for_labels <- tree_for_plot$data %>% 
    left_join(df_aa_muts, by = c('label' = 'strain')) %>% 
    filter(! is.na(muts_for_tree))
  
  return(list(
    'sub_tree' = sub_tree,
    'metadata' = metadata_tips,
    'df_for_labels' = df_for_labels
  ))
}

list_tree_node_US_cluster_1 <- get_tree_metadata_from_node(name_nodes_to_plot[4])
list_tree_node_US_cluster_2 <- get_tree_metadata_from_node(name_nodes_to_plot[5])
list_tree_node_Angola_UK <- get_tree_metadata_from_node(name_nodes_to_plot[6]) # Small Angola / UK cluster
list_tree_node_global_cluster <- get_tree_metadata_from_node(name_nodes_to_plot[7]) # Denmark / Thailand / UK / USA cluster


color_palette_states <- c('Massachusetts' = '#591F00FF',
                          'New Hampshire' = '#A5420BFF',
                          'New Jersey' = '#EE6100FF', 
                          'New York' = '#FFAD0AFF', 
                          'Rhode Island' = '#FFCE9FFF',
                          'Maryland' = '#D72000FF', 
                          'Virginia' = '#FF7E79FF',
                          'Tennessee' = '#FFC2C1FF',
                          'Georgia' = '#4B1D3FFF',
                          'Florida' = "#956F89FF",
                          'Indiana' = '#004A49FF',
                          'Minnesota' = '#17BEBBFF',
                          'California' = '#006634FF',
                          'Arizona' = "#30985BFF",
                          'Colorado' = '#95C3A2FF',
                          'Texas' = 'gray50')


library(paletteer)
palette_africa <- paletteer_dynamic(palette = "cartography::blue.pal", n = 20)
palette_asia <- paletteer_dynamic(palette = "cartography::red.pal", n = 20)
palette_north_america <- paletteer_dynamic(palette = "cartography::sand.pal", n = 20)
palette_south_america <- paletteer_dynamic(palette = "cartography::brown.pal", n = 20)
palette_europe <- paletteer_dynamic(palette = "cartography::green.pal", n = 20)
palette_oceania <- paletteer_dynamic(palette = "cartography::purple.pal", n = 20)

rectangle_start_US_cluster_1 <- 15

color_male <- 'gray70'
color_female <- 'black'

age_levels <- c( "[0,10)", "[10,20)", "[20,30)",
                 "[30,40)", "[40,50)", "[50,60)",
                 "[60,70)", "[70,80)", "[80,120]")


tree_US_cluster_1 <- ggtree(list_tree_node_US_cluster_1$sub_tree) %<+%
  (list_tree_node_US_cluster_1$metadata %>% mutate(subdivision_1 = factor(subdivision_1, levels = names(color_palette_states)))) +
  geom_rootedge(color = col_adaptive, size = 1) +
  geom_tippoint(aes(colour = subdivision_1)) +
  geom_richtext(
    data = list_tree_node_US_cluster_1$df_for_labels,
    aes(x = x, y = y, label = muts_for_tree_colored),
    hjust = 1., vjust = 0., nudge_y = 1, nudge_x = -0.5,
    size = 3.,
    fill = colorspace::lighten(col_adaptive, 0.9), 
    label.color = col_adaptive
  ) +
  geom_tile(aes(x = rectangle_start_US_cluster_1, y = y, fill = as.Date(date)), width = 0.45, height = 0.9) +
  scale_fill_viridis_c(na.value = NA, trans = "date", name = 'Collection date') +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = rectangle_start_US_cluster_1 + 0.5, y = y, fill = sex),
            width = 0.45, height = 0.9, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Female" = color_female, "Male" = color_male),
                    breaks = c('Male', 'Female'),
                    na.value = NA, name = "Sex") +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = rectangle_start_US_cluster_1 + 1., y = y, fill = age_group),
            width = 0.45, height = 0.9, inherit.aes = FALSE) +
  scale_fill_manual(
    na.value = NA, 
    name = "Age",
    breaks = age_levels,
    values = brewer.pal(9, "RdYlBu")
  ) +
  geom_treescale(x = -10, y = 0, width = 1., offset = 1.) +
  scale_colour_manual(values = color_palette_states, name = 'State') +
  coord_fixed(ratio = 0.5, clip = "off") +
  theme(legend.background = element_blank(),
        legend.position = "right",
        legend.box = "horizontal")

rectangle_start_US_cluster_2 <- 9

tree_US_cluster_2 <- ggtree(list_tree_node_US_cluster_2$sub_tree) %<+%
  (list_tree_node_US_cluster_2$metadata %>% mutate(subdivision_1 = factor(subdivision_1, levels = names(color_palette_states)))) +
  geom_rootedge(color = col_adaptive, size = 1) +
  geom_tippoint(aes(colour = subdivision_1)) +
  geom_richtext(
    data = list_tree_node_US_cluster_2$df_for_labels,
    aes(x = x, y = y, label = muts_for_tree_colored),
    hjust = 1., vjust = 0., nudge_x = -0.3, nudge_y = 0.4,
    size = 3.,
    fill = colorspace::lighten(col_adaptive, 0.9), 
    label.color = col_adaptive,
  ) +
  geom_treescale(x = -10, y = 0, width = 1., offset = 0.1) +
  scale_colour_manual(values = color_palette_states, name = 'State') +
  geom_tile(aes(x = rectangle_start_US_cluster_2, y = y, fill = as.Date(date)), width = 0.45, height = 0.9) +
  scale_fill_viridis_c(na.value = NA, trans = "date", name = 'Collection date') +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = rectangle_start_US_cluster_2 + 0.5, y = y, fill = sex),
            width = 0.45, height = 0.9, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Female" = color_female, "Male" = color_male),
                    breaks = c('Male', 'Female'),
                    na.value = NA, name = "Sex") +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = rectangle_start_US_cluster_2 + 1., y = y, fill = age_group),
            width = 0.45, height = 0.9, inherit.aes = FALSE) +
  scale_fill_manual(
    na.value = NA, 
    name = "Age",
    breaks = age_levels,
    values = brewer.pal(9, "RdYlBu")
  ) +
  coord_fixed(ratio = 0.5, clip = "off") +
  theme(legend.background = element_blank(),
        legend.position = "right",
        legend.box = "horizontal")


rectangle_start_Angola_UK <- 5

tree_cluster_Angola_UK <- ggtree(list_tree_node_Angola_UK$sub_tree) %<+%
  (list_tree_node_Angola_UK$metadata) +
  geom_rootedge(color = col_adaptive, size = 1) +
  geom_tiplab(aes(label = label_age_sex_date)) +
  geom_tippoint(aes(colour = country)) +
  geom_richtext(
    data = list_tree_node_Angola_UK$df_for_labels,
    aes(x = x, y = y, label = muts_for_tree_colored),
    hjust = 1., vjust = 0., nudge_x = - 0.3, nudge_y = 0.1,
    size = 3.,
    fill = colorspace::lighten(col_adaptive, 0.9), 
    label.color = col_adaptive
  ) +
  geom_treescale(x = -10, y = 0, width = 1., offset = 0.1) +
  scale_colour_manual(values = c(palette_africa[10], palette_europe[15]),
                      breaks = c('Angola', 'United Kingdom'),
                      name = 'Country') +
  geom_tile(aes(x = rectangle_start_Angola_UK, y = y, fill = as.Date(date)), width = 0.45, height = 0.9) +
  scale_fill_viridis_c(na.value = NA, trans = "date", name = 'Collection date') +
  coord_fixed(ratio = 0.5, clip = "off")  +
  theme(legend.background = element_blank(),
        legend.position = "right",
        legend.box = "horizontal")


rectangle_start_global_cluster <- 11

tree_global_cluster <- ggtree(list_tree_node_global_cluster$sub_tree) %<+%
  (list_tree_node_global_cluster$metadata) +
  geom_rootedge(color = col_adaptive, size = 1) +
  geom_tippoint(aes(colour = country)) +
  geom_richtext(
    data = list_tree_node_global_cluster$df_for_labels,
    aes(x = x, y = y, label = muts_for_tree_colored),
    hjust = 1., vjust = 0., nudge_x = -0.5, nudge_y = 1.,
    size = 3.,
    fill = colorspace::lighten(col_adaptive, 0.9), 
    label.color = col_adaptive
  ) +
  geom_treescale(x = -10, y = 0, width = 1., offset = 0.1) +
  scale_colour_manual(values = c(palette_europe[5], palette_europe[15], palette_asia[10], palette_north_america[10]),
                      breaks = c('Denmark', 'United Kingdom', 'Thailand', 'USA'),
                      name = 'Country') +
  geom_tile(aes(x = rectangle_start_global_cluster, y = y, fill = as.Date(date)), width = 0.45, height = 0.9) +
  scale_fill_viridis_c(na.value = NA, trans = "date", name = 'Collection date') +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = rectangle_start_global_cluster + 0.5, y = y, fill = sex),
            width = 0.45, height = 0.9, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Female" = color_female, "Male" = color_male),
                    breaks = c('Male', 'Female'),
                    na.value = NA, name = "Sex") +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = rectangle_start_global_cluster + 1., y = y, fill = age_group),
            width = 0.45, height = 0.9, inherit.aes = FALSE) +
  scale_fill_manual(
    na.value = NA, 
    name = "Age",
    breaks = age_levels,
    values = brewer.pal(9, "RdYlBu")
  ) +
  coord_fixed(ratio = 0.5, clip = "off")  +
  theme(legend.background = element_blank(),
        legend.position = "right",
        legend.box = "horizontal")

plot(tree_global_cluster)

pdf('../figure_example_saltational_clusters/tree_cluster_US_1.pdf',
    width = 10, height = 5)
plot(tree_US_cluster_1)
dev.off()
pdf('../figure_example_saltational_clusters/tree_cluster_US_2.pdf',
    width = 10, height = 5)
plot(tree_US_cluster_2)
dev.off()
pdf('../figure_example_saltational_clusters/tree_cluster_UK_Angola.pdf',
    width = 8, height = 5)
plot(tree_cluster_Angola_UK)
dev.off()
pdf('../figure_example_saltational_clusters/tree_global_cluster.pdf',
    width = 12, height = 5)
plot(tree_global_cluster)
dev.off()
