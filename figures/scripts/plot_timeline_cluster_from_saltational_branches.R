library(tidyverse)
library(viridis)
library(RColorBrewer)
library(ggpubr)
# Load name of saltational branches
vec_name_node_saltational <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv') %>% unlist() %>% as.character()

# Load number of descendants
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% vec_name_node_saltational) %>% 
  mutate(earliest_date_descendants = as.Date(earliest_date_descendants),
         latest_date_descendant = as.Date(latest_date_descendant))

# Load timeline of descendants
df_name_terminal_descendants <- read_tsv('../../results/saltational_branches/terminal_descendants_CV_1.8_BF_20.tsv')

# Load SRA metadata
metadata_tips <- read_tsv('../../results/saltational_branches/metadata_taxonium_CV_1.8_BF_20.tsv')
metadata_clean_dates <- read_tsv('../../results/sra_metadata/sra_metadata_just_clean_dates.tsv')
chronumental_dates <- read_tsv('../../results/chronumental/chronumental_inferred_dates.tsv')

df_for_plot_cluster_timeline <- df_name_terminal_descendants %>% 
  left_join(metadata_clean_dates %>% rename(date_descendant = date), by = c('strain_descendant' = 'strain')) %>% 
  left_join(chronumental_dates %>% rename(date_saltational = predicted_date), by = c('strain_saltational' = 'strain')) %>% 
  mutate(month_descendant = as.Date(paste0(format(date_descendant, '%Y-%m'), '-01'), format = '%Y-%m-%d'),
         month_saltational = as.Date(paste0(format(date_saltational, '%Y-%m'), '-01'), format = '%Y-%m-%d')) %>% 
  filter(date_saltational <= date_descendant) %>% 
  group_by(strain_saltational, month_descendant, date_saltational,month_saltational) %>% 
  summarise(n_terminal_descendants = n()) %>% 
  mutate(strain_saltational = factor(strain_saltational))
  
plt_cluster_through_time <- df_for_plot_cluster_timeline %>% 
  ggplot(aes(x = month_descendant, y = reorder(strain_saltational, date_saltational))) +
  geom_tile(aes(fill = n_terminal_descendants > 0, colour = n_terminal_descendants > 0)) +
  scale_x_date() +
  scale_y_discrete(name = 'Saltational cluster',
                   expand = expansion(mult = c(0.01, 0.01))) +
  scale_fill_manual(values = c('gray22', 'white'), breaks = c(T, F), guide = 'none') +
  scale_colour_manual(values = c('gray22', 'white'), breaks = c(T, F), guide = 'none') +
  theme_classic() +
  theme(axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white'))

label_125 <- function(x) {
  x <- ifelse(abs(x - round(x)) < 1e-5, round(x), x)
  m <- as.integer(signif(x / 10^floor(log10(x) + 1e-8), 12))
  ifelse(m %in% as.integer(c(1, 2, 5)), x, "")
}

label_1 <- function(x) {
  x <- ifelse(abs(x - round(x)) < 1e-5, round(x), x)
  m <- as.integer(signif(x / 10^floor(log10(x) + 1e-8), 12))
  ifelse(m %in% as.integer(c(1)), x, "")
}

plt_ecdf_cluster_duration <- df_name_terminal_descendants %>% 
  left_join(metadata_clean_dates %>% rename(date_descendant = date), by = c('strain_descendant' = 'strain')) %>% 
  left_join(chronumental_dates %>% rename(date_saltational = predicted_date), by = c('strain_saltational' = 'strain')) %>% 
  mutate(month_descendant = as.Date(paste0(format(date_descendant, '%Y-%m'), '-01'), format = '%Y-%m-%d'),
         month_saltational = as.Date(paste0(format(date_saltational, '%Y-%m'), '-01'), format = '%Y-%m-%d')) %>% 
  filter(date_saltational <= date_descendant) %>% 
  mutate(time_from_saltational = as.numeric(date_descendant - as.Date(date_saltational), units = 'days')) %>% 
  group_by(strain_saltational) %>% 
  summarise(max_time_from_saltational = max(time_from_saltational),
            min_time_from_saltational = min(time_from_saltational)) %>% 
  ggplot(aes(x = max_time_from_saltational - min_time_from_saltational)) +
  geom_step(aes(y = log10(1. - ..y..)), stat = 'ecdf', color = 'gray22') +
  scale_x_continuous(trans = 'log', 
                     labels = label_1,
                     breaks = c(seq(1e0, 9e0, 1e0), seq(1e1, 9e1, 1e1), seq(1e2, 9e2, 1e2), 1e3, 2e3),
                     'Cluster duration\n(days)') +
  scale_y_continuous(breaks = log10(c(seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 1e0, 1e-1))), 
                     labels = c(expression(10^{-3}), rep('', 8), expression(10^{-2}), rep('', 8), 
                                expression(10^{-1}), rep('', 8), expression(10^{0})),
                     name = 'Inverse ECDF') +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 11))
  

panel_cluster_duration <- 
  ggarrange(plt_cluster_through_time, plt_ecdf_cluster_duration+ coord_fixed(ratio = 2.),
            labels = 'AUTO', widths = c(1., 0.7))

pdf('../figure_descendant_analysis/panel_cluster_duration.pdf',
    height = 4.5, width = 7.)
plot(panel_cluster_duration)
dev.off()
png('../figure_descendant_analysis/panel_cluster_duration.png',
    height = 4.5, width = 7., res = 350, units = 'in')
plot(panel_cluster_duration)
dev.off()


df_name_terminal_descendants %>% 
  left_join(metadata_clean_dates %>% rename(date_descendant = date), by = c('strain_descendant' = 'strain')) %>% 
  left_join(chronumental_dates %>% rename(date_saltational = predicted_date), by = c('strain_saltational' = 'strain')) %>% 
  mutate(month_descendant = as.Date(paste0(format(date_descendant, '%Y-%m'), '-01'), format = '%Y-%m-%d'),
         month_saltational = as.Date(paste0(format(date_saltational, '%Y-%m'), '-01'), format = '%Y-%m-%d')) %>% 
  filter(date_saltational <= date_descendant) %>% 
  mutate(time_from_saltational = as.numeric(date_descendant - as.Date(date_saltational), units = 'days')) %>% 
  group_by(strain_saltational) %>% 
  summarise(max_time_from_saltational = max(time_from_saltational),
            min_time_from_saltational = min(time_from_saltational)) %>% 
  mutate(cluster_duration = max_time_from_saltational - min_time_from_saltational) %>% 
  summarise(median_cluster_duration = median(cluster_duration),
            percentile_90_cluster_duration = quantile(cluster_duration, 0.9),
            percentile_95_cluster_duration = quantile(cluster_duration, 0.95))
