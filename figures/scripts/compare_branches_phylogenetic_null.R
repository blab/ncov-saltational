rm(list = ls())

### Compare the saltational branches detected in the observed tree with those detected
### in replicates of the null, where each branch keeps its number of mutations in each
### gene and only the synonymous / non-synonymous split is redrawn based on each gene's
### tree-wide dN/dS.

library(tidyverse)
library(ggpubr)

source('config_plot.R')

col_observed <- col_adaptive
col_null <- 'steelblue4'
vec_types <- c('Observed', 'Null')
vec_cols <- setNames(c(col_observed, col_null), vec_types)

dir_branch_counts <- '../../results/null_model/branch_counts_CV_1.8_BF_20'
vec_replicates <- 1:20

#############################
### LOADING MODEL RESULTS ###
#############################

# Saltational branches of the observed tree, as defined by the pipeline
name_saltational_branches <-
  read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv',
             header = F) %>% unlist() %>% as.character()

# Mutation counts per branch and gene, for the observed tree and for each null replicate.
# The observed file is produced by the null model itself, so both arms are counted
# identically; the observed saltational branches are taken from the pipeline.
read_branch_counts <- function(file_name){
  read_tsv(file.path(dir_branch_counts, file_name), show_col_types = FALSE) %>%
    mutate(n_muts_total = n_nonsyn_total + n_syn_total,
           n_muts_S = n_nonsyn_S + n_syn_S)
}

df_observed <- read_branch_counts('branch_counts_observed.tsv.gz') %>%
  filter(strain %in% name_saltational_branches) %>%
  mutate(type = 'Observed', replicate = NA_integer_)

df_null <- Reduce('bind_rows', lapply(vec_replicates, FUN = function(i_replicate){
  read_branch_counts(sprintf('branch_counts_replicate_%d.tsv.gz', i_replicate)) %>%
    filter(is_saltational) %>%
    mutate(type = 'Null', replicate = i_replicate)
}))

df_branches <- bind_rows(df_observed, df_null) %>%
  mutate(type = factor(type, levels = vec_types))

##############################################
### (i) NUMBER OF SALTATIONAL BRANCHES     ###
##############################################
df_n_branches <- df_branches %>%
  group_by(type, replicate) %>%
  summarise(n_branches = n(), .groups = 'drop')

n_observed <- df_n_branches$n_branches[df_n_branches$type == 'Observed']

plt_n_saltational <- df_n_branches %>%
  filter(type == 'Null') |> 
  ggplot(aes(x = 0, y = n_branches)) +
  stat_summary(fun.data = function(x) {
    r <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    return(r)
  }, 
  geom = "boxplot", 
  fill = col_null, alpha = 0.6, colour = col_null) +
  geom_point(aes(y = n_observed, x = 0), colour = col_observed, shape = 18, size = 4) +
  scale_y_continuous(name = 'Number of saltational branches') +
  theme_classic() +
  theme(axis.text = element_text(size = 11),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

##############################################
### (ii) MUTATIONS ON SALTATIONAL BRANCHES ###
##############################################
vec_metrics <- c('n_muts_total' = 'All mutations',
                 'n_nonsyn_total' = 'Nonsyn. mutations',
                 'n_muts_S' = 'Spike mutations ',
                 'n_nonsyn_S' = 'Spike nonsyn. mutations')

df_long <- df_branches %>%
  select(type, replicate, all_of(names(vec_metrics))) %>%
  pivot_longer(cols = all_of(names(vec_metrics)), names_to = 'metric', values_to = 'n_muts') %>%
  mutate(metric = factor(vec_metrics[metric], levels = vec_metrics))

# Proportion of saltational branches carrying at least x mutations (1 - ECDF), evaluated
# on a common grid so that the spread across null replicates can be shown as an interval.
vec_grid <- 0:60

df_survival <- df_long %>%
  group_by(type, replicate, metric) %>%
  reframe(n_muts_grid = vec_grid,
          prop_at_least = sapply(vec_grid, FUN = function(k) mean(n_muts >= k))) %>%
  rename(n_muts = n_muts_grid) %>%
  group_by(type, metric, n_muts) %>%
  summarise(prop_mean = mean(prop_at_least),
            prop_lower = min(prop_at_least),
            prop_upper = max(prop_at_least), .groups = 'drop')

plt_distribution <- df_survival %>%
  filter(prop_mean > 0) %>%
  ggplot(aes(x = n_muts, colour = type, fill = type)) +
  geom_ribbon(aes(ymin = prop_lower, ymax = prop_upper), alpha = 0.3, colour = NA) +
  geom_line(aes(y = prop_mean), linewidth = 0.7) +
  facet_wrap(. ~ metric, nrow = 2) +
  scale_x_continuous(name = 'Number of mutations on a saltational branch (n)') +
  scale_y_continuous(name = 'Proportion of saltational branches\nwith at least n mutations',
                     trans = 'log10',
                     labels = scales::label_number(drop0trailing = TRUE)) +
  scale_colour_manual(name = '', values = vec_cols) +
  scale_fill_manual(name = '', values = vec_cols) +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.position = 'bottom')

panel_compare_null_saltational <- 
  ggarrange(ggarrange(ggplot() + theme_void(),
                      plt_n_saltational, 
                      ggplot() + theme_void(),
                      nrow = 3, heights = c(0.1, 1., .3)),
            ggplot() + theme_void(),
            plt_distribution, ncol = 3, widths = c(0.3, 0.04, 1.),
            labels = c('A','B'))

png('../figure_model_illustration/compare_model_null.png',
    height = 4.5, width = 6.5, res = 350, units = 'in')
plot(panel_compare_null_saltational)
dev.off()
