library(tidyverse)
library(glmmTMB)
library(scales)
library(MASS)
source('utils_descendants_analysis.R')
source('config_plot.R')

# Load saltational branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', header = F) %>%
  unlist() %>% as.character()

# Load AA mutations count information
df_aa_muts_count_adaptive_branches <- get_mut_counts_per_gene_adapt_branches(
  '../../results/mutations/nodes_aa_muts_summary_coded.tsv', 
  name_saltational_branches)

# Load descendant information
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches)

# Create dataframe with number of 
df_saltational_descendants <- df_descendants %>% 
  mutate(n_terminal_descendants = ifelse(n_terminal_descendants == 1, 0, n_terminal_descendants)) %>% 
  mutate(has_descendants = n_terminal_descendants > 0) %>% 
  left_join(df_aa_muts_count_adaptive_branches) %>% 
  mutate(n_terminal_descendants_counts = ifelse(n_terminal_descendants > 0, n_terminal_descendants - 1, n_terminal_descendants))
  
#############
# Assess whether Poisson, Negative Binomial
# or Zero-inflated negative binomial
# provides a lower AIC
#############
vec_genes <- names(palette_genes)[-length(names(palette_genes))]
name_cols_syn_non_syn <- c(paste0('n_non_syn_', vec_genes), paste0('n_syn_', vec_genes))
AIC_poisson <- get_AIC_poisson_reg(df_saltational_descendants, 
                                   name_cols_syn_non_syn,
                                   threshold_n_descendants_inclusion = 5e4)
AIC_negbin <- get_AIC_negbin_reg(df_saltational_descendants, 
                                 name_cols_syn_non_syn,
                                 threshold_n_descendants_inclusion = 5e4)

likelihood_ratio_test_poisson_negbin(df_saltational_descendants, 
                                     name_cols_syn_non_syn,
                                     threshold_n_descendants_inclusion = 5e4)

## Run the negative binomial model across thresholds
vec_thresholds_n_descendants <- c(1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, Inf)
df_negbin_param <- Reduce('bind_rows', lapply(vec_thresholds_n_descendants, FUN = function(curr_threshold){
  run_negbin_reg(df_saltational_descendants, 
                 name_cols_syn_non_syn,
                 threshold_n_descendants_inclusion = curr_threshold) %>% 
    mutate(threshold_n_descendants_inclusion = curr_threshold)
}))

## Sensitivity analysis - Impact of imposing a threshold on the number of descendants over model parameters
plot_sensitivity_threshold_param_estim <- function(df_negbin_param){
  plt_sensitivity_estimates_large_branches <- df_negbin_param %>% 
    mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
    ggplot(aes(x = as.factor(threshold_n_descendants_inclusion), 
               y = exp(param_estim), group = interaction(gene),
               colour = gene)) +
    geom_hline(yintercept = 1., linetype = 'dashed', color = 'darkgrey') + 
    geom_point() +
    geom_linerange(aes(ymin = exp(lower_param), ymax = exp(upper_param))) +
    scale_x_discrete(name = 'Threshold for branch exclusion', 
                     limits = rev,
                     breaks = c(1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, Inf), 
                     labels = c('1,000', '5,000','10,000', '50,000', '100,000', '500,000', '1,000,000', 'No threshold')
    ) +
    scale_y_continuous(minor_breaks =
                         minor_breaks_log(smallest = 1, detail = 1),
                       breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4),
                       labels = c(expression(10^{-4}), expression(10^{-3}), expression(10^{-2}),
                                  expression(10^{-1}), expression(10^{0}), expression(10^{1}),
                                  expression(10^{2}), expression(10^{3}), expression(10^{4})),
                       trans = 'log',
                       name = 'Fold change in number of descendants\nfrom one additional mutation') +
    facet_wrap(gene ~ mut_type, ncol = 4) +
    scale_colour_manual(values = palette_genes, guide = 'none')  +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1., vjust = 0.5), ,
          strip.background = element_rect(fill = 'white'))
  
  return(plt_sensitivity_estimates_large_branches)
}

plt_sensitivity_threshold <- plot_sensitivity_threshold_param_estim(df_negbin_param)

png('../figure_descendant_analysis/sensitivity_fold_change_threshold.png',
    height = 8., width = 6., res = 350, units = 'in')
plot(plt_sensitivity_threshold)
dev.off()
pdf('../figure_descendant_analysis/sensitivity_fold_change_threshold.pdf',
    height = 8., width = 6.)
plot(plt_sensitivity_threshold)
dev.off()


## Display parameter estimates
plot_negbbin <- function(df_negbin, threshold_to_plot){
  label_1 <- function(x) {
    m <- signif(x / 10^floor(log10(x)), 12)
    ifelse(m %in% c(1), x, "")
  }
  
  plot_zero_infl_negbin <- df_negbin %>% 
    filter(abs(threshold_n_descendants_inclusion - threshold_to_plot) < 1e-6) %>% 
    mutate(gene = factor(gene, levels = names(palette_genes)),
           is_signif = (lower_param > 0) | (upper_param < 0)) %>%
    ggplot(aes(x = gene, y = exp(param_estim), shape = mut_type, colour = gene)) +
    geom_hline(yintercept = 1., linetype = 'dashed', color = 'gray50') +
    geom_point(position = position_dodge(0.5)) +
    geom_linerange(aes(ymin = exp(lower_param), ymax = exp(upper_param)),
                   position = position_dodge(0.5))  +
    scale_x_discrete(name = 'Gene', limits = rev) +
    scale_y_continuous(minor_breaks = minor_breaks_log(smallest = 1, detail = 1),
                       breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4),
                       labels = c(expression(10^{-4}), expression(10^{-3}), expression(10^{-2}),
                                  expression(10^{-1}), expression(10^{0}), expression(10^{1}),
                                  expression(10^{2}), expression(10^{3}), expression(10^{4})),
                       trans = 'log',
                       name = 'Fold change in number of descendants\nfrom one additional mutation') +
    scale_colour_manual(values = palette_genes, guide = 'none') +
    scale_shape_manual(name = '', values = c(16, 17)) +
    scale_alpha_manual(name = '', breaks = c(T, F), values = c(1.0, 0.5)) +
    coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_text(size = 11))
  
  return(plot_zero_infl_negbin)
}

plt_param_estim_negbin <- plot_negbbin(df_negbin_param, 5e4)

pdf('../figure_descendant_analysis/fold_change_nb_descendants_1_extra_mutation_negbin.pdf',
    height = 4., width = 5.)
plot(plt_param_estim_negbin)
dev.off()
png('../figure_descendant_analysis/fold_change_nb_descendants_1_extra_mutation_negbin.png',
    height = 4., width = 5., res = 350, units = 'in')
plot(plt_param_estim_negbin)
dev.off()

## Display parameter estimates as a function of mutation frequency
plt_fold_change_func_prop_branches_with_mut <- function(df_negbin_param, 
                                                        df_saltational_descendants, 
                                                        threshold_to_plot){
  scatter_plot_fold_change <- df_saltational_descendants %>% 
    dplyr::select(starts_with("n_non_syn_"), starts_with("n_syn_")) %>% 
    summarise(across(everything(), mean)) %>% 
    pivot_longer(cols = everything(), names_to = 'name_col', values_to = 'mean_nb_muts') %>% 
    group_by(name_col) %>% 
    mutate(gene = strsplit(name_col, '_')[[1]][length(strsplit(name_col, '_')[[1]])],
           mut_type = ifelse(startsWith(name_col, 'n_syn'), 'Synonymous', 'Nonsynonymous')) %>% 
    ungroup() %>% 
    dplyr::select(gene, mut_type, mean_nb_muts) %>% 
    left_join(df_negbin_param, by = c('gene', 'mut_type')) %>% 
    mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
    filter(threshold_n_descendants_inclusion == threshold_to_plot) %>% 
    ggplot(aes(x = mean_nb_muts, y = exp(param_estim), colour = gene, shape = mut_type)) +
    geom_hline(yintercept = 1., linetype = 'dashed') +
    geom_point() +
    geom_linerange(aes(ymin = exp(lower_param), ymax = exp(upper_param))) +
    scale_x_continuous(minor_breaks =
                         minor_breaks_log(smallest = 1, detail = 1),
                       breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4),
                       labels = c(expression(10^{-4}), expression(10^{-3}), expression(10^{-2}),
                                  expression(10^{-1}), expression(10^{0}), expression(10^{1}),
                                  expression(10^{2}), expression(10^{3}), expression(10^{4})),
                       trans = 'log',
                       name = 'Mean number of mutations\nper saltational branch') +
    scale_y_continuous(minor_breaks =
                         minor_breaks_log(smallest = 1, detail = 1),
                       breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4),
                       labels = c(expression(10^{-4}), expression(10^{-3}), expression(10^{-2}),
                                  expression(10^{-1}), expression(10^{0}), expression(10^{1}),
                                  expression(10^{2}), expression(10^{3}), expression(10^{4})),
                       trans = 'log',
                       limits = c(NA, 10),
                       expand = expansion(mult = c(0.01, 0.)),
                       name = 'Fold change in number of\ndescendants from one additional mutation') +
    scale_colour_manual(values = palette_genes, name = 'Gene') +
    scale_shape_manual(name = '', values = c(16, 17), breaks = c('Nonsynonymous', 'Synonymous'),
                       labels = c('Nonsynonymous', 'Synonymous')) +
    theme_bw() +
    coord_fixed() +
    theme(axis.text = element_text(size = 11)) +
    guides(shape = guide_legend(order = 2, ncol = 1),
           colour = guide_legend(order = 1, ncol = 1)) +
    theme(legend.position = 'bottom', legend.direction = 'horizontal') +
    guides(color = guide_legend(ncol = 4))
  
  return(scatter_plot_fold_change)
}

plt_fold_change_func_prop_salt_branches <- plt_fold_change_func_prop_branches_with_mut(df_negbin_param, 
                                                                                       df_saltational_descendants, 
                                                                                       5e4)

pdf('../figure_descendant_analysis/fold_change_func_prop_salt_branches_with_mut.pdf',
    height = 5., width = 5.)
plot(plt_fold_change_func_prop_salt_branches)
dev.off()
png('../figure_descendant_analysis/fold_change_func_prop_salt_branches_with_mut.png',
    height = 5., width = 7.5, res = 350, units = 'in')
plot(plt_fold_change_func_prop_salt_branches)
dev.off()

# Get values for text
df_negbin_param %>% 
  filter(gene %in% c('S', 'ORF1a', 'ORF1b'), 
         mut_type == 'Nonsynonymous', 
         threshold_n_descendants_inclusion == 5e4) %>% 
  mutate(lower_param = exp(lower_param),
         upper_param = exp(upper_param),
         param_estim = exp(param_estim))
