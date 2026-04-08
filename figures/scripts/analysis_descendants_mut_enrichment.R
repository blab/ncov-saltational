library(tidyverse)
library(pscl)
library(scales)
source('utils_descendants_analysis.R')
source('config_plot.R')

# Load saltational branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.7_BF_20.tsv', header = F) %>%
  unlist() %>% as.character()

# Load AA mutations count information
df_aa_muts_count_adaptive_branches <- get_mut_counts_per_gene_adapt_branches(
  '../../results/mutations/nodes_aa_muts_summary_coded.tsv', 
  name_saltational_branches)

# Load descendant information
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv') %>% 
  filter(strain %in% name_saltational_branches)

# Investigate determinants of having descendants
df_saltational_descendants <- df_descendants %>% 
  mutate(n_terminal_descendants = ifelse(n_terminal_descendants == 1, 0, n_terminal_descendants)) %>% 
  mutate(has_descendants = n_terminal_descendants > 0) %>% 
  left_join(df_aa_muts_count_adaptive_branches) %>% 
  mutate(n_terminal_descendants_counts = ifelse(n_terminal_descendants > 0, n_terminal_descendants - 1, n_terminal_descendants))
  

### LOGISTIC MODEL

# Univariate logistic regression having descendants func having mutations
vec_genes <- c('ORF1a', 'ORF1b', 'S', 'ORF3a', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'ORF9b', 'N', 'M', 'E')
names_col <- c(paste0('n_muts_', vec_genes), paste0('n_non_syn_', vec_genes), paste0('n_syn_', vec_genes))

df_log_reg_having_mut_having_descendants <- Reduce('bind_rows', lapply(names_col, FUN = function(curr_col){
  curr_log_reg <- glm(data = df_saltational_descendants %>% mutate(has_mut = get(curr_col) > 0),
                      formula = has_descendants ~ has_mut, 
                      family = 'binomial')
  confint_param <- confint(curr_log_reg)['has_mutTRUE', ]
  tibble(col_name = curr_col,
         estim_param = curr_log_reg$coefficients['has_mutTRUE'],
         lower_param = confint(curr_log_reg)['has_mutTRUE', 1],
         upper_param = confint(curr_log_reg)['has_mutTRUE', 2])
}))

# Univariate logistic regression having descendants func number of mutations
df_log_reg_n_mut_having_descendants <- Reduce('bind_rows', lapply(names_col, FUN = function(curr_col){
  curr_log_reg <- glm(data = df_saltational_descendants %>% mutate(n_mut = get(curr_col)), formula = has_descendants ~ n_mut, family = 'binomial')
  confint_param <- confint(curr_log_reg)['n_mut', ]
  tibble(col_name = curr_col,
         estim_param = curr_log_reg$coefficients['n_mut'],
         lower_param = confint(curr_log_reg)['n_mut', 1],
         upper_param = confint(curr_log_reg)['n_mut', 2])
}))

# Multivariate logistic regression having descendants nb mutations
multi_log_reg_having_descendants <- glm(has_descendants ~ .,
  data = bind_cols(df_saltational_descendants %>% 
                     dplyr::select(has_descendants, starts_with("n_non_syn_"), starts_with("n_syn_"))),
  family = 'binomial'
)

summary(multi_log_reg_having_descendants)
confint(multi_log_reg_having_descendants)

### COUNT MODEL
name_cols_syn_non_syn <- c(paste0('n_non_syn_', vec_genes), paste0('n_syn_', vec_genes))
vec_thresholds_n_descendants <- c(1e5)
df_negbin_zero_inflated_joint <- Reduce('bind_rows', lapply(vec_thresholds_n_descendants, FUN = function(curr_threshold){
  run_zero_infl_negbin_reg(df_saltational_descendants, 
                           name_cols_syn_non_syn,
                           threshold_n_descendants_inclusion = curr_threshold) %>% 
    mutate(threshold_n_descendants_inclusion = curr_threshold)
}))


plot_negbbin_zero_infl <- function(df_negbin_zero_inflated_joint){
  label_1 <- function(x) {
    m <- signif(x / 10^floor(log10(x)), 12)
    ifelse(m %in% c(1), x, "")
  }
  
  plot_zero_infl_negbin <- df_negbin_zero_inflated_joint %>% 
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

plt_param_estim_negbin <- plot_negbbin_zero_infl(df_negbin_zero_inflated_joint)

pdf('../figure_descendant_analysis/fold_change_nb_descendants_1_extra_mutation_zero_infl_negbin_joint.pdf',
    height = 4., width = 5.)
plot(plt_param_estim_negbin)
dev.off()
png('../figure_descendant_analysis/fold_change_nb_descendants_1_extra_mutation_zero_infl_negbin_joint.png',
    height = 4., width = 5., res = 350, units = 'in')
plot(plt_param_estim_negbin)
dev.off()
