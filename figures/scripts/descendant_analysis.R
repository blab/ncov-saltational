library(tidyverse)
source('config_plot.R')

vec_genes <- names(palette_genes)
vec_genes <- vec_genes[vec_genes != 'Genome']

# Load adaptive branches names
name_adaptive_branches <- read.table('../../results/adaptive_branches/concat_name_adaptive_branches_CV_1.7_BF_31.6227766.tsv', 
                                     header = F) %>% unlist() %>% as.character()

# Load aa muts count information
df_aa_muts_count_adaptive_branches <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  filter(is_coded, node_id %in% name_adaptive_branches) %>% 
  rename(strain = node_id) %>% 
  group_by(strain) %>% 
  summarise(n_muts = n(), n_non_syn = sum(aa_initial != aa_new),
            n_muts_S = sum(gene == 'S'), n_non_syn_S = sum((aa_initial != aa_new) & (gene == 'S')),
            n_muts_ORF1a = sum(gene == 'ORF1a'), n_non_syn_ORF1a = sum((aa_initial != aa_new) & (gene == 'ORF1a')),
            n_muts_ORF1b = sum(gene == 'ORF1b'), n_non_syn_ORF1b = sum((aa_initial != aa_new) & (gene == 'ORF1b')),
            n_muts_ORF3a = sum(gene == 'ORF3a'), n_non_syn_ORF3a = sum((aa_initial != aa_new) & (gene == 'ORF3a')),
            n_muts_ORF6 = sum(gene == 'ORF6'), n_non_syn_ORF6 = sum((aa_initial != aa_new) & (gene == 'ORF6')),
            n_muts_ORF7a = sum(gene == 'ORF7a'), n_non_syn_ORF7a = sum((aa_initial != aa_new) & (gene == 'ORF7a')),
            n_muts_ORF7b = sum(gene == 'ORF7b'), n_non_syn_ORF7b = sum((aa_initial != aa_new) & (gene == 'ORF7b')),
            n_muts_ORF8 = sum(gene == 'ORF8'), n_non_syn_ORF8 = sum((aa_initial != aa_new) & (gene == 'ORF8')),
            n_muts_ORF9b = sum(gene == 'ORF9b'), n_non_syn_ORF9b = sum((aa_initial != aa_new) & (gene == 'ORF9b')),
            n_muts_M = sum(gene == 'M'), n_non_syn_M = sum((aa_initial != aa_new) & (gene == 'M')),
            n_muts_N = sum(gene == 'N'), n_non_syn_N = sum((aa_initial != aa_new) & (gene == 'N')),
            n_muts_E = sum(gene == 'E'), n_non_syn_E = sum((aa_initial != aa_new) & (gene == 'E')),
            n_syn_S = n_muts_S - n_non_syn_S,
            n_syn_M = n_muts_M - n_non_syn_M,
            n_syn_N = n_muts_N - n_non_syn_N,
            n_syn_E = n_muts_E - n_non_syn_E,
            n_syn_ORF1a = n_muts_ORF1a - n_non_syn_ORF1a,
            n_syn_ORF1b = n_muts_ORF1b - n_non_syn_ORF1b,
            n_syn_ORF3a = n_muts_ORF3a - n_non_syn_ORF3a,
            n_syn_ORF6 = n_muts_ORF6 - n_non_syn_ORF6,
            n_syn_ORF7a = n_muts_ORF7a - n_non_syn_ORF7a,
            n_syn_ORF7b = n_muts_ORF7b - n_non_syn_ORF7b,
            n_syn_ORF8 = n_muts_ORF8 - n_non_syn_ORF8,
            n_syn_ORF9b = n_muts_ORF9b - n_non_syn_ORF9b)

# Load name of all branches
df_child_parent <- read_tsv('../../results/bte/df_parent_child_relationship.tsv') %>% 
  filter(parent_name != 'None') %>% 
  mutate(is_adaptive = strain %in% name_adaptive_branches)

# Load descendant information
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv')


col_adaptive <- 'orange2'
col_non_adaptive <- 'darkgrey'

# Proportion of branches with descendants
barplot_prop_branches_with_descendants <- df_descendants %>% 
  mutate(n_terminal_descendants = ifelse(n_terminal_descendants == 1, 0, n_terminal_descendants)) %>% 
  mutate(is_adaptive = strain %in% name_adaptive_branches) %>% 
  group_by(is_adaptive) %>% summarise(prop_with_descendants = sum(n_terminal_descendants > 0)/n()) %>% 
  ggplot(aes(x = !is_adaptive, fill = is_adaptive, y = prop_with_descendants)) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(name = '', breaks = c(F, T), labels = c('Adaptive', 'Non-adaptive')) +
  scale_y_continuous(name = 'Proportion of branches with descendants',
                     expand = expansion(mult = c(0., 0.05)))  +
  scale_fill_manual(breaks = c(T, F), values = c(col_adaptive, col_non_adaptive), guide = 'none') +
  theme_classic() +
  theme(axis.text = element_text(size = 11))

pdf('../descendant_analysis/barplot_prop_branches_with_descendants.pdf', height = 3.8, width = 2.5)
plot(barplot_prop_branches_with_descendants)
dev.off()
png('../descendant_analysis/barplot_prop_branches_with_descendants.png', height = 3.8, width = 2.5, res = 350, units = 'in')
plot(barplot_prop_branches_with_descendants)
dev.off()

# Distribution of descendants among adaptive branches with descendants
barplot_number_descendants_adaptive_branches_with_descendants <- df_descendants %>% 
  mutate(n_terminal_descendants = ifelse(n_terminal_descendants == 1, 0, n_terminal_descendants)) %>% 
  filter(n_terminal_descendants > 0) %>% 
  mutate(is_adaptive = strain %in% name_adaptive_branches) %>% 
  filter(is_adaptive) %>% 
  mutate(n_terminal_descendants_crop = ifelse(n_terminal_descendants >= 10, 10, n_terminal_descendants)) %>% 
  group_by(n_terminal_descendants_crop) %>% 
  summarise(n_adaptive_branches = n()) %>% 
  ggplot(aes(x = n_terminal_descendants_crop, y = n_adaptive_branches)) +
  geom_bar(stat = 'identity', fill = 'gray22') +
  scale_x_continuous(name = 'Number of descendants', breaks = 0:10, labels = c(0:9, '10+')) +
  scale_y_continuous(name = 'Number of adaptive branches', expand = expansion(mult = c(0., 0.05))) +
  theme_classic() +
  theme(axis.text = element_text(size = 11))

pdf('../descendant_analysis/barplot_number_descendants_adaptive_branches_with_descendants.pdf', height = 3.8, width = 5.)
plot(barplot_number_descendants_adaptive_branches_with_descendants)
dev.off()
png('../descendant_analysis/barplot_number_descendants_adaptive_branches_with_descendants.png',
    height = 3.8, width = 5., res = 350, units = 'in')
plot(barplot_number_descendants_adaptive_branches_with_descendants)
dev.off()


# Look at determinants of having descendants
df_adaptive_descendants <- df_descendants %>% 
  mutate(n_terminal_descendants = ifelse(n_terminal_descendants == 1, 0, n_terminal_descendants)) %>% 
  mutate(has_descendants = n_terminal_descendants > 0) %>% 
  filter(strain %in% name_adaptive_branches) %>% 
  left_join(df_aa_muts_count_adaptive_branches)

df_adaptive_descendants %>% 
  select(has_descendants, n_terminal_descendants, 
         n_muts_S, n_muts_ORF1a, n_muts_ORF1b, n_muts_ORF3a, n_muts_ORF6, n_muts_ORF7a, n_muts_ORF7b,
         n_muts_ORF8, n_muts_ORF9b, n_muts_N, n_muts_M, n_muts_E) %>% 
  pivot_longer(cols = 3:14, values_to = 'n_muts', names_to = 'gene', names_prefix = 'n_muts_') %>% 
  mutate(has_muts = n_muts > 0) %>% 
  group_by(gene) %>% 
  summarise(n_branches = n(), n_branches_with_muts = sum(has_muts),
            n_branches_without_descendants_with_muts = sum(! has_descendants & has_muts),
            n_branches_with_descendants_with_muts = sum(has_descendants & has_muts),
            n_branches_without_descendants_without_muts = sum(! has_descendants & !has_muts),
            n_branches_with_descendants_without_muts  = sum(has_descendants & !has_muts)) %>% 
  group_by(gene) %>% 
  summarise(prop_with_muts_with_desc = n_branches_with_descendants_with_muts/(n_branches_with_descendants_with_muts + n_branches_without_descendants_with_muts),
            prop_without_muts_with_desc = n_branches_with_descendants_without_muts/(n_branches_with_descendants_without_muts + n_branches_without_descendants_without_muts)) %>% 
  pivot_longer(cols = 2:3, values_to = 'prop') %>% 
  ggplot(aes(x = gene, y = prop, colour = name)) +
  geom_point(position = position_dodge(0.8))
  


names_col <- c(paste0('n_muts_', vec_genes), paste0('n_non_syn_', vec_genes), paste0('n_syn_', vec_genes))

df_log_reg_mut_having_descendants <- Reduce('bind_rows', lapply(names_col, FUN = function(curr_col){
  curr_log_reg <- glm(data = df_adaptive_descendants %>% mutate(has_mut = get(curr_col) > 0), formula = has_descendants ~ has_mut, family = 'binomial')
  confint_param <- confint(curr_log_reg)['has_mutTRUE', ]
  tibble(col_name = curr_col,
         estim_param = curr_log_reg$coefficients['has_mutTRUE'],
         lower_param = confint(curr_log_reg)['has_mutTRUE', 1],
         upper_param = confint(curr_log_reg)['has_mutTRUE', 2])
}))

df_log_reg_n_mut_having_descendants <- Reduce('bind_rows', lapply(names_col, FUN = function(curr_col){
  curr_log_reg <- glm(data = df_adaptive_descendants %>% mutate(n_mut = get(curr_col)), formula = has_descendants ~ n_mut, family = 'binomial')
  confint_param <- confint(curr_log_reg)['n_mut', ]
  tibble(col_name = curr_col,
         estim_param = curr_log_reg$coefficients['n_mut'],
         lower_param = confint(curr_log_reg)['n_mut', 1],
         upper_param = confint(curr_log_reg)['n_mut', 2])
}))

plt_OR_having_descendants_func_having_mut_in_gene <- df_log_reg_mut_having_descendants %>% 
  group_by(col_name) %>% 
  mutate(mut_type = case_when(startsWith(col_name, 'n_muts_') ~ 'All', 
                              startsWith(col_name, 'n_non_syn') ~ 'Non-synonymous', 
                              startsWith(col_name, 'n_syn') ~ 'Synonymous'),
         gene = case_when(startsWith(col_name, 'n_muts_') ~ strsplit(col_name, 'n_muts_')[[1]][2],
                          startsWith(col_name, 'n_non_syn_') ~ strsplit(col_name, 'n_non_syn_')[[1]][2],
                          startsWith(col_name, 'n_syn_') ~ strsplit(col_name, 'n_syn_')[[1]][2])) %>% 
  filter(mut_type == 'Non-synonymous') %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  ggplot(aes(x = gene, y = exp(estim_param), colour = gene, shape = mut_type)) +
  geom_hline(yintercept = 1.0, linetype = 'dashed', color = 'darkgrey') +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = exp(lower_param), ymax = exp(upper_param)), position = position_dodge(0.4)) +
  scale_x_discrete(name = 'Gene') +
  scale_y_continuous(name = 'OR', trans = 'log', breaks = c(0.1, 0.2, 0.5, 1.0, 2, 5)) +
  scale_colour_manual(breaks = names(palette_genes), values = palette_genes, name = '') +
  #scale_shape_manual(name = '', values = c(15, 16, 17)) +
  theme_classic() +
  theme(axis.text = element_text(size = 11))

png('../descendant_analysis/OR_having_descendants_with_mut_in_gene.png',
    height = 3.5, width = 6.5, res = 350, units = 'in')
plot(plt_OR_having_descendants_func_having_mut_in_gene + theme(legend.position = 'none'))
dev.off()
pdf('../descendant_analysis/OR_having_descendants_with_mut_in_gene.pdf',
    height = 3.5, width = 6.5)
plot(plt_OR_having_descendants_func_having_mut_in_gene + theme(legend.position = 'none'))
dev.off()

plt_OR_having_descendants_func_n_mut_in_gene <- df_log_reg_mut_having_descendants %>% 
  group_by(col_name) %>% 
  mutate(mut_type = case_when(startsWith(col_name, 'n_muts_') ~ 'All', 
                              startsWith(col_name, 'n_non_syn') ~ 'Non-synonymous', 
                              startsWith(col_name, 'n_syn') ~ 'Synonymous'),
         gene = case_when(startsWith(col_name, 'n_muts_') ~ strsplit(col_name, 'n_muts_')[[1]][2],
                          startsWith(col_name, 'n_non_syn_') ~ strsplit(col_name, 'n_non_syn_')[[1]][2],
                          startsWith(col_name, 'n_syn_') ~ strsplit(col_name, 'n_syn_')[[1]][2])) %>% 
  filter(mut_type == 'Non-synonymous') %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  ggplot(aes(x = gene, y = exp(estim_param), colour = gene, shape = mut_type)) +
  geom_hline(yintercept = 1.0, linetype = 'dashed', color = 'darkgrey') +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = exp(lower_param), ymax = exp(upper_param)), position = position_dodge(0.4)) +
  scale_x_discrete(name = 'Gene') +
  scale_y_continuous(name = 'OR', trans = 'log', breaks = c(0.1, 0.2, 0.5, 1.0, 2, 5)) +
  scale_colour_manual(breaks = names(palette_genes), values = palette_genes, name = '') +
  #scale_shape_manual(name = '', values = c(15, 16, 17)) +
  theme_classic() +
  theme(axis.text = element_text(size = 11))

png('../descendant_analysis/OR_having_descendants_with_one_additional_AA_change_in_gene.png',
    height = 3.5, width = 6.5, res = 350, units = 'in')
plot(plt_OR_having_descendants_func_n_mut_in_gene + theme(legend.position = 'none'))
dev.off()
pdf('../descendant_analysis/OR_having_descendants_with_one_additional_AA_change_in_gene.pdf',
    height = 3.5, width = 6.5)
plot(plt_OR_having_descendants_func_n_mut_in_gene + theme(legend.position = 'none'))
dev.off()




df_poisson_reg_having_mut_n_descendants <- Reduce('bind_rows', lapply(names_col, FUN = function(curr_col){
  curr_poisson_reg <- glm(data = df_adaptive_descendants %>% 
                           mutate(n_descendants = ifelse(n_terminal_descendants > 0, n_terminal_descendants - 1, n_terminal_descendants),
                                  has_mut = get(curr_col) > 0),
                         formula = n_descendants ~ has_mut, family = 'poisson')
  confint_param <- confint(curr_poisson_reg)['has_mutTRUE', ]
  
  tibble(col_name = curr_col,
         estim_param = curr_poisson_reg$coefficients['has_mutTRUE'],
         lower_param = confint(curr_poisson_reg)['has_mutTRUE', 1],
         upper_param = confint(curr_poisson_reg)['has_mutTRUE', 2])
}))

plt_param_estim_more_descendants_having_mut_in_gene <- df_poisson_reg_having_mut_n_descendants %>% 
  group_by(col_name) %>% 
  mutate(mut_type = case_when(startsWith(col_name, 'n_muts_') ~ 'All', 
                              startsWith(col_name, 'n_non_syn') ~ 'Non-synonymous', 
                              startsWith(col_name, 'n_syn') ~ 'Synonymous'),
         gene = case_when(startsWith(col_name, 'n_muts_') ~ strsplit(col_name, 'n_muts_')[[1]][2],
                          startsWith(col_name, 'n_non_syn_') ~ strsplit(col_name, 'n_non_syn_')[[1]][2],
                          startsWith(col_name, 'n_syn_') ~ strsplit(col_name, 'n_syn_')[[1]][2])) %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  ggplot(aes(x = gene, y = exp(estim_param), colour = gene, shape = mut_type)) +
  geom_hline(yintercept = 1., linetype = 'dashed', color = 'darkgrey') +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = exp(lower_param), ymax = exp(upper_param)), position = position_dodge(0.4)) +
  scale_x_discrete(name = 'Gene') +
  scale_y_continuous(name = 'Incidence rate ratio', trans = 'log', 
                     breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1e3, 1e4),
                     labels = c(expression(10^{-4}), expression(10^{-3}), expression(10^{-2}), expression(10^{-1}),
                                expression(10^{0}), expression(10^{1}), expression(10^{2}), expression(10^{3}),
                                expression(10^{4}))) +
  scale_colour_manual(breaks = names(palette_genes), values = palette_genes, name = '') +
  scale_shape_manual(name = '', values = c(15, 16, 17)) +
  theme_classic() +
  theme(axis.text = element_text(size = 11))


### Comparing Negative Binomial vs Poisson model
df_comp_poisson_negbin_having_mut_n_descendants <- Reduce('bind_rows', lapply(names_col, FUN = function(curr_col){
  curr_poisson_reg <- glm(data = df_adaptive_descendants %>% 
                            mutate(n_descendants = ifelse(n_terminal_descendants > 0, n_terminal_descendants - 1, n_terminal_descendants),
                                   has_mut = get(curr_col) > 0),
                          formula = n_descendants ~ has_mut, family = 'poisson')
  curr_negbin_reg <- MASS::glm.nb(data = df_adaptive_descendants %>% 
                                    mutate(n_descendants = ifelse(n_terminal_descendants > 0, n_terminal_descendants - 1, n_terminal_descendants),
                                           has_mut = get(curr_col) > 0),
                                  formula = n_descendants ~ has_mut)
  curr_AIC <- AIC(curr_poisson_reg, curr_negbin_reg)
  tibble(col_name = curr_col,
         AIC_poisson = curr_AIC['curr_poisson_reg', 'AIC'],
         AIC_NB = curr_AIC['curr_negbin_reg', 'AIC'])
}))

df_comp_poisson_negbin_having_mut_n_descendants %>% View()
  
## Reimplement but with a NegBin
df_NB_reg_having_mut_n_descendants <- Reduce('bind_rows', lapply(names_col, FUN = function(curr_col){
  curr_negbin_reg <- MASS::glm.nb(data = df_adaptive_descendants %>% 
                                    mutate(n_descendants = ifelse(n_terminal_descendants > 0, n_terminal_descendants - 1, n_terminal_descendants),
                                           has_mut = get(curr_col) > 0),
                                  formula = n_descendants ~ has_mut)
  confint_param <- confint(curr_negbin_reg)['has_mutTRUE', ]
  
  tibble(col_name = curr_col,
         estim_param = curr_negbin_reg$coefficients['has_mutTRUE'],
         lower_param = confint(curr_negbin_reg)['has_mutTRUE', 1],
         upper_param = confint(curr_negbin_reg)['has_mutTRUE', 2])
}))

plt_param_estim_NB_more_descendants_having_mut_in_gene <- df_NB_reg_having_mut_n_descendants %>% 
  group_by(col_name) %>% 
  mutate(mut_type = case_when(startsWith(col_name, 'n_muts_') ~ 'All', 
                              startsWith(col_name, 'n_non_syn') ~ 'Non-synonymous', 
                              startsWith(col_name, 'n_syn') ~ 'Synonymous'),
         gene = case_when(startsWith(col_name, 'n_muts_') ~ strsplit(col_name, 'n_muts_')[[1]][2],
                          startsWith(col_name, 'n_non_syn_') ~ strsplit(col_name, 'n_non_syn_')[[1]][2],
                          startsWith(col_name, 'n_syn_') ~ strsplit(col_name, 'n_syn_')[[1]][2])) %>% 
  filter(mut_type == 'Non-synonymous') %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  ggplot(aes(x = gene, y = exp(estim_param), colour = gene)) +
  geom_hline(yintercept = 1., linetype = 'dashed', color = 'darkgrey') +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = exp(lower_param), ymax = exp(upper_param)), position = position_dodge(0.4)) +
  scale_x_discrete(name = 'Gene') +
  scale_y_continuous(name = 'Fold change in number of descendants\nwith 1 additional AA change', trans = 'log', 
                     breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1e3, 1e4),
                     labels = c(expression(10^{-4}), expression(10^{-3}), expression(10^{-2}), expression(10^{-1}),
                                expression(10^{0}), expression(10^{1}), expression(10^{2}), expression(10^{3}),
                                expression(10^{4}))) +
  scale_colour_manual(breaks = names(palette_genes), values = palette_genes, name = '') +
  #scale_shape_manual(name = '', values = c(15, 16, 17)) +
  theme_classic() +
  theme(axis.text = element_text(size = 11)) 

pdf('../descendant_analysis/fold_change_nb_descendants_1_AA_change_univariate_NB_model.pdf',
    height = 3.5, width = 6.8)
plot(plt_param_estim_NB_more_descendants_having_mut_in_gene + theme(legend.position = 'none'))
dev.off()
png('../descendant_analysis/fold_change_nb_descendants_1_AA_change_univariate_NB_model.png',
    height = 3.5, width = 6.8, res = 350, units = 'in')
plot(plt_param_estim_NB_more_descendants_having_mut_in_gene + theme(legend.position = 'none'))
dev.off()

## N descendants func n mut in gene - Trying a zero-inflated model
library(pscl)
df_zeroinfl_NB_reg_n_mut_n_descendants <- Reduce('bind_rows', lapply(names_col, FUN = function(curr_col){
  curr_negbin_reg <- pscl::zeroinfl(data = df_adaptive_descendants %>% 
                                      mutate(n_descendants = ifelse(n_terminal_descendants > 0, n_terminal_descendants - 1, n_terminal_descendants),
                                             n_mut = get(curr_col)),
                                    formula = n_descendants ~ n_mut | 1, dist = 'negbin')
  
  confint_param <- confint(curr_negbin_reg)
  
  tibble(col_name = curr_col,
         estim_param = curr_negbin_reg$coefficients$count['n_mut'],
         lower_param = confint(curr_negbin_reg)['count_n_mut', 1],
         upper_param = confint(curr_negbin_reg)['count_n_mut', 2])
}))

plt_param_estim_zero_infl_NB_more_descendants_having_mut_in_gene <- df_zeroinfl_NB_reg_n_mut_n_descendants %>% 
  group_by(col_name) %>% 
  mutate(mut_type = case_when(startsWith(col_name, 'n_muts_') ~ 'All', 
                              startsWith(col_name, 'n_non_syn') ~ 'Non-synonymous', 
                              startsWith(col_name, 'n_syn') ~ 'Synonymous'),
         gene = case_when(startsWith(col_name, 'n_muts_') ~ strsplit(col_name, 'n_muts_')[[1]][2],
                          startsWith(col_name, 'n_non_syn_') ~ strsplit(col_name, 'n_non_syn_')[[1]][2],
                          startsWith(col_name, 'n_syn_') ~ strsplit(col_name, 'n_syn_')[[1]][2])) %>%  
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  filter(mut_type == 'Non-synonymous') %>% 
  ggplot(aes(x = gene, y = exp(estim_param), colour = gene)) +
  geom_hline(yintercept = 1., linetype = 'dashed', color = 'darkgrey') +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = exp(lower_param), ymax = exp(upper_param)), position = position_dodge(0.4)) +
  scale_x_discrete(name = 'Gene') +
  scale_y_continuous(name = 'Fold change in nb of descendants\n from\n1 additional AA change', trans = 'log', 
                     breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1e3, 1e4),
                     labels = c(expression(10^{-4}), expression(10^{-3}), expression(10^{-2}), expression(10^{-1}),
                                expression(10^{0}), expression(10^{1}), expression(10^{2}), expression(10^{3}),
                                expression(10^{4}))) +
  scale_colour_manual(breaks = names(palette_genes), values = palette_genes, name = '') +
  #scale_shape_manual(name = '', values = c(15, 16, 17)) +
  theme_classic() +
  theme(axis.text = element_text(size = 11)) 

pdf('../descendant_analysis/fold_change_nb_descendants_1_AA_change_univariate_zero_infl_NB_model.pdf',
    height = 3.5, width = 6.8)
plot(plt_param_estim_zero_infl_NB_more_descendants_having_mut_in_gene + theme(legend.position = 'none'))
dev.off()
png('../descendant_analysis/fold_change_nb_descendants_1_AA_change_univariate_zero_infl_NB_model.png',
    height = 3.5, width = 6.8, res = 350, units = 'in')
plot(plt_param_estim_zero_infl_NB_more_descendants_having_mut_in_gene + theme(legend.position = 'none'))
dev.off()

# Joint model on number of non-synonymous mutations
negbin_model_joint_n_non_syn_muts <- pscl::zeroinfl(
  data = df_adaptive_descendants %>% 
    mutate(n_descendants = ifelse(n_terminal_descendants > 0, n_terminal_descendants - 1, n_terminal_descendants)),
  formula = n_descendants ~ n_non_syn_S + n_non_syn_ORF1a + n_non_syn_ORF1b + n_non_syn_ORF3a + 
    n_non_syn_ORF6 + n_non_syn_ORF7a + n_non_syn_ORF7b + n_non_syn_ORF8 + n_non_syn_ORF9b + 
    n_non_syn_E + n_non_syn_M + n_non_syn_N | 1, dist = 'negbin')

plt_fold_change_n_desc_func_extra_AA_change_joint <- confint(negbin_model_joint_n_non_syn_muts)[startsWith(names(coef(negbin_model_joint_n_non_syn_muts)), 'count_n_non_syn_'), ] %>% 
  as_tibble() %>% 
  mutate(param_estim = coef(negbin_model_joint_n_non_syn_muts)[startsWith(names(coef(negbin_model_joint_n_non_syn_muts)), 'count_n_non_syn_')]) %>% 
  rename(lower_param = `2.5 %`, upper_param = `97.5 %`) %>% 
  mutate(name_col = names(coef(negbin_model_joint_n_non_syn_muts))[startsWith(names(coef(negbin_model_joint_n_non_syn_muts)), 'count_n_non_syn_')]) %>% 
  group_by(name_col) %>% 
  mutate(gene = strsplit(name_col, 'count_n_non_syn_')[[1]][2]) %>% ungroup() %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  ggplot(aes(x = gene, y = exp(param_estim), colour = gene)) +
  geom_hline(yintercept = 1., linetype = 'dashed', color = 'darkgrey') +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = exp(lower_param), ymax = exp(upper_param)), position = position_dodge(0.4)) +
  scale_x_discrete(name = 'Gene') +
  scale_y_continuous(name = 'Fold change in nb of descendants\nin adaptive branches from\n1 additional AA change', trans = 'log', 
                     breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1e3, 1e4),
                     labels = c(expression(10^{-4}), expression(10^{-3}), expression(10^{-2}), expression(10^{-1}),
                                expression(10^{0}), expression(10^{1}), expression(10^{2}), expression(10^{3}),
                                expression(10^{4}))) +
  scale_colour_manual(breaks = names(palette_genes), values = palette_genes, name = '') +
  scale_shape_manual(name = '', values = c(15, 16, 17)) +
  theme_classic() +
  theme(axis.text = element_text(size = 11)) 

pdf('../descendant_analysis/fold_change_nb_descendants_1_AA_change_multivariate_zero_infl_NB_model.pdf',
    height = 3.5, width = 6.8)
plot(plt_fold_change_n_desc_func_extra_AA_change_joint + theme(legend.position = 'none'))
dev.off()

png('../descendant_analysis/fold_change_nb_descendants_1_AA_change_multivariate_zero_infl_NB_model.png',
    height = 3.5, width = 6.8, res = 350, units = 'in')
plot(plt_fold_change_n_desc_func_extra_AA_change_joint + theme(legend.position = 'none'))
dev.off()


negbin_model_joint_n_across_muts <-pscl::zeroinfl(
  data = df_adaptive_descendants %>% 
    mutate(n_descendants = ifelse(n_terminal_descendants > 0, n_terminal_descendants - 1, n_terminal_descendants)),
  formula = n_descendants ~ n_non_syn_S + n_non_syn_ORF1a + n_non_syn_ORF1b + n_non_syn_ORF3a + 
    n_non_syn_ORF6 + n_non_syn_ORF7a + n_non_syn_ORF7b + n_non_syn_ORF8 + n_non_syn_ORF9b + 
    n_non_syn_E + n_non_syn_M + n_non_syn_N + n_syn_S + n_syn_ORF1a + n_syn_ORF1b + n_syn_ORF3a + 
    n_syn_ORF6 + n_syn_ORF7a + n_syn_ORF7b + n_syn_ORF8 + n_syn_ORF9b + 
    n_syn_E + n_syn_M + n_syn_N | 1, dist = 'negbin')


plt_fold_change_n_desc_func_extra_mut_joint <- confint(negbin_model_joint_n_across_muts)[startsWith(names(coef(negbin_model_joint_n_across_muts)), 'count_n_'), ] %>% 
  as_tibble() %>% 
  mutate(param_estim = coef(negbin_model_joint_n_across_muts)[startsWith(names(coef(negbin_model_joint_n_across_muts)), 'count_n_')]) %>% 
  rename(lower_param = `2.5 %`, upper_param = `97.5 %`) %>% 
  mutate(name_col = names(coef(negbin_model_joint_n_across_muts))[startsWith(names(coef(negbin_model_joint_n_across_muts)), 'count_n_')]) %>% 
  group_by(name_col) %>% 
  mutate(name_col = strsplit(name_col, 'count_n_')[[1]][2],
         gene = strsplit(name_col, '_')[[1]][length(strsplit(name_col, '_')[[1]])],
         mut_type = ifelse(startsWith(name_col, 'syn_'), 'Synonymous', 'Non-synonymous')) %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  ggplot(aes(x = gene, y = exp(param_estim), colour = gene, shape = mut_type)) +
  geom_hline(yintercept = 1., linetype = 'dashed', color = 'darkgrey') +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = exp(lower_param), ymax = exp(upper_param)), position = position_dodge(0.4)) +
  scale_x_discrete(name = 'Gene') +
  scale_y_continuous(name = 'Fold change in nb of descendants\nfrom 1 additional mutation', trans = 'log', 
                     breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1e3, 1e4),
                     labels = c(expression(10^{-4}), expression(10^{-3}), expression(10^{-2}), expression(10^{-1}),
                                expression(10^{0}), expression(10^{1}), expression(10^{2}), expression(10^{3}),
                                expression(10^{4}))) +
  scale_colour_manual(breaks = names(palette_genes), values = palette_genes, name = '', guide = 'none') +
  scale_shape_manual(name = '', values = c(15, 16, 17)) +
  theme_classic() +
  theme(axis.text = element_text(size = 11),
        legend.position = 'top')



pdf('../descendant_analysis/fold_change_nb_descendants_1_mut_multivariate_zero_infl_NB_model.pdf',
    height = 3.5, width = 6.8)
plot(plt_fold_change_n_desc_func_extra_mut_joint)
dev.off()

png('../descendant_analysis/fold_change_nb_descendants_1_mut_multivariate_zero_infl_NB_model.png',
    height = 3.8, width = 6.8, res = 350, units = 'in')
plot(plt_fold_change_n_desc_func_extra_mut_joint)
dev.off()

