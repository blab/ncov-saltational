library(tidyverse)

source('config_plot.R')

df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  filter(is_coded) %>% rename(strain = node_id) %>% 
  group_by(strain) %>% 
  summarise(n_non_syn = sum(aa_initial != aa_new),
            n_syn = sum(aa_initial == aa_new))

df_non_syn_count_by_branch <- read_tsv('../../results/bte/df_parent_child_relationship.tsv') %>% 
  filter(strain != 'node_1') %>% 
  left_join(df_aa_muts) %>% 
  mutate(n_non_syn = replace_na(n_non_syn, 0),
         n_syn = replace_na(n_syn, 0)) %>% 
  select(strain, n_non_syn, n_syn)

ecdf_mut_counts_by_branch <- df_non_syn_count_by_branch %>% 
  ggplot() +
  geom_step(aes(x = n_non_syn, y = after_stat(log10(1. - y)),
                color = 'Nonsynonymous'), stat = 'ecdf') +
  geom_step(aes(x = n_syn, y = after_stat(log10(1. - y)),
                color = 'Synonymous'), stat = 'ecdf') +
  scale_x_continuous(name = 'Mutation count',
                     breaks = seq(0, 65, 5)) +
  scale_y_continuous(name = 'Complementary empirical\ncumulative distribution function',
                     breaks = seq(0, -7, -1),
                     labels = c(expression(10^{0}), expression(10^{-1}), 
                                expression(10^{-2}), expression(10^{-3}), 
                                expression(10^{-4}), expression(10^{-5}), 
                                expression(10^{-6}), expression(10^{-7})))  +
  scale_colour_manual(name = 'Mutation type',
                      breaks = c('Nonsynonymous', 'Synonymous'),
                      values = c(col_non_synonymous, col_synonymous)) +
  theme_classic() +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.7, 0.8),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 11))

png('../figure_model_illustration/ecdf_mut_counts_by_branch.png',
    height = 3., width = 5., res = 350, units = 'in')
plot(ecdf_mut_counts_by_branch)
dev.off()
