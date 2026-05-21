library(tidyverse)
library(ggpubr)
source('config_plot.R')

df_OR_AA_mut <- read_tsv('../../results/mutation_enrichment/df_AA_mut_enrichment_saltational_CV_1.8_BF_20.tsv') %>% 
  select(gene, mutation, threshold_for_inclusion, p_value_BH, OR, OR_lower, OR_upper)
df_OR_AA_mut_background <- read_tsv('../../results/mutation_enrichment/df_AA_mut_enrichment_saltational_opportunity_CV_1.8_BF_20.tsv') %>% 
  select(gene, mutation, threshold_for_inclusion, p_value_BH, OR, OR_lower, OR_upper) %>% 
  rename(p_value_BH_account_background = p_value_BH,
         OR_account_background = OR,
         OR_lower_account_background = OR_lower,
         OR_upper_account_background = OR_upper)


plt_OR_comparison_mut_opportunity <- df_OR_AA_mut %>% 
  left_join(df_OR_AA_mut_background, by = c('gene', 'mutation', 'threshold_for_inclusion')) %>% 
  filter(threshold_for_inclusion == 1) %>% 
  mutate(gene = factor(gene, levels = names(palette_genes))) %>% 
  ggplot(aes(x = OR, y = OR_account_background)) +
  geom_point(aes(colour = gene), alpha = 0.2) +
  stat_cor(method = 'pearson') +
  scale_x_log10(trans = 'log10',  labels = scales::label_math(10^.x, format = log10),
                expand = expansion(mult = c(0.05, 0.05)), name = 'OR not accounting for mutation opportuntiy') +
  scale_y_log10(trans = 'log10',  labels = scales::label_math(10^.x, format = log10), 
                expand = expansion(mult = c(0.05, 0.05)), name = 'OR accounting for mutation opportuntiy') +
  scale_colour_manual(values = palette_genes, name = 'Gene') +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1.))) +
  theme(axis.text = element_text(size = 11)) +
  coord_fixed()

png('../figure_mutation_enrichment/plt_compare_OR_mut_accounting_opportunity.png',
    height = 3.6, width = 5., res = 350, units = 'in')
plot(plt_OR_comparison_mut_opportunity)
dev.off()
pdf('../figure_mutation_enrichment/plt_compare_OR_mut_accounting_opportunity.pdf',
    height = 3.6, width = 5.)
plot(plt_OR_comparison_mut_opportunity)
dev.off()
