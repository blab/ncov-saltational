library(tidyverse)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(binom)
source('utils_plot_chronic.R')
source('config_plot.R')

# Load recurrent mutations in deers / minks / persistent infections
frequent_mutations_summary <- get_frequent_mutations_summary('../../data/frequent_mutations/frequent_mutations_minks_deers.csv',
                                                             '../../data/frequent_mutations/frequent_mutations_chronic.csv')

# Load saltational branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
                                        header = F) %>% unlist() %>% as.character()

# Load OR of AA mut occurring in saltational branches
df_OR_AA_mut <- get_OR_AA_mut('../../results/mutation_enrichment/df_AA_mut_enrichment_saltational_opportunity_CV_1.8_BF_20.tsv',
                              2,
                              frequent_mutations_summary)

# Load mutations counts from chronic infections vs background sequences
df_OR_chronic_saltational <- compute_OR_muts_chronic_infections_saltational('../../data/frequent_mutations/ryan_epci_hqcs_AA_ratio.csv',
                                                                            '../../results/mutations/nodes_aa_muts_summary_coded.tsv',
                                                                            name_saltational_branches)

# Load proportion of 
df_inferred_countries <- read_tsv('../../results/bte/inferred_countries.tsv') %>% mutate(scale = 'Country')
df_prop_saltational_branches_per_country <- get_prop_salt_branches_per_location(name_saltational_branches, df_inferred_countries)

df_prop_salt_branches_HIV_burden_per_country <- get_prop_salt_branches_HIV_burden_per_location(name_saltational_branches, 
                                                                                               df_inferred_countries, 
                                                                                               5000, 
                                                                                               '../../data/country_demographics.csv',
                                                                                               '../../data/ihme_burden_country.csv')

#############
# Figure 4A #
# OR as a function of mutation frequency coloured by recurrent mutation profile
#############
plt_OR_func_mut_frequency_coloured_by_profile_type <- plot_OR_func_nb_mut_colour_recurrent_mut(df_OR_AA_mut)

#############
# Figure 4B: 
# OR of mutation occurring in saltational branches as a
# function of OR of mutation occurring in chronic infections
#############
plt_OR_salt_func_OR_chronic <- plot_OR_chronic_func_OR_saltational(df_OR_chronic_saltational)

#############
# Figure 4C #
# Proportion of branches that are saltational by country as a function 
# of HIV burden in that country
#############
plt_prop_salt_branches_func_HIV_burden <- plot_prop_salt_branches_func_HIV_burden(df_prop_salt_branches_HIV_burden_per_country)

full_panel_chronic <- (plt_OR_func_mut_frequency_coloured_by_profile_type + plt_OR_salt_func_OR_chronic + plot_layout(widths = c(0.9, 1.0))) /
  (plot_spacer() + plt_prop_salt_branches_func_HIV_burden + plot_spacer() + plot_layout(widths = c(0.5, 1., 0.5))) + 
  plot_layout(heights = c(1., 1.)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = 'bold', size = 20))


# Save results
pdf('../../figures/figure_chronic_infections.pdf', height = 7.5, width = 10)
plot(full_panel_chronic)
dev.off()

png('../../figures/figure_chronic_infections.png', height = 7.5, width = 10, res = 350, units = 'in')
plot(full_panel_chronic)
dev.off()


######## Get Spearman correlation coefficients
cor.test(df_prop_salt_branches_HIV_burden_per_country$modif_saltational_branches, 
         df_prop_salt_branches_HIV_burden_per_country$hiv_deaths_per_1000,
         method = 'spearman')
cor.test(df_prop_salt_branches_HIV_burden_per_country$prop_saltational_branches[df_prop_salt_branches_HIV_burden_per_country$inferred_location != 'South Africa'], 
         df_prop_salt_branches_HIV_burden_per_country$hiv_deaths_per_1000[df_prop_salt_branches_HIV_burden_per_country$inferred_location != 'South Africa'],
         method = 'spearman')

df_OR_AA_mut %>% 
  filter(is_signif) %>% 
  summarise(n_signif_deers = sum(is_within_profile_deers),
            n_signif_minks = sum(is_within_profile_minks),
            n_signif_persistent = sum(is_within_profile_persistent))
