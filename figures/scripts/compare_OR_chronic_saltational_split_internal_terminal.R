library(tidyverse)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(binom)
source('utils_plot_chronic.R')
source('config_plot.R')

name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv',
                                        header = F) %>% unlist() %>% as.character()

# Load mutations counts from chronic infections vs background sequences
df_OR_chronic_saltational <- compute_OR_muts_chronic_infections_saltational('../../data/frequent_mutations/ryan_epci_hqcs_AA_ratio.csv',
                                                                            '../../results/mutations/nodes_aa_muts_summary_coded.tsv',
                                                                            name_saltational_branches)
df_OR_chronic_saltational_internal <- compute_OR_muts_chronic_infections_saltational_internal('../../data/frequent_mutations/ryan_epci_hqcs_AA_ratio.csv',
                                                                                              '../../results/mutations/nodes_aa_muts_summary_coded.tsv',
                                                                                              name_saltational_branches)
df_OR_chronic_saltational_terminal <- compute_OR_muts_chronic_infections_saltational_terminal('../../data/frequent_mutations/ryan_epci_hqcs_AA_ratio.csv',
                                                                                              '../../results/mutations/nodes_aa_muts_summary_coded.tsv',
                                                                                              name_saltational_branches)

plt_OR_salt_func_OR_chronic <- plot_OR_chronic_func_OR_saltational(df_OR_chronic_saltational)
plt_OR_salt_func_OR_chronic_internal <- plot_OR_chronic_func_OR_saltational(df_OR_chronic_saltational_internal)
plt_OR_salt_func_OR_chronic_terminal <- plot_OR_chronic_func_OR_saltational(df_OR_chronic_saltational_terminal)

panel_comp_OR_split_branch_type <- 
  ggarrange(plt_OR_salt_func_OR_chronic + ggtitle('All saltational branches'),
            plt_OR_salt_func_OR_chronic_internal + ggtitle('Internal saltational branches'),
            plt_OR_salt_func_OR_chronic_terminal + ggtitle('Terminal saltational branches'),
            ncol = 3, labels = 'AUTO',
            common.legend = T,
            legend = 'bottom')

pdf('../figure_mutation_enrichment/panel_comp_OR_chronic_salt_by_branch_type.pdf',
    height = 4., width = 10.)
plot(panel_comp_OR_split_branch_type)
dev.off()

png('../figure_mutation_enrichment/panel_comp_OR_chronic_salt_by_branch_type.png',
    height = 4., width = 10., res = 350, units = 'in')
plot(panel_comp_OR_split_branch_type)
dev.off()
