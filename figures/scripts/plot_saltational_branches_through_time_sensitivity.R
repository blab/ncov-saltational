### This script reproduces figures exploring the
### impact of sample size on the number of saltational branches
### per year

library(tidyverse)
library(binom)

#############################
### LOADING MODEL RESULTS ###
#############################

# Load name of saltational branches
name_saltational_branches <- 
  read.table('../../results/sensitivity_analysis/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
             header = F) %>% unlist() %>% as.character()


# Load list of all branches
df_child_parent <- read_tsv('../../results/sensitivity_analysis/df_parent_child_relationship.tsv') %>% 
  filter(parent_name != 'None') %>% # Remove root node (not a branch)
  mutate(is_saltational = strain %in% name_saltational_branches)


# Load SRA metadata
metadata <- read_tsv('../../results/sra_metadata/sra_metadata_just_clean_dates.tsv') %>% 
  mutate(year = as.Date(paste0(substr(date, start = 1, stop = 4), '-01-01')),
         bool_missing_year = is.na(year)) %>% 
  filter(strain %in% df_child_parent$strain)

# Get names of tips with missing year
name_tips_missing_year <- metadata$strain[metadata$bool_missing_year]

# Load chronumental data for tips with missing year and internal nodes
chronumental <- read_tsv('../../results/sensitivity_analysis/chronumental/chronumental_inferred_dates.tsv') %>% 
  filter(startsWith(strain, 'node_') | strain %in% name_tips_missing_year) %>% 
  mutate(date_chronumental = as.Date(predicted_date),
         year_chronumental = as.Date(paste0(substr(date_chronumental, start = 1, stop = 4), '-01-01')))

######################################
### DETERMINE THE DATE OF BRANCHES ###
######################################
## We define the time of a branches as the time of sequence collection 
## defining the end of the branch if the branch is terminal and the
## sequence collection date is known. Otherwise, use Chronumental 
## inferred dates. 
df_branches_with_dates <- df_child_parent %>% 
  left_join(metadata %>% select(strain, year, bool_missing_year), by = 'strain') %>% 
  left_join(chronumental %>% select(strain, date_chronumental, year_chronumental), by = 'strain') %>% 
  mutate(bool_missing_year = replace_na(bool_missing_year, T),
         year_final = as.Date(ifelse(bool_missing_year, as.character(year_chronumental), as.character(year)))) %>% 
  mutate(is_saltational = strain %in% name_saltational_branches,
         is_internal = startsWith(strain, 'node_'))

# Summarize branch count per year per type
df_branches_with_dates_branch_count_per_year <- df_branches_with_dates %>% 
  group_by(year_final) %>% 
  summarise(n_saltational = sum(is_saltational), n_branches = n(), 
            n_saltational_internal = sum(is_internal & is_saltational), n_internal = sum(is_internal),
            n_saltational_terminal = sum(!is_internal & is_saltational), n_terminal = sum(!is_internal)) %>% 
  filter(year_final < as.Date('2025-01-01')) %>% # Remove branches that were inferred after 2025 (all )
  mutate(prop_saltational = n_saltational/n_branches,
         prop_saltational_terminal = n_saltational_terminal/n_terminal,
         prop_saltational_internal = n_saltational_internal / n_internal)

# Compute Wilson 95% confidence intervals aroud the proportion of saltational branches
df_CI_prop_saltational_branches <- binom.confint(
  x = df_branches_with_dates_branch_count_per_year$n_saltational,
  n = df_branches_with_dates_branch_count_per_year$n_branches,
  methods = "wilson"
) %>% 
  select(lower, upper) %>% 
  rename(lower_prop_saltational = lower, upper_prop_saltational = upper)

df_CI_prop_saltational_branches_internal <- binom.confint(
  x = df_branches_with_dates_branch_count_per_year$n_saltational_internal,
  n = df_branches_with_dates_branch_count_per_year$n_internal,
  methods = "wilson"
) %>% 
  select(lower, upper) %>% 
  rename(lower_prop_saltational_internal  = lower, upper_prop_saltational_internal = upper)

df_CI_prop_saltational_branches_terminal <- binom.confint(
  x = df_branches_with_dates_branch_count_per_year$n_saltational_terminal,
  n = df_branches_with_dates_branch_count_per_year$n_terminal,
  methods = "wilson"
) %>% 
  select(lower, upper) %>% 
  rename(lower_prop_saltational_terminal  = lower, upper_prop_saltational_terminal = upper)

df_for_CI <- df_branches_with_dates_branch_count_per_year %>%
  bind_cols(df_CI_prop_saltational_branches) %>% 
  bind_cols(df_CI_prop_saltational_branches_internal) %>%
  bind_cols(df_CI_prop_saltational_branches_terminal)

#######################################
### VISUALIZE BRANCHES THROUGH TIME ###
#######################################
# Color used in the plot
col_internal <- '#f7b274'
col_terminal <- '#a7bafa'
col_all <- 'gray35'

plot_prop_salt_branches_per_year <- function(){
  plt_prop_salt_branches <- df_for_CI %>% 
    ggplot() +
    geom_point(data = df_for_CI, aes(x = year_final - 40, y = prop_saltational, color = 'all')) +
    geom_linerange(data = df_for_CI, 
                   aes(x = year_final - 40, ymin = lower_prop_saltational, ymax = upper_prop_saltational,
                       color = 'all')) +
    geom_point(data = df_for_CI,
               aes(x = year_final, 
                   y = prop_saltational_terminal, color = 'terminal')) +
    geom_linerange(data = df_for_CI, 
                   aes(x = year_final, ymin = lower_prop_saltational_terminal,
                       ymax = upper_prop_saltational_terminal,
                       color = 'terminal')) +
    geom_point(data = df_for_CI,
               aes(x = year_final + 40,
                   y = prop_saltational_internal, 
                   color = 'internal')) +
    geom_linerange(data = df_for_CI, 
                   aes(x = year_final + 40, 
                       ymin = lower_prop_saltational_internal,
                       ymax = upper_prop_saltational_internal,
                       color = 'internal')) +
    scale_colour_manual(breaks = c('internal', 'terminal', 'all'),
                        values = colorspace::darken(c(col_internal, col_terminal, col_all), 0.),
                        labels = c('Internal', 'Terminal', 'All'), 
                        name = 'Branch type') +
    scale_x_date(name = '') +
    scale_y_continuous(name = 'Proportion of branches\nthat are saltational', labels = scales::percent) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    coord_cartesian(ylim = c(0, 0.0035))
  
  return(plt_prop_salt_branches)
}
plt_prop_salt_branches <- plot_prop_salt_branches_per_year()

png('../figure_spatio_temporal_pattern/prop_saltational_branches_per_year_downsample_constant_by_year.png',
    height = 2.8, width = 4.5, res = 350, units = 'in')
plot(plt_prop_salt_branches)
dev.off()
pdf('../figure_spatio_temporal_pattern/prop_saltational_branches_per_year_downsample_constant_by_year.pdf',
    height = 2.8, width = 4.5)
plot(plt_prop_salt_branches)
dev.off()
