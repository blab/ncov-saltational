library(tidyverse)
library(binom)

# Load adaptive branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.7_BF_20.tsv', 
                                        header = F) %>% unlist() %>% as.character()

# Load name of all branches
df_child_parent <- read_tsv('../../results/bte/df_parent_child_relationship.tsv') %>% 
  filter(parent_name != 'None') %>% 
  mutate(is_adaptive = strain %in% name_adaptive_branches)

# Load SRA metadata
metadata <- read_tsv('../../results/sra_metadata/sra_metadata_just_clean_dates.tsv') %>% 
  mutate(year = as.Date(paste0(substr(date, start = 1, stop = 4), '-01-01')),
         bool_missing_year = is.na(year))
name_tips_missing_year <- metadata$strain[metadata$bool_missing_year]

# Load chronumental data for tips with missing year and internal nodes
chronumental <- read_tsv('../../results/chronumental/chronumental_inferred_dates.tsv') %>% 
  filter(startsWith(strain, 'node_') | strain %in% name_tips_missing_year) %>% 
  mutate(date_chronumental = as.Date(predicted_date),
         year_chronumental = as.Date(paste0(substr(date_chronumental, start = 1, stop = 4), '-01-01')))

# Use SRA dates if available and chronumental dates for internal nodes and tips with missing dates
df_branches_with_dates <- df_child_parent %>% 
  left_join(metadata %>% select(strain, year, bool_missing_year), by = 'strain') %>% 
  left_join(chronumental %>% select(strain, date_chronumental, year_chronumental), by = 'strain') %>% 
  mutate(bool_missing_year = replace_na(bool_missing_year, T),
         year_final = as.Date(ifelse(bool_missing_year, as.character(year_chronumental), as.character(year)))) %>% 
  mutate(is_adaptive = strain %in% name_saltational_branches,
         is_internal = startsWith(strain, 'node_'))

# 
df_branches_with_dates_branch_count_per_year <- df_branches_with_dates %>% 
  group_by(year_final) %>% 
  summarise(n_adaptive = sum(is_adaptive), n_branches = n(), 
            n_adaptive_internal = sum(is_internal & is_adaptive), n_internal = sum(is_internal),
            n_adaptive_terminal = sum(!is_internal & is_adaptive), n_terminal = sum(!is_internal)) %>% 
  filter(year_final < as.Date('2025-01-01')) %>% 
  mutate(prop_adaptive = n_adaptive/n_branches,
         prop_adaptive_terminal = n_adaptive_terminal/n_terminal,
         prop_adaptive_internal = n_adaptive_internal / n_internal)
  
# Get CIs around proportion for plot
df_CI_prop_adaptive_branches <- binom.confint(
  x = df_branches_with_dates_branch_count_per_year$n_adaptive,
  n = df_branches_with_dates_branch_count_per_year$n_branches,
  methods = "wilson"
) %>% 
  select(lower, upper) %>% 
  rename(lower_prop_adaptive = lower, upper_prop_adaptive = upper)

df_CI_prop_adaptive_branches_internal <- binom.confint(
  x = df_branches_with_dates_branch_count_per_year$n_adaptive_internal,
  n = df_branches_with_dates_branch_count_per_year$n_internal,
  methods = "wilson"
) %>% 
  select(lower, upper) %>% 
  rename(lower_prop_adaptive_internal  = lower, upper_prop_adaptive_internal = upper)

df_CI_prop_adaptive_branches_terminal <- binom.confint(
  x = df_branches_with_dates_branch_count_per_year$n_adaptive_terminal,
  n = df_branches_with_dates_branch_count_per_year$n_terminal,
  methods = "wilson"
) %>% 
  select(lower, upper) %>% 
  rename(lower_prop_adaptive_terminal  = lower, upper_prop_adaptive_terminal = upper)

df_for_CI <- df_branches_with_dates_branch_count_per_year %>%
  bind_cols(df_CI_prop_adaptive_branches) %>% 
  bind_cols(df_CI_prop_adaptive_branches_internal) %>%
  bind_cols(df_CI_prop_adaptive_branches_terminal)

# Display results
col_internal <- '#f7b274'
col_terminal <- '#a7bafa'
col_all <- 'gray35'

scaler_nb_branches <- 5000
plt_nb_branches_per_year <- df_branches_with_dates_branch_count_per_year %>% 
  select(year_final, n_adaptive_internal, n_adaptive_terminal) %>% 
  pivot_longer(cols = 2:3, names_to = 'branch_type', names_prefix = 'n_adaptive_',
               values_to = 'n_adaptive') %>% 
  ggplot(aes(x = year_final)) +
  geom_bar(aes(y = n_adaptive, fill = branch_type), stat = 'identity') +
  geom_point(data = df_for_CI, 
             aes(x = year_final, y = n_branches / scaler_nb_branches), 
             col = 'gray20', shape = 18, size = 3) +
  scale_x_date(name = '') +
  scale_y_continuous(name = 'Number of saltational events', breaks = seq(0, 1000, 100),
                     expand = expansion(mult = c(0.01, 0.05)),
                     sec.axis = sec_axis(~ . * scaler_nb_branches,
                                         name = 'Total number of branches',
                                         breaks = seq(0, 3e6, 1e6), 
                                         labels = c('0', '1M', '2M', '3M'))) +
  scale_fill_manual(breaks = c('internal', 'terminal', 'all'),
                    values = colorspace::lighten(c(col_internal, col_terminal, col_all), 0.),
                    labels = c('Internal', 'Terminal', 'All'), 
                    name = 'Branch type') +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

plt_prop_salt_branches <- df_for_CI %>% 
  ggplot() +
  geom_point(data = df_for_CI, aes(x = year_final - 40, y = prop_adaptive, color = 'all')) +
  geom_linerange(data = df_for_CI, 
                 aes(x = year_final - 40, ymin = lower_prop_adaptive, ymax = upper_prop_adaptive,
                     color = 'all')) +
  geom_point(data = df_for_CI,
             aes(x = year_final, 
                 y = prop_adaptive_terminal, color = 'terminal')) +
  geom_linerange(data = df_for_CI, 
                 aes(x = year_final, ymin = lower_prop_adaptive_terminal,
                     ymax = upper_prop_adaptive_terminal,
                     color = 'terminal')) +
  geom_point(data = df_for_CI,
             aes(x = year_final + 40,
                 y = prop_adaptive_internal, 
                 color = 'internal')) +
  geom_linerange(data = df_for_CI, 
                 aes(x = year_final + 40, 
                     ymin = lower_prop_adaptive_internal,
                     ymax = upper_prop_adaptive_internal,
                     color = 'internal')) +
  scale_colour_manual(breaks = c('internal', 'terminal', 'all'),
                    values = colorspace::darken(c(col_internal, col_terminal, col_all), 0.),
                    labels = c('Internal', 'Terminal', 'All'), 
                    name = 'Branch type') +
  scale_x_date(name = '') +
  scale_y_continuous(name = 'Proportion of saltational branches', labels = scales::percent) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.0025))

## Save results
plot(plt_nb_branches_per_year)
plot(plt_prop_salt_branches)

pdf('../figure_spatio_temporal_pattern/nb_branches_per_year.pdf',
    height = 2.8, width = 5.)
plot(plt_nb_branches_per_year)
dev.off()
pdf('../figure_spatio_temporal_pattern/prop_saltational_branches_per_year.pdf',
    height = 2.8, width = 4.5)
plot(plt_prop_salt_branches)
dev.off()
