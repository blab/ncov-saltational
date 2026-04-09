library(tidyverse)
library(binom)

# Load adaptive branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.7_BF_20.tsv', 
                                        header = F) %>% unlist() %>% as.character()

# Load inferred county / subregion / continent
df_inferred_countries <- read_tsv('../../results/bte/inferred_countries.tsv') %>% mutate(scale = 'Country')
df_inferred_subregions <- read_tsv('../../results/bte/inferred_subregions.tsv') %>% mutate(scale = 'Region')
df_inferred_continents <- read_tsv('../../results/bte/inferred_continents.tsv') %>% mutate(scale = 'Continent')

df_inferred_locations <- bind_rows(df_inferred_continents, df_inferred_subregions, df_inferred_countries)

# 
df_prop_salt_per_location <- df_inferred_locations %>% 
  filter(strain != 'node_1') %>% 
  mutate(is_saltational = strain %in% name_saltational_branches) %>% 
  group_by(inferred_location, scale, is_saltational) %>% 
  summarise(n_branches = n()) %>% 
  group_by(inferred_location, scale) %>% 
  summarise(n_saltational_branches = ifelse(sum(is_saltational) == 0, 0, n_branches[is_saltational]),
            n_tot_branches = sum(n_branches)) %>% 
  mutate(prop_saltational_branches = n_saltational_branches/n_tot_branches)


df_CI_prop_salt_branches <- binom.confint(
  x = df_prop_salt_per_location$n_saltational_branches,
  n = df_prop_salt_per_location$n_tot_branches,
  methods = "wilson"
) %>% 
  select(lower, upper) %>% 
  rename(lower_prop_salt = lower, upper_prop_salt = upper)



threshold_branches <- 5000
df_prop_salt_per_location %>% bind_cols(df_CI_prop_salt_branches) %>% 
  filter(inferred_location != 'Multiple descendant locations', n_tot_branches >= threshold_branches) %>% 
  ggplot(aes(x = inferred_location, y = prop_saltational_branches)) +
  geom_point() +
  geom_linerange(aes(ymin = lower_prop_salt, ymax = upper_prop_salt)) +
  facet_wrap(. ~ scale, scale = 'free') +
  coord_flip()





df_prop_salt_per_location_per_year %>% group_by(scale, inferred_location) %>% 
  summarise(n_saltational_branches = sum(n_saltational_branches), 
            n_tot_branches = sum(n_tot_branches)) %>% 
  mutate(prop_saltational_branches = n_saltational_branches/n_tot_branches) %>% 
  filter(scale == 'Continent', 
         inferred_location != 'Multiple descendant locations') %>% 
  ggplot(aes(x = inferred_location, y = prop_saltational_branches, colour = inferred_location, group = inferred_location)) +
  geom_point() 

df_prop_salt_per_location_per_year %>% group_by(scale, inferred_location) %>% 
  summarise(n_saltational_branches = sum(n_saltational_branches), 
            n_tot_branches = sum(n_tot_branches)) %>% 
  mutate(prop_saltational_branches = n_saltational_branches/n_tot_branches) %>% 
  filter(scale == 'Country', 
         inferred_location != 'Multiple descendant locations', n_tot_branches >= 10000) %>%
  ggplot(aes(x = inferred_location, y = prop_saltational_branches, colour = inferred_location, group = inferred_location)) +
  geom_point() +
  coord_flip()



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
