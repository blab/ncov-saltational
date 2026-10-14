### This script reproduces figures depicting the 
### temporal distribution of saltational branches

library(tidyverse)
library(binom)

#############################
### LOADING MODEL RESULTS ###
#############################

# Load name of saltational branches
name_saltational_branches <- 
  read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
             header = F) %>% unlist() %>% as.character()

# Load list of all branches
df_child_parent <- read_tsv('../../results/bte/df_parent_child_relationship.tsv') %>% 
  filter(parent_name != 'None') %>% # Remove root node (not a branch)
  mutate(is_saltational = strain %in% name_saltational_branches)

## Total proportion of saltational events
length(name_saltational_branches)/nrow(df_child_parent) * 100

# Load SRA metadata
metadata <- read_tsv('../../results/sra_metadata/sra_metadata_just_clean_dates.tsv') %>% 
  mutate(year = as.Date(paste0(substr(date, start = 1, stop = 4), '-01-01')),
         bool_missing_year = is.na(year))
# Get names of tips with missing year
name_tips_missing_year <- metadata$strain[metadata$bool_missing_year]

# Load chronumental data for tips with missing year and internal nodes
chronumental <- read_tsv('../../results/chronumental/chronumental_inferred_dates.tsv') %>% 
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

# Reproduce Figure 2A
plot_number_branches_per_year <- function(){
  scaler_nb_branches <- 4000
  
  plt_nb_branches_per_year <- df_branches_with_dates_branch_count_per_year %>% 
    select(year_final, n_saltational_internal, n_saltational_terminal) %>% 
    pivot_longer(cols = 2:3, names_to = 'branch_type', names_prefix = 'n_saltational_',
                 values_to = 'n_saltational') %>% 
    ggplot(aes(x = year_final)) +
    geom_bar(aes(y = n_saltational, fill = branch_type), stat = 'identity') +
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
  
  return(plt_nb_branches_per_year)
}
plt_nb_branches_per_year <- plot_number_branches_per_year()

pdf('../figure_spatio_temporal_pattern/nb_branches_per_year.pdf',
    height = 2.8, width = 5.)
plot(plt_nb_branches_per_year)
dev.off()
png('../figure_spatio_temporal_pattern/nb_branches_per_year.png',
    height = 2.8, width = 5., res = 350, units = 'in')
plot(plt_nb_branches_per_year)
dev.off()

# Reproduce Figure 2B
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
    coord_cartesian(ylim = c(0, 0.0058))
  
  return(plt_prop_salt_branches)
}
plt_prop_salt_branches <- plot_prop_salt_branches_per_year()

pdf('../figure_spatio_temporal_pattern/prop_saltational_branches_per_year.pdf',
    height = 2.8, width = 4.5)
plot(plt_prop_salt_branches)
dev.off()
png('../figure_spatio_temporal_pattern/prop_saltational_branches_per_year.png',
    height = 2.8, width = 4.5, res = 350, units = 'in')
plot(plt_prop_salt_branches)
dev.off()

#################################################
### VISUALIZE BRANCHES THROUGH TIME AND SPACE ###
#################################################
# Figure not included in the manuscript
## Load mapping between country and region
df_country_mapping <- read_csv('../../params/mapping/country_to_continent_mapping.csv')

## Load inferred continent (ancestral reconstruction)
df_inferred_continents <- read_tsv('../../results/bte/inferred_continents.tsv') %>% mutate(scale = 'Continent')

## Get proportion of saltational branches per year per continent
df_branch_count_per_year_per_continent <- df_branches_with_dates %>% 
  left_join(df_inferred_continents %>% select(strain, inferred_location)) %>% 
  group_by(year_final, inferred_location) %>% 
  summarise(n_saltational = sum(is_saltational), n_branches = n(), 
            n_saltational_internal = sum(is_internal & is_saltational), n_internal = sum(is_internal),
            n_saltational_terminal = sum(!is_internal & is_saltational), n_terminal = sum(!is_internal)) %>% 
  filter(year_final < as.Date('2025-01-01')) %>% 
  mutate(prop_saltational = n_saltational/n_branches,
         prop_saltational_terminal = n_saltational_terminal/n_terminal,
         prop_saltational_internal = n_saltational_internal / n_internal) %>% 
  filter(inferred_location != 'Multiple descendant locations') # Remove branches w/o inferred continent

# Get CI around these proportions
df_CI_per_year_per_continent <- binom.confint(
  x = df_branch_count_per_year_per_continent$n_saltational,
  n = df_branch_count_per_year_per_continent$n_branches,
  methods = "wilson"
) %>% 
  select(lower, upper) %>% 
  rename(lower_prop_saltational = lower, upper_prop_saltational = upper)

## Define colors for figure
palette_africa <- paletteer_dynamic(palette = "cartography::blue.pal", n = 20)
palette_asia <- paletteer_dynamic(palette = "cartography::red.pal", n = 20)
palette_north_america <- paletteer_dynamic(palette = "cartography::sand.pal", n = 20)
palette_south_america <- paletteer_dynamic(palette = "cartography::brown.pal", n = 20)
palette_europe <- paletteer_dynamic(palette = "cartography::green.pal", n = 20)
palette_oceania <- paletteer_dynamic(palette = "cartography::purple.pal", n = 20)

get_plot_prop_branches_per_continent_per_year <- function(){
  plt_prop_branches_per_continent_per_year <- df_branch_count_per_year_per_continent %>% 
    bind_cols(df_CI_per_year_per_continent) %>% 
    ggplot(aes(x = year_final, y = prop_saltational, colour = inferred_location,
               group = inferred_location)) +
    geom_point(position = position_dodge(0.5)) +
    geom_linerange(aes(ymin = lower_prop_saltational, 
                       ymax = upper_prop_saltational),
                   position = position_dodge(0.5)) +
    scale_x_date(name = '') +
    scale_y_continuous(name = 'Proportion of branches\nthat are saltational', labels = scales::percent) +
    scale_colour_manual(breaks = c('North America', 'South America', 'Europe', 'Africa', 'Asia', 'Oceania'),
                        values = c(palette_north_america[10], palette_south_america[10], 
                                   palette_europe[10], palette_africa[10], 
                                   palette_asia[10], palette_oceania[10]),
                        guide = 'none') +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    coord_cartesian(ylim = c(0, 0.01)) +
    facet_wrap(. ~ inferred_location)
  
  return(plt_prop_branches_per_continent_per_year)
}

plt_prop_branches_per_continent_per_year <- get_plot_prop_branches_per_continent_per_year()

# pdf('../figure_spatio_temporal_pattern/prop_saltational_branches_per_year_per_continent.pdf',
#     height = 4., width = 7.5)
# plot(plt_prop_branches_per_continent_per_year)
# dev.off()
# png('../figure_spatio_temporal_pattern/prop_saltational_branches_per_year_per_continent.png',
#     height = 4., width = 7.5, res = 350, units = 'in')
# plot(plt_prop_branches_per_continent_per_year)
# dev.off()
