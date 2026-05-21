### This script reproduces figures depicting the 
### geographic distribution of saltational branches

library(tidyverse)
library(binom)
library(paletteer)
library(ggpubr)
library(ggrepel)

#############################
### LOADING MODEL RESULTS ###
#############################

# Load name of saltational branches
name_saltational_branches <- 
  read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
             header = F) %>% unlist() %>% as.character()

# Load country mapping info
df_country_mapping <- read_csv('../../params/mapping/country_to_continent_mapping.csv') %>% 
  select(country, continent, subregion) %>% distinct()

# Load inferred county / subregion / continent and merge into one dataframe
df_inferred_countries <- read_tsv('../../results/bte/inferred_countries.tsv') %>% mutate(scale = 'Country')
df_inferred_subregions <- read_tsv('../../results/bte/inferred_subregions.tsv') %>% mutate(scale = 'Region')
df_inferred_continents <- read_tsv('../../results/bte/inferred_continents.tsv') %>% mutate(scale = 'Continent')

df_inferred_locations <- bind_rows(df_inferred_continents, df_inferred_subregions, df_inferred_countries)

# Get proportion of saltational branches by location
df_prop_salt_per_location <- df_inferred_locations %>% 
  filter(strain != 'node_1') %>% 
  mutate(is_saltational = strain %in% name_saltational_branches) %>% 
  group_by(inferred_location, scale, is_saltational) %>% 
  summarise(n_branches = n()) %>% 
  group_by(inferred_location, scale) %>% 
  summarise(n_saltational_branches = ifelse(sum(is_saltational) == 0, 0, n_branches[is_saltational]),
            n_tot_branches = sum(n_branches)) %>% 
  mutate(prop_saltational_branches = n_saltational_branches/n_tot_branches)

# Get uncertainty around these proportions
df_CI_prop_salt_branches <- binom.confint(
  x = df_prop_salt_per_location$n_saltational_branches,
  n = df_prop_salt_per_location$n_tot_branches,
  methods = "wilson"
) %>% 
  select(lower, upper) %>% 
  rename(lower_prop_salt = lower, upper_prop_salt = upper)

########################################
### VISUALIZE BRANCHES THROUGH SPACE ###
########################################
# Define colors
palette_africa <- paletteer_dynamic(palette = "cartography::blue.pal", n = 20)
palette_asia <- paletteer_dynamic(palette = "cartography::red.pal", n = 20)
palette_north_america <- paletteer_dynamic(palette = "cartography::sand.pal", n = 20)
palette_south_america <- paletteer_dynamic(palette = "cartography::brown.pal", n = 20)
palette_europe <- paletteer_dynamic(palette = "cartography::green.pal", n = 20)
palette_oceania <- paletteer_dynamic(palette = "cartography::purple.pal", n = 20)

# Define threshold used to display results
# (only geographies with at least threshold_branches in the tree are plotted)
threshold_branches <- 5000

# Reproduce Figure 2C (by continent)
get_prop_branches_per_continent <- function(){
  plt_continents <- df_prop_salt_per_location %>% bind_cols(df_CI_prop_salt_branches) %>% 
    filter(inferred_location != 'Multiple descendant locations', n_tot_branches >= threshold_branches,
           scale == 'Continent') %>%
    ggplot(aes(x = reorder(inferred_location, prop_saltational_branches), y = prop_saltational_branches, colour = inferred_location)) +
    geom_point() +
    geom_linerange(aes(ymin = lower_prop_salt, ymax = upper_prop_salt)) +
    scale_y_continuous(name = 'Proportion of branches\nthat are saltational',
                       labels = scales::percent) +
    scale_x_discrete(name = '') +
    scale_colour_manual(breaks = c('North America', 'South America', 'Europe', 'Africa', 'Asia', 'Oceania'),
                        values = c(palette_north_america[10], palette_south_america[10], 
                                   palette_europe[10], palette_africa[10], 
                                   palette_asia[10], palette_oceania[10]),
                        guide = 'none') +
    coord_flip() +
    theme_bw()
  
  return(plt_continents)
}
plt_continents <- get_prop_branches_per_continent()

# Reproduce Figure 2D (by subregion)
get_prop_branches_per_subregion <- function(){
  plt_regions <- df_prop_salt_per_location %>% bind_cols(df_CI_prop_salt_branches) %>% 
    filter(inferred_location != 'Multiple descendant locations', n_tot_branches >= threshold_branches,
           scale == 'Region') %>% 
    left_join(df_country_mapping, by = c('inferred_location' = 'subregion')) %>% 
    ggplot(aes(x = reorder(inferred_location, prop_saltational_branches), y = prop_saltational_branches, colour = continent)) +
    geom_point() +
    geom_linerange(aes(ymin = lower_prop_salt, ymax = upper_prop_salt)) +
    scale_y_continuous(name = 'Proportion of branches\nthat are saltational',
                       labels = scales::percent) +
    scale_x_discrete(name = '') +
    scale_colour_manual(breaks = c('North America', 'South America', 'Europe', 'Africa', 'Asia', 'Oceania'),
                        values = c(palette_north_america[10], palette_south_america[10], 
                                   palette_europe[10], palette_africa[10], 
                                   palette_asia[10], palette_oceania[10]),
                        guide = 'none') +
    coord_flip() +
    theme_bw()
  
  return(plt_regions)
}
plt_regions <- get_prop_branches_per_subregion()

# Reproduce Figure 2E (by country)
get_prop_branches_per_country <- function(){
  plt_countries <- df_prop_salt_per_location %>% bind_cols(df_CI_prop_salt_branches) %>% 
    filter(inferred_location != 'Multiple descendant locations', 
           n_tot_branches >= threshold_branches,
           scale == 'Country') %>%
  left_join(df_country_mapping, by = c('inferred_location' = 'country')) %>% 
    mutate(inferred_location = ifelse(inferred_location == 'United Kingdom', 'UK', inferred_location)) %>% 
    ggplot(aes(x = reorder(inferred_location, prop_saltational_branches), y = prop_saltational_branches, colour = continent)) +
    geom_point() +
    geom_linerange(aes(ymin = lower_prop_salt, ymax = upper_prop_salt)) +
    scale_y_continuous(name = 'Proportion of branches\nthat are saltational',
                       labels = scales::percent) +
    scale_x_discrete(name = '') +
    scale_colour_manual(breaks = c('North America', 'South America', 'Europe', 'Africa', 'Asia', 'Oceania'),
                        values = c(palette_north_america[10], palette_south_america[10], 
                                   palette_europe[10], palette_africa[10], 
                                   palette_asia[10], palette_oceania[10]),
                        guide = 'none') +
    coord_flip() +
    theme_bw()
  
  return(plt_countries)
}
plt_countries <- get_prop_branches_per_country()


# Save Figures 2C-E
panel_prop_saltational_branches_by_location <- 
  ggarrange(plt_continents + ggtitle('Continent'),
            plt_regions + ggtitle('Region'),
            plt_countries + ggtitle('Country'), 
            ncol = 3)

pdf('../figure_spatio_temporal_pattern/prop_saltational_branches_by_location.pdf',
    height = 3., width = 9.)
plot(panel_prop_saltational_branches_by_location)
dev.off()
png('../figure_spatio_temporal_pattern/prop_saltational_branches_by_location.png',
    height = 3., width = 9., res = 350, units = 'in')
plot(panel_prop_saltational_branches_by_location)
dev.off()
     
#######################################
### COMPARISON WITH HIV PREVALENCE  ###
### (Ecological comparison)
#######################################
# Get countries to be included in the plot (need to meet the threshold_branches threshold)
vec_countries_to_keep <- df_prop_salt_per_location %>% bind_cols(df_CI_prop_salt_branches) %>% 
  filter(inferred_location != 'Multiple descendant locations', n_tot_branches >= threshold_branches,
         scale == 'Country') %>% 
  left_join(df_country_mapping, by = c('inferred_location' = 'country')) %>% 
  select(inferred_location) %>% unlist() %>% as.character() %>% unique()


# Load population size by country (Worldbank)
countries_demographics <- read_csv('../../data/country_demographics.csv')

# Load yearly HIV deaths (estimated from IHME)
countries_ihme_deaths <- read_csv('../../data/ihme_burden_country.csv') %>% 
  filter(cause == 'HIV/AIDS') %>% 
  group_by(country) %>% summarise(n_deaths = sum(n_deaths)) %>% 
  rename(n_deaths_hiv = n_deaths) %>% 
  mutate(country = ifelse(country == 'United Kingdom', 'UK', country),
         country = ifelse(country == 'United States of America', 'USA', country))


df_n_samples_in_tree_per_countries <- df_inferred_countries %>% 
  filter(!startsWith(strain, 'node_')) %>% 
  filter(inferred_location %in% vec_countries_to_keep) %>% 
  mutate(inferred_location = ifelse(inferred_location == 'United Kingdom', 'UK', inferred_location)) %>%
  group_by(inferred_location) %>% summarise(n_tips = n())

# Create dataframe used for comparison
df_for_ecological_plot <- countries_demographics %>% 
  left_join(countries_ihme_deaths) %>% 
  left_join(df_prop_salt_per_location %>% bind_cols(df_CI_prop_salt_branches) %>% 
              filter(inferred_location != 'Multiple descendant locations', n_tot_branches >= threshold_branches,
                     scale == 'Country') %>% 
              left_join(df_country_mapping, by = c('inferred_location' = 'country')) %>% 
              mutate(inferred_location = ifelse(inferred_location == 'United Kingdom', 'UK', inferred_location)) %>% 
              select(inferred_location, n_saltational_branches, n_tot_branches, prop_saltational_branches, lower_prop_salt, upper_prop_salt, continent),
            by = c('country' = 'inferred_location'))


# Get correlation between proportion of saltational branches and HIV deaths
df_for_ecological_plot %>% 
  summarise(cor_spearman = cor(n_deaths_hiv/population_size_2022_worldbank * 1000, 
                               (n_saltational_branches)/n_tot_branches, 
                               method = 'spearman'),
            cor_pearson = cor(n_deaths_hiv/population_size_2022_worldbank * 1000, 
                              (n_saltational_branches)/n_tot_branches, 
                              method = 'pearson'),
            cor_log_spearman = cor(log(n_deaths_hiv/population_size_2022_worldbank * 1000), 
                                   log((n_saltational_branches + 1)/n_tot_branches), 
                                   method = 'spearman'),
            cor_log_pearson = cor(log(n_deaths_hiv/population_size_2022_worldbank * 1000), 
                                  log((n_saltational_branches + 1)/n_tot_branches), 
                                  method = 'pearson'),
  )

## Visualize correlation
get_plot_cor_HIV_deaths_prop_salt_branches <- function(){
  label_125_percent <- function(x) {
    m <- signif(x / 10^floor(log10(x)), 12)   # mantissa
    labs <- scales::label_percent(accuracy = 1e-2)(x)
    ifelse(m %in% c(1, 2, 5), labs, "")
  }
  
  correlation_prop_saltational_branches_hiv_deaths <- df_for_ecological_plot %>% 
    ggplot(aes(x = n_deaths_hiv/population_size_2022_worldbank * 1000, 
               y = (n_saltational_branches + 1)/n_tot_branches)) +
    geom_point(aes(colour = continent)) +
    geom_text_repel(aes(label = country, colour = continent), max.overlaps = Inf) +
    scale_colour_manual(breaks = c('North America', 'South America', 'Europe', 'Africa', 'Asia', 'Oceania'),
                        values = c(palette_north_america[10], palette_south_america[10], 
                                   palette_europe[10], palette_africa[10], 
                                   palette_asia[10], palette_oceania[10]),
                        guide = 'none') +
    scale_x_continuous(trans = 'log', name = 'Yearly HIV deaths (per 1,000)',
                       breaks = c(seq(1e-4, 9e-4, 1e-4), seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 9e-1, 1e-1),
                                  seq(1e0, 9e0, 1e0)),
                       labels = c(expression(10^{-4}), rep('', 8),
                                  expression(10^{-3}), rep('', 8),
                                  expression(10^{-2}), rep('', 8),
                                  expression(10^{-1}), rep('', 8),
                                  expression(10^{0}), rep('', 8))) + 
    scale_y_continuous(trans = 'log', name = 'Proportion of branches that are saltational',
                       breaks = c(seq(1e-4, 9e-4, 1e-4), seq(1e-3, 9e-3, 1e-3)),
                       labels = label_125_percent) +
    stat_cor(method = 'pearson', label.x = -2., label.y = -9) +
    stat_cor(method = 'spearman', label.x = -2., label.y = -9.3, cor.coef.name = 'rho') +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size = 11)) 
  
  return(correlation_prop_saltational_branches_hiv_deaths)
}
correlation_prop_saltational_branches_hiv_deaths <- get_plot_cor_HIV_deaths_prop_salt_branches()

pdf('../figure_spatio_temporal_pattern/cor_prop_salt_branches_hiv_deaths.pdf',
    height = 4.5, width = 5.)
plot(correlation_prop_saltational_branches_hiv_deaths)
dev.off()
png('../figure_spatio_temporal_pattern/cor_prop_salt_branches_hiv_deaths.png',
    height = 4.5, width = 5., res = 350, units = 'in')
plot(correlation_prop_saltational_branches_hiv_deaths)
dev.off()

## Visualize correlation - Sensitivity analysis without South Africa
get_plot_cor_HIV_deaths_prop_salt_branches_no_south_africa <- function(){
  label_125_percent <- function(x) {
    m <- signif(x / 10^floor(log10(x)), 12)   # mantissa
    labs <- scales::label_percent(accuracy = 1e-2)(x)
    ifelse(m %in% c(1, 2, 5), labs, "")
  }
  
  correlation_prop_saltational_branches_hiv_deaths <- df_for_ecological_plot %>% 
    filter(country != 'South Africa') %>% 
    ggplot(aes(x = n_deaths_hiv/population_size_2022_worldbank * 1000, 
               y = (n_saltational_branches + 1)/n_tot_branches)) +
    geom_point(aes(colour = continent)) +
    geom_text_repel(aes(label = country, colour = continent), max.overlaps = Inf) +
    scale_colour_manual(breaks = c('North America', 'South America', 'Europe', 'Africa', 'Asia', 'Oceania'),
                        values = c(palette_north_america[10], palette_south_america[10], 
                                   palette_europe[10], palette_africa[10], 
                                   palette_asia[10], palette_oceania[10]),
                        guide = 'none') +
    scale_x_continuous(trans = 'log', name = 'Yearly HIV deaths (per 1,000)',
                       breaks = c(seq(1e-4, 9e-4, 1e-4), seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 9e-1, 1e-1),
                                  seq(1e0, 9e0, 1e0)),
                       labels = c(expression(10^{-4}), rep('', 8),
                                  expression(10^{-3}), rep('', 8),
                                  expression(10^{-2}), rep('', 8),
                                  expression(10^{-1}), rep('', 8),
                                  expression(10^{0}), rep('', 8))) + 
    scale_y_continuous(trans = 'log', name = 'Proportion of branches that are saltational',
                       breaks = c(seq(1e-4, 9e-4, 1e-4), seq(1e-3, 9e-3, 1e-3)),
                       labels = label_125_percent) +
    stat_cor(method = 'pearson', label.x = -5., label.y = -9) +
    stat_cor(method = 'spearman', label.x = -5., label.y = -9.3, cor.coef.name = 'rho') +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size = 11)) 
  
  return(correlation_prop_saltational_branches_hiv_deaths)
}
correlation_prop_saltational_branches_hiv_deaths_no_south_africa <- get_plot_cor_HIV_deaths_prop_salt_branches_no_south_africa()

pdf('../figure_spatio_temporal_pattern/cor_prop_salt_branches_hiv_deaths_no_SA.pdf',
    height = 4.5, width = 5.)
plot(correlation_prop_saltational_branches_hiv_deaths_no_south_africa)
dev.off()
png('../figure_spatio_temporal_pattern/cor_prop_salt_branches_hiv_deaths_no_SA.png',
    height = 4.5, width = 5., res = 350, units = 'in')
plot(correlation_prop_saltational_branches_hiv_deaths_no_south_africa)
dev.off()

# 
# correlation_prop_saltational_branches_hiv_prev <- df_for_ecological_plot %>% 
#   ggplot(aes(x = hiv_prevalence_percent_ihme_2022 / 100, 
#              y = (n_saltational_branches + 1)/n_tot_branches)) +
#   geom_point(aes(colour = continent)) +
#   geom_text_repel(aes(label = country, colour = continent), max.overlaps = Inf) +
#   scale_colour_manual(breaks = c('North America', 'South America', 'Europe', 'Africa', 'Asia', 'Oceania'),
#                       values = c(palette_north_america[10], palette_south_america[10], 
#                                  palette_europe[10], palette_africa[10], 
#                                  palette_asia[10], palette_oceania[10]),
#                       guide = 'none') +
#   scale_x_continuous(trans = 'log', name = 'HIV prevalence',
#                      breaks = c(seq(1e-4, 9e-4, 1e-4), seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 9e-1, 1e-1),
#                                 seq(1e0, 9e0, 1e0)),
#                      labels = label_125_percent) + 
#   scale_y_continuous(trans = 'log', name = 'Proportion of branches that are saltational',
#                      breaks = c(seq(1e-4, 9e-4, 1e-4), seq(1e-3, 9e-3, 1e-3)),
#                      labels = label_125_percent) +
#   stat_cor(method = 'pearson', label.x = -4., label.y = -9) +
#   stat_cor(method = 'spearman', label.x = -4., label.y = -9.3, cor.coef.name = 'rho') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 11)) 
# 
# pdf('../figure_spatio_temporal_pattern/cor_prop_salt_branches_hiv_prev.pdf',
#     height = 4.5, width = 6.)
# plot(correlation_prop_saltational_branches_hiv_prev)
# dev.off()
# png('../figure_spatio_temporal_pattern/cor_prop_salt_branches_hiv_prev.png',
#     height = 4.5, width = 6., res = 350, units = 'in')
# plot(correlation_prop_saltational_branches_hiv_prev)
# dev.off()
