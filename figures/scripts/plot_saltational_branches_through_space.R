library(tidyverse)
library(binom)
library(paletteer)
library(ggpubr)
library(ggrepel)
# Load adaptive branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
                                        header = F) %>% unlist() %>% as.character()


# Load country mapping info
df_country_mapping <- read_csv('../../params/mapping/country_to_continent_mapping.csv') %>% 
  select(country, continent, subregion) %>% distinct()

# Load inferred county / subregion / continent
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


df_CI_prop_salt_branches <- binom.confint(
  x = df_prop_salt_per_location$n_saltational_branches,
  n = df_prop_salt_per_location$n_tot_branches,
  methods = "wilson"
) %>% 
  select(lower, upper) %>% 
  rename(lower_prop_salt = lower, upper_prop_salt = upper)

## Define colors
palette_africa <- paletteer_dynamic(palette = "cartography::blue.pal", n = 20)
palette_asia <- paletteer_dynamic(palette = "cartography::red.pal", n = 20)
palette_north_america <- paletteer_dynamic(palette = "cartography::sand.pal", n = 20)
palette_south_america <- paletteer_dynamic(palette = "cartography::brown.pal", n = 20)
palette_europe <- paletteer_dynamic(palette = "cartography::green.pal", n = 20)
palette_oceania <- paletteer_dynamic(palette = "cartography::purple.pal", n = 20)

df_country_mapping_with_colors <- df_country_mapping %>% 
  mutate(col_continent = case_when(continent == 'Asia' ~ palette_asia[10],
                                   continent == 'Africa' ~ palette_africa[10],
                                   continent == 'South America' ~ palette_south_america[10],
                                   continent == 'North America' ~ palette_north_america[10],
                                   continent == 'Europe' ~ palette_europe[10],
                                   continent == 'Oceania' ~ palette_oceania[10],
                                   T ~ NA)) %>% 
  select(country, continent, subregion, col_continent)


# Plot proportion of branches with saltational branches in locations
threshold_branches <- 5000
plt_continents <- df_prop_salt_per_location %>% bind_cols(df_CI_prop_salt_branches) %>% 
  filter(inferred_location != 'Multiple descendant locations', n_tot_branches >= threshold_branches,
         scale == 'Continent') %>%
  ggplot(aes(x = reorder(inferred_location, prop_saltational_branches), y = prop_saltational_branches, colour = inferred_location)) +
  geom_point() +
  geom_linerange(aes(ymin = lower_prop_salt, ymax = upper_prop_salt)) +
  scale_y_continuous(name = 'Proportion of\nsaltational branches',
                     labels = scales::percent) +
  scale_x_discrete(name = '') +
  scale_colour_manual(breaks = c('North America', 'South America', 'Europe', 'Africa', 'Asia', 'Oceania'),
                      values = c(palette_north_america[10], palette_south_america[10], 
                                 palette_europe[10], palette_africa[10], 
                                 palette_asia[10], palette_oceania[10]),
                      guide = 'none') +
  coord_flip() +
  theme_bw()

plt_regions <- df_prop_salt_per_location %>% bind_cols(df_CI_prop_salt_branches) %>% 
  filter(inferred_location != 'Multiple descendant locations', n_tot_branches >= threshold_branches,
         scale == 'Region') %>% 
  left_join(df_country_mapping_with_colors, by = c('inferred_location' = 'subregion')) %>% 
  ggplot(aes(x = reorder(inferred_location, prop_saltational_branches), y = prop_saltational_branches, colour = continent)) +
  geom_point() +
  geom_linerange(aes(ymin = lower_prop_salt, ymax = upper_prop_salt)) +
  scale_y_continuous(name = 'Proportion of\nsaltational branches',
                     labels = scales::percent) +
  scale_x_discrete(name = '') +
  scale_colour_manual(breaks = c('North America', 'South America', 'Europe', 'Africa', 'Asia', 'Oceania'),
                      values = c(palette_north_america[10], palette_south_america[10], 
                                 palette_europe[10], palette_africa[10], 
                                 palette_asia[10], palette_oceania[10]),
                      guide = 'none') +
  coord_flip() +
  theme_bw()

plt_countries <- df_prop_salt_per_location %>% bind_cols(df_CI_prop_salt_branches) %>% 
  filter(inferred_location != 'Multiple descendant locations', n_tot_branches >= threshold_branches,
         scale == 'Country') %>% 
  left_join(df_country_mapping_with_colors, by = c('inferred_location' = 'country')) %>% 
  mutate(inferred_location = ifelse(inferred_location == 'United Kingdom', 'UK', inferred_location)) %>% 
  ggplot(aes(x = reorder(inferred_location, prop_saltational_branches), y = prop_saltational_branches, colour = continent)) +
  geom_point() +
  geom_linerange(aes(ymin = lower_prop_salt, ymax = upper_prop_salt)) +
  scale_y_continuous(name = 'Proportion of\nsaltational branches',
                     labels = scales::percent) +
  scale_x_discrete(name = '') +
  scale_colour_manual(breaks = c('North America', 'South America', 'Europe', 'Africa', 'Asia', 'Oceania'),
                      values = c(palette_north_america[10], palette_south_america[10], 
                                 palette_europe[10], palette_africa[10], 
                                 palette_asia[10], palette_oceania[10]),
                      guide = 'none') +
  coord_flip() +
  theme_bw()


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
     


#### Ecological comparison
vec_countries_to_keep <- df_prop_salt_per_location %>% bind_cols(df_CI_prop_salt_branches) %>% 
  filter(inferred_location != 'Multiple descendant locations', n_tot_branches >= threshold_branches,
         scale == 'Country') %>% 
  left_join(df_country_mapping, by = c('inferred_location' = 'country')) %>% 
  select(inferred_location) %>% unlist() %>% as.character() %>% unique()

df_n_samples_in_tree_per_countries <- df_inferred_countries %>% 
  filter(!startsWith(strain, 'node_')) %>% 
  filter(inferred_location %in% vec_countries_to_keep) %>% 
  mutate(inferred_location = ifelse(inferred_location == 'United Kingdom', 'UK', inferred_location)) %>%
  group_by(inferred_location) %>% summarise(n_tips = n())

countries_demographics <- read_csv('../../data/country_demographics.csv')
countries_ihme_deaths <- read_csv('../../data/ihme_burden_country.csv') %>% 
  group_by(country) %>% summarise(n_deaths = sum(n_deaths)) %>% 
  rename(n_deaths_hiv_lymphoma_leukemia = n_deaths) %>% 
  mutate(country = ifelse(country == 'United Kingdom', 'UK', country),
         country = ifelse(country == 'United States of America', 'USA', country))


get_proba_cluster_detected <- function(p_seq, R0){
  proba_det <- 1. + lambertWp(-(1. - p_seq) * R0 * exp(-R0)) / R0
  return(proba_det)
}

df_for_ecological_plot <- countries_demographics %>% 
  left_join(countries_ihme_deaths) %>% 
  left_join(df_n_samples_in_tree_per_countries, by = c('country' = 'inferred_location')) %>% 
  left_join(df_prop_salt_per_location %>% bind_cols(df_CI_prop_salt_branches) %>% 
              filter(inferred_location != 'Multiple descendant locations', n_tot_branches >= threshold_branches,
                     scale == 'Country') %>% 
              left_join(df_country_mapping_with_colors, by = c('inferred_location' = 'country')) %>% 
              mutate(inferred_location = ifelse(inferred_location == 'United Kingdom', 'UK', inferred_location)) %>% 
              select(inferred_location, n_saltational_branches, n_tot_branches, prop_saltational_branches, lower_prop_salt, upper_prop_salt, continent),
            by = c('country' = 'inferred_location')) %>% 
  mutate(prop_seq_proxy = n_tips /population_size_2022_worldbank,
         proxy_proba_cluster_detected = get_proba_cluster_detected(prop_seq_proxy, 0.8),
         proba_cluster_detected_0.01 = get_proba_cluster_detected(0.01, 0.8),
         prop_salt_branches_corrected = prop_saltational_branches * proba_cluster_detected_0.01 / 0.01 / proxy_proba_cluster_detected * prop_seq_proxy,
         lower_prop_salt_branches_corrected = lower_prop_salt * proba_cluster_detected_0.01 / 0.01 / proxy_proba_cluster_detected * prop_seq_proxy,
         upper_prop_salt_branches_corrected = upper_prop_salt * proba_cluster_detected_0.01 / 0.01 / proxy_proba_cluster_detected * prop_seq_proxy)


label_125_percent <- function(x) {
  m <- signif(x / 10^floor(log10(x)), 12)   # mantissa
  labs <- scales::label_percent(accuracy = 1e-2)(x)
  ifelse(m %in% c(1, 2, 5), labs, "")
}
label_125 <- function(x) {
  m <- signif(x / 10^floor(log10(x)), 12)   # mantissa
  labs <- scales::label_number(accuracy = 1e-4)(x)
  ifelse(m %in% c(1), labs, "")
}


df_for_ecological_plot %>% 
  summarise(cor_spearman = cor(hiv_deaths_ihme_2022/population_size_2022_worldbank * 1000, 
                               (n_saltational_branches)/n_tot_branches, 
                               method = 'spearman'),
            cor_pearson = cor(hiv_deaths_ihme_2022/population_size_2022_worldbank * 1000, 
                               (n_saltational_branches)/n_tot_branches, 
                               method = 'pearson'),
            
            cor_log_spearman = cor(log(hiv_deaths_ihme_2022/population_size_2022_worldbank * 1000), 
                                   log((n_saltational_branches + 1)/n_tot_branches), 
                                   method = 'spearman'),
            cor_log_pearson = cor(log(hiv_deaths_ihme_2022/population_size_2022_worldbank * 1000), 
                                  log((n_saltational_branches + 1)/n_tot_branches), 
                                  method = 'pearson'),
            )


df_for_ecological_plot %>% 
  summarise(cor_spearman = cor(n_deaths_hiv_lymphoma_leukemia/population_size_2022_worldbank * 1000, 
                               (n_saltational_branches)/n_tot_branches, 
                               method = 'spearman'),
            cor_pearson = cor(n_deaths_hiv_lymphoma_leukemia/population_size_2022_worldbank * 1000, 
                              (n_saltational_branches)/n_tot_branches, 
                              method = 'pearson'),
            
            cor_log_spearman = cor(log(n_deaths_hiv_lymphoma_leukemia/population_size_2022_worldbank * 1000), 
                                   log((n_saltational_branches + 1)/n_tot_branches), 
                                   method = 'spearman'),
            cor_log_pearson = cor(log(n_deaths_hiv_lymphoma_leukemia/population_size_2022_worldbank * 1000), 
                                  log((n_saltational_branches + 1)/n_tot_branches), 
                                  method = 'pearson'),
  )



df_for_ecological_plot %>% 
  filter(country != 'South Africa') %>% 
  summarise(cor_spearman = cor(hiv_deaths_ihme_2022/population_size_2022_worldbank * 1000, 
                               (n_saltational_branches)/n_tot_branches, 
                               method = 'spearman'),
            cor_pearson = cor(hiv_deaths_ihme_2022/population_size_2022_worldbank * 1000, 
                              (n_saltational_branches)/n_tot_branches, 
                              method = 'pearson'),
            
            cor_log_spearman = cor(log(hiv_deaths_ihme_2022/population_size_2022_worldbank * 1000), 
                                   log((n_saltational_branches + 1)/n_tot_branches), 
                                   method = 'spearman'),
            cor_log_pearson = cor(log(hiv_deaths_ihme_2022/population_size_2022_worldbank * 1000), 
                                  log((n_saltational_branches + 1)/n_tot_branches), 
                                  method = 'pearson')
  )

correlation_prop_saltational_branches_hiv_deaths <- df_for_ecological_plot %>% 
  ggplot(aes(x = hiv_deaths_ihme_2022/population_size_2022_worldbank * 1000, 
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

correlation_prop_saltational_branches_hiv_deaths_no_south_africa <- df_for_ecological_plot %>% 
  filter(country != 'South Africa') %>% 
  ggplot(aes(x = hiv_deaths_ihme_2022/population_size_2022_worldbank * 1000, 
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


pdf('../figure_spatio_temporal_pattern/cor_prop_salt_branches_hiv_deaths.pdf',
    height = 4.5, width = 5.)
plot(correlation_prop_saltational_branches_hiv_deaths)
dev.off()
png('../figure_spatio_temporal_pattern/cor_prop_salt_branches_hiv_deaths.png',
    height = 4.5, width = 5., res = 350, units = 'in')
plot(correlation_prop_saltational_branches_hiv_deaths)
dev.off()

pdf('../figure_spatio_temporal_pattern/cor_prop_salt_branches_hiv_deaths_no_SA.pdf',
    height = 4.5, width = 5.)
plot(correlation_prop_saltational_branches_hiv_deaths_no_south_africa)
dev.off()
png('../figure_spatio_temporal_pattern/cor_prop_salt_branches_hiv_deaths_no_SA.png',
    height = 4.5, width = 5., res = 350, units = 'in')
plot(correlation_prop_saltational_branches_hiv_deaths_no_south_africa)
dev.off()


correlation_prop_saltational_branches_hiv_prev <- df_for_ecological_plot %>% 
  ggplot(aes(x = hiv_prevalence_percent_ihme_2022 / 100, 
             y = (n_saltational_branches + 1)/n_tot_branches)) +
  geom_point(aes(colour = continent)) +
  geom_text_repel(aes(label = country, colour = continent), max.overlaps = Inf) +
  scale_colour_manual(breaks = c('North America', 'South America', 'Europe', 'Africa', 'Asia', 'Oceania'),
                      values = c(palette_north_america[10], palette_south_america[10], 
                                 palette_europe[10], palette_africa[10], 
                                 palette_asia[10], palette_oceania[10]),
                      guide = 'none') +
  scale_x_continuous(trans = 'log', name = 'HIV prevalence',
                     breaks = c(seq(1e-4, 9e-4, 1e-4), seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 9e-1, 1e-1),
                                seq(1e0, 9e0, 1e0)),
                     labels = label_125_percent) + 
  scale_y_continuous(trans = 'log', name = 'Proportion of branches that are saltational',
                     breaks = c(seq(1e-4, 9e-4, 1e-4), seq(1e-3, 9e-3, 1e-3)),
                     labels = label_125_percent) +
  stat_cor(method = 'pearson', label.x = -4., label.y = -9) +
  stat_cor(method = 'spearman', label.x = -4., label.y = -9.3, cor.coef.name = 'rho') +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 11)) 

pdf('../figure_spatio_temporal_pattern/cor_prop_salt_branches_hiv_prev.pdf',
    height = 4.5, width = 6.)
plot(correlation_prop_saltational_branches_hiv_prev)
dev.off()
png('../figure_spatio_temporal_pattern/cor_prop_salt_branches_hiv_prev.png',
    height = 4.5, width = 6., res = 350, units = 'in')
plot(correlation_prop_saltational_branches_hiv_prev)
dev.off()
