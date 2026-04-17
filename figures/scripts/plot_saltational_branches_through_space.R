library(tidyverse)
library(binom)
library(paletteer)
library(ggpubr)

# Load adaptive branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
                                        header = F) %>% unlist() %>% as.character()


# Load country mapping info
df_country_mapping <- read_csv('../../params/mapping/country_to_continent_mapping.csv')

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
     
