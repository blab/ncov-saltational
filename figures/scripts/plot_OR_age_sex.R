library(tidyverse)
library(binom)

# Load saltational branches names
name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', header = F) %>%
  unlist() %>% as.character()

# Load all terminal branches
df_child_parent <- read_tsv('../../results/bte/df_parent_child_relationship.tsv') %>% 
  filter(! startsWith(strain, 'node_'))

# Load metadata
metadata_tips <- read_tsv('../../results/saltational_branches/metadata_taxonium_CV_1.8_BF_20.tsv')

metadata_tips %>% filter(! is.na(sex)) %>% 
  mutate(is_saltational = strain %in% name_saltational_branches) %>% 
  group_by(sex, is_saltational) %>% 
  summarise(n_branches = n())

# Global - not significant (first glance at the continent - doesn't look significant either)
fisher.test(
  x = matrix(c(199, 227, 431489, 488799), ncol = 2, nrow = 2)
)

# 
df_prop_branches_per_age_group <- metadata_tips %>% 
  filter(! is.na(age)) %>% 
  mutate(age_group = cut(age, breaks = c(seq(0, 80, 10), Inf), include.lowest = T, right = F)) %>% 
  mutate(is_saltational = strain %in% name_saltational_branches) %>% 
  group_by(age_group) %>% 
  summarise(n_branches = n(), n_saltational = sum(is_saltational), n_non_saltational = sum(! is_saltational)) %>% 
  mutate(prop_saltational = n_saltational/n_branches)

df_CI_df_prop_branches_per_age_group <- binom.confint(
  x = df_prop_branches_per_age_group$n_saltational,
  n = df_prop_branches_per_age_group$n_branches,
  methods = "wilson"
) %>% 
  select(lower, upper) %>% 
  rename(lower_prop_salt = lower, upper_prop_salt = upper)


df_prop_branches_per_age_group %>% 
  bind_cols(df_CI_df_prop_branches_per_age_group) %>% 
  ggplot(aes(x = age_group, y = prop_saltational)) +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = lower_prop_salt, ymax = upper_prop_salt), position = position_dodge(0.4)) +
  coord_cartesian(ylim = c(0., NA))


df_prop_branches_per_age_group_per_continent <- metadata_tips %>% 
  filter(! is.na(age)) %>% 
  mutate(age_group = cut(age, breaks = c(0, 18, 60, Inf), include.lowest = T, right = F)) %>% 
  mutate(is_saltational = strain %in% name_saltational_branches) %>% 
  group_by(age_group, continent) %>% 
  summarise(n_branches = n(), n_saltational = sum(is_saltational), n_non_saltational = sum(! is_saltational)) %>% 
  mutate(prop_saltational = n_saltational/n_branches)

df_CI_df_prop_branches_per_age_group_per_continent <- binom.confint(
  x = df_prop_branches_per_age_group_per_continent$n_saltational,
  n = df_prop_branches_per_age_group_per_continent$n_branches,
  methods = "wilson"
) %>% 
  select(lower, upper) %>% 
  rename(lower_prop_salt = lower, upper_prop_salt = upper)



df_prop_branches_per_age_group_per_continent %>% 
  bind_cols(df_CI_df_prop_branches_per_age_group_per_continent) %>% 
  ggplot(aes(x = age_group, y = prop_saltational, colour = continent, group = continent)) +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = lower_prop_salt, ymax = upper_prop_salt), position = position_dodge(0.4)) +
  coord_cartesian(ylim = c(0., 0.005))


tmp_metadata_tips <- metadata_tips %>% 
  filter(! is.na(age), ! is.na(sex), continent == 'North America') %>% 
  mutate(age_group = cut(age, breaks = c(0, 80, Inf), include.lowest = T, right = F))

log_reg <- glm(data = tmp_metadata_tips, formula = is_saltational ~ age_group*sex)
summary(log_reg)
