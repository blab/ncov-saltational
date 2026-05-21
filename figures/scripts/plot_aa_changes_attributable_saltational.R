library(tidyverse)

# Load name of saltational branches
vec_name_saltational_branches <- read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv') %>%
  unlist() %>% as.character()

# Load AA mutation
df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  rename(strain = node_id) %>% 
  filter(is_coded, aa_initial != aa_new) %>% 
  mutate(is_saltational = strain %in% vec_name_saltational_branches) %>% 
  summarise(n_aa_changes = n(),
            n_aa_changes_salt = sum(is_saltational))

# Load cum AA changes from root and only retain tips
df_cum_aa_changes <- read_tsv('../../results/saltational_branches/df_cum_aa_changes_from_root_CV_1.8_BF_20.tsv') %>% 
  rename(strain = node_id) %>% 
  filter(! startsWith(strain, 'node_'))

# Load metadata (date) of tips
metadata_clean_dates <- read_tsv('../../results/sra_metadata/sra_metadata_just_clean_dates.tsv') %>% 
  filter(strain %in% df_cum_aa_changes$strain)

df_for_plot <- df_cum_aa_changes %>% 
  left_join(metadata_clean_dates, by = 'strain') %>% 
  filter(! is.na(date)) %>% 
  mutate(time = as.numeric(date - as.Date('2020-01-01'), unit = 'days'))

# Cumulative AA changes through time
plt_cum_aa_changes_through_time <- df_for_plot %>% 
  slice_sample(prop = 0.01) %>% 
  ggplot(aes(x = date)) +
  geom_point(aes(y = cum_aa_changes_from_root, color = 'All'), alpha = 0.1) +
  geom_point(aes(y = cum_saltational_aa_changes_from_root, color = 'From saltational'), alpha = 0.1) +
  scale_x_date(name = 'Sequence collection date') +
  scale_y_continuous(name = 'Cumulative AA changes from root') +
  scale_color_manual(values = c('darkgrey', 'orange2'), name = '') +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1.)))

png('../figure_prop_aa_change_from_salt/plt_cum_aa_changes_through_time.png',
    height = 3.5, width = 5., res = 350, units = 'in')
plot(plt_cum_aa_changes_through_time)
dev.off()


plt_prop_aa_changes_attributable_to_salt_2024 <- df_for_plot %>% 
  filter(date >= as.Date('2024-01-01')) %>% 
  mutate(y_to_plot = cum_saltational_aa_changes_from_root/cum_aa_changes_from_root) %>% 
  summarise(y05 = quantile(y_to_plot, 0.025),
            y25 = quantile(y_to_plot, 0.25),
            y50 = quantile(y_to_plot, 0.5),
            y75 = quantile(y_to_plot, 0.75),
            y95 = quantile(y_to_plot, 0.975)) %>% 
  ggplot(aes(x = 1)) + 
  geom_boxplot(
    aes(ymin = y05, lower = y25, middle = y50, upper = y75, ymax = y95),
    stat = "identity"
  ) +
  scale_x_continuous(breaks = c(), name = '') +
  scale_y_continuous(name = 'Proportion of AA changes attributable\nto saltational evolution',
                     breaks = seq(0., 0.4, 0.01)) +
  theme_classic() +
  coord_flip()

png('../figure_prop_aa_change_from_salt/plt_prop_aa_changes_attributable_to_salt_2024.png',
    height = 1.5, width = 6., res = 350, units = 'in')
plot(plt_prop_aa_changes_attributable_to_salt_2024)
dev.off()









loess_cum_aa_changes <-  loess(cum_aa_changes_from_root ~ time, data = df_for_plot)
loess_cum_aa_changes_salt <-  loess(cum_saltational_aa_changes_from_root ~ time, data = df_for_plot)

df_for_plot %>% 
  mutate(pred_loess_cum_aa_changes = predict(loess_cum_aa_changes),
         pred_loess_cum_aa_changes_salt = predict(loess_cum_aa_changes_salt)) %>% 
  ggplot(aes(x = time + as.Date('2020-01-01'))) +
  geom_line(aes(y = pred_loess_cum_aa_changes, color = 'All')) +
  geom_line(aes(y = pred_loess_cum_aa_changes_salt, color = 'Salt')) 

df_for_plot %>% 
  filter(date >= as.Date('2024-01-01')) %>% 
  ggplot(aes(x = date, y = cum_saltational_aa_changes_from_root)) +
  geom_point(alpha = 0.2)


  

df_for_plot %>% 
  filter(date >= as.Date('2024-01-01')) %>% 
  ggplot(aes(y = cum_saltational_aa_changes_from_root)) +
  geom_boxplot()

df_for_plot %>% 
  filter(date >= as.Date('2024-01-01')) %>% 
  ggplot(aes(y = cum_aa_changes_from_root)) +
  geom_boxplot()

    filter(date >= as.Date('2024-01-01'))
get_minus_log_lik_2024 <- function(p){
  tmp_df <- df_for_plot %>% 
    filter(date >= as.Date('2024-01-01'))
  
  log_lik_contrib <- 
    log(p) * sum(tmp_df$cum_saltational_aa_changes_from_root) + 
    log(1. - p) * sum(tmp_df$cum_aa_changes_from_root - tmp_df$cum_saltational_aa_changes_from_root)
  return(- log_lik_contrib)
}
  
sapply(seq(0.0, 1., 0.1), FUN = get_log_lik_2024)

chisq_stat_95 <- qchisq(0.95, df = 1)

MLE_estim <- optim(par = c(0.5), fn = get_minus_log_lik_2024)
max_minus_log_lik <- MLE_estim$value
seq_p <- seq(0., 1., 1e-4)
vec_minus_log_lik <- sapply(seq_p, FUN = get_minus_log_lik_2024)
plot(seq_p, 2* (vec_minus_log_lik - max_minus_log_lik))

seq_p[2 * (vec_minus_log_lik - max_minus_log_lik) - chisq_stat_95 < 0]


