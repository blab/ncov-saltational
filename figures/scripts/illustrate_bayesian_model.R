library(tidyverse)
library(extraDistr)
library(zipfR)
source('../../scripts/bayesian-model/utils_model.R')
source('config_plot.R')

### This script reproduces figures used to illustrate how the 
### Bayesian model infers the posterior distribution of dN, dS and omega

# Load prior parametrisation
df_prior <- read_tsv('../../results/input_model/df_prior_parametrization.tsv')

# Function to run model on Spike
get_posterior_from_example <- function(df_prior, range_S_values, range_N_values, CV_prior){
  
  # Define prior parameters
  prior_param <- list('alpha_N' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dN[df_prior$gene == 'S']), cv = CV_prior)['alpha'], 
                      'beta_N' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dN[df_prior$gene == 'S']), cv = CV_prior)['beta'],
                      'alpha_S' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dS[df_prior$gene == 'S']), cv = CV_prior)['alpha'], 
                      'beta_S' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dS[df_prior$gene == 'S']), cv = CV_prior)['beta'])
  
  # Define site information
  mean_non_syn_sites <- round(df_prior$mean_nonsyn_sites[df_prior$gene == 'S'])
  mean_syn_sites <- round(df_prior$mean_syn_sites[df_prior$gene == 'S'])
  
  df_mut_combinations <- as_tibble(expand.grid(N = range_N_values, S = range_S_values))
  
  vec_dN <- seq(0, 3e-4, 1e-7)
  vec_dS <- seq(0, 3e-4, 1e-7)
  vec_omega <- seq(0.0, 6., 0.01)
  
  df_posterior <- Reduce('bind_rows', lapply(1:nrow(df_mut_combinations), FUN = function(i_row){
    data_gene <- list('N' = df_mut_combinations$N[i_row], 
                      'S' = df_mut_combinations$S[i_row], 
                      'LN' = mean_non_syn_sites,
                      'LS' = mean_syn_sites)
    
    mod <- create_model(data_gene, prior_param)
    mod$get_quantile_posterior_dS(c(0.5))
    mod$get_quantile_posterior_dS(c(0.975))
    mod$get_quantile_posterior_dS(c(0.025))
    
    mod$get_quantile_posterior_dN(c(0.5))
    mod$get_quantile_posterior_dN(c(0.975))
    mod$get_quantile_posterior_dN(c(0.025))
    
    df_dN <- tibble(value = vec_dN,
                    posterior = mod$get_posterior_density_dN(vec_dN),
                    prior = mod$get_prior_density_dN(vec_dN)) %>% 
      mutate(param = 'dN',
             i_row = i_row)
    
    df_dS <- tibble(value = vec_dS,
                    posterior = mod$get_posterior_density_dS(vec_dS),
                    prior = mod$get_prior_density_dS(vec_dS)) %>% 
      mutate(param = 'dS',
             i_row = i_row)
    
    df_omega <- tibble(value = vec_omega,
                       posterior = exp(mod$get_log_posterior_density_omega(vec_omega)),
                       prior = exp(mod$get_log_prior_density_omega(vec_omega))) %>% 
      mutate(param = 'omega',
             i_row = i_row)
    
    bind_rows(df_dN, df_dS, df_omega) %>% 
      mutate('N' = df_mut_combinations$N[i_row], 
             'S' = df_mut_combinations$S[i_row])
  }))
  
  return(df_posterior)
}

get_omega_estim_from_example <- function(df_prior, range_S_values, range_N_values, CV_prior){
  
  # Define prior parameters
  prior_param <- list('alpha_N' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dN[df_prior$gene == 'S']), cv = CV_prior)['alpha'], 
                      'beta_N' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dN[df_prior$gene == 'S']), cv = CV_prior)['beta'],
                      'alpha_S' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dS[df_prior$gene == 'S']), cv = CV_prior)['alpha'], 
                      'beta_S' = get_gamma_param_from_mean_CV(mu = exp(df_prior$log_mean_dS[df_prior$gene == 'S']), cv = CV_prior)['beta'])
  
  # Define site information
  mean_non_syn_sites <- round(df_prior$mean_nonsyn_sites[df_prior$gene == 'S'])
  mean_syn_sites <- round(df_prior$mean_syn_sites[df_prior$gene == 'S'])
  
  df_mut_combinations <- as_tibble(expand.grid(N = range_N_values, S = range_S_values))
  
  
  df_omega_info <- Reduce('bind_rows', lapply(1:nrow(df_mut_combinations), FUN = function(i_row){
    data_gene <- list('N' = df_mut_combinations$N[i_row], 
                      'S' = df_mut_combinations$S[i_row], 
                      'LN' = mean_non_syn_sites,
                      'LS' = mean_syn_sites)
    
    mod <- create_model(data_gene, prior_param)
    
    mod$posterior_proba_omega_greater_1
    mod$prior_proba_omega_greater_1
    
    tibble(posterior_proba_omega_greater_1 = mod$posterior_proba_omega_greater_1,
           prior_proba_omega_greater_1 = mod$prior_proba_omega_greater_1,
           'N' = df_mut_combinations$N[i_row], 
           'S' = df_mut_combinations$S[i_row])
  }))
  
  return(df_omega_info)
}

# Posterior distribution of dN and dS for different combinations of N and S on Spike
labs <- c(
  `2`  = "paste(S[i], ' = ', 1, ', ', N[i], ' = ', 2)",
  `10` = "paste(S[i], ' = ', 1, ', ', N[i], ' = ', 10)"
)
plt_posterior_rates <- get_posterior_from_example(df_prior, range_S_values = c(1), 
                                                  range_N_values = c(2, 10), CV_prior = 1.8) %>% 
  filter(param %in% c('dN', 'dS')) %>% 
  mutate(N_label = paste0("paste(N[i], ' = ', ", N, ")")) %>% 
  ggplot(aes(x = value, y = posterior,
             group = interaction(N, S, param), 
             color = as.factor(param),
             linetype = as.factor(N))) +
  geom_line()  +
  facet_wrap(. ~ N, labeller = as_labeller(labs, label_parsed)) +
  scale_x_continuous(name = 'Rate value', breaks = seq(0, 3e-4, 1e-4),
                     labels = c('0', 
                                expression(paste(1, '\u00b7', 10^{-4})), 
                                expression(paste(2, '\u00b7', 10^{-4})),
                                expression(paste(3, '\u00b7', 10^{-4})))) +
  scale_y_continuous(name = 'Posterior density', breaks = seq(0, 2e4, 1e4),
                     labels = c(0, 
                                expression(paste(1, '\u00b7', 10^{4})), 
                                expression(paste(2, '\u00b7', 10^{4})))) +
  scale_colour_manual(breaks = c('dS', 'dN'), values = c(col_synonymous, col_non_synonymous), name = 'Rate') +
  scale_linetype_manual(name = expression(N[i]), values = c(1, 2), guide = 'none') +
  theme_bw() + 
  theme(strip.background = element_rect(fill = 'gray22'),
        strip.text = element_text(colour = 'white', size = 11),
        axis.text = element_text(size = 11))

pdf('../figure_model_illustration/plt_posterior_dNdS_rates.pdf', height = 2.5, width = 6.)
plot(plt_posterior_rates)
dev.off()
png('../figure_model_illustration/plt_posterior_dNdS_rates.png', height = 2.5, width = 6.,
    res = 350, units = 'in')
plot(plt_posterior_rates)
dev.off()

## Plot posterior distribution of omega
plt_posterior_omega <- get_posterior_from_example(df_prior, range_S_values = c(1), 
                                                  range_N_values = c(2, 10), CV_prior = 1.7) %>% 
  filter(param %in% c('omega')) %>% 
  ggplot(aes(x = value, y = posterior, linetype = as.factor(N))) +
  geom_line(aes(y = prior, color = 'Prior')) +
  geom_line(aes(color = 'Posterior')) +
  scale_linetype_manual(name = expression(N[i]), values = c(1, 2)) +
  scale_color_manual(breaks = c('Prior', 'Posterior'), values = c(col_prior, col_posterior), name = '') +
  scale_x_continuous(name = expression(omega)) +
  scale_y_continuous(name = 'Posterior density',
                     expand = expansion(mult = c(0., 0.02))) +
  theme_classic() +
  coord_cartesian(ylim = c(0, .9), xlim = c(0, 5.)) +
  theme(panel.grid = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.7, 0.7),
        legend.background = element_blank(),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11))

pdf('../figure_model_illustration/plt_posterior_omega.pdf', height = 2.5, width = 3.)
plot(plt_posterior_omega)
dev.off()
png('../figure_model_illustration/plt_posterior_omega.png', height = 2.5, width = 3.,
    res = 350, units = 'in')
plot(plt_posterior_omega)
dev.off()

## Plot probability of omega being greater than 1
labs_split <- c(
  `Prior` = "Prior",
  `2`  = "paste(N[i], ' = ', 2)",
  `10` = "paste(N[i], ' = ', 10)"
)

plt_proba_omega_above_1 <- get_omega_estim_from_example(df_prior, range_S_values = c(1), 
                                                        range_N_values = c(2, 10), 
                                                        CV_prior = 1.8) %>% 
  pivot_longer(cols = 1:2, names_pattern = "(.*)_proba_omega_greater_1", values_to = 'proba_omega_greater_1') %>% 
  filter(name == 'posterior' | N == 2) %>% 
  mutate(x_axis = ifelse(name == 'prior', 'Prior', as.character(N)),
         x_axis = factor(x_axis, levels = c('Prior', '2', '10'))) %>% 
  ggplot(aes(x = x_axis, y = proba_omega_greater_1, fill = name, color = name, linetype = as.factor(N))) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(labels = function(x) parse(text = labs_split[x])) +
  scale_y_continuous(name = expression(paste('Probability of ', omega > 1)),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0., 0.02)),
                     breaks = seq(0., 1., 0.2)) +
  scale_colour_manual(breaks = c('prior', 'posterior'), values = c(col_prior, col_posterior), guide = 'none') +
  scale_fill_manual(breaks = c('prior', 'posterior'), values = c(col_prior, 'gray92'), guide = 'none') +
  scale_linetype_manual(name = expression(N[i]), values = c(1, 2), guide = 'none') +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11))

pdf('../figure_model_illustration/plt_proba_omega_above_1.pdf', height = 2.5, width = 2.5)
plot(plt_proba_omega_above_1)
dev.off()
png('../figure_model_illustration/plt_proba_omega_above_1.pdf', 
    height = 2.5, width = 2.5, res = 350, units = 'in')
plot(plt_proba_omega_above_1)
dev.off()

## Illustrate Bayes Factor values across (N,S values) for Spike
df_BF_Spike <- get_omega_estim_from_example(df_prior, range_S_values = 0:8, 
                                            range_N_values = 0:20,
                                            CV_prior = 1.8) %>% 
  mutate(posterior_odds = posterior_proba_omega_greater_1 / (1. - posterior_proba_omega_greater_1),
         prior_odds = prior_proba_omega_greater_1 / (1. - prior_proba_omega_greater_1),
         BF = posterior_odds / prior_odds,
         log10_BF = log10(BF)) %>%
  mutate(log10_BF_bins = case_when(log10_BF < 0 ~ '< 0',
                                   log10_BF < 0.5 ~ '0-0.5 (Barely worth mentioning)',
                                   log10_BF < 1 ~ '0.5-1.0 (Substantial)',
                                   log10_BF < 1.5 ~ '1.0-1.5 (Strong)',
                                   log10_BF <  2. ~ '1.5-2.0 (Very strong)', 
                                   T ~ '2+ (Decisive)'),
         log10_BF_bins = factor(log10_BF_bins, levels = c('< 0',
                                                          '0-0.5 (Barely worth mentioning)', 
                                                          '0.5-1.0 (Substantial)',
                                                          '1.0-1.5 (Strong)',
                                                          '1.5-2.0 (Very strong)',
                                                          '2+ (Decisive)')))


df_BF_Spike_Kass_Raftery <- get_omega_estim_from_example(df_prior, range_S_values = 0:8, 
                                                        range_N_values = 0:20,
                                                        CV_prior = 1.8) %>% 
  mutate(posterior_odds = posterior_proba_omega_greater_1 / (1. - posterior_proba_omega_greater_1),
         prior_odds = prior_proba_omega_greater_1 / (1. - prior_proba_omega_greater_1),
         BF = posterior_odds / prior_odds) %>%
  mutate(BF_bins_KR = case_when(BF < 1. ~ 'Negative (< 1)',
                                BF < 3 ~ 'Weak (1-3)',
                                BF < 20 ~ 'Positive (3-20)',
                                BF < 150 ~ 'Strong (20-150)',
                                T ~ 'Very strong (> 150)'),
         BF_bins_KR = factor(BF_bins_KR, levels = c('Negative (< 1)',
                                                    'Weak (1-3)',
                                                    'Positive (3-20)',
                                                    'Strong (20-150)',
                                                    'Very strong (> 150)')))


heatmap_BF_func_SN <- df_BF_Spike %>% 
  ggplot(aes(x = N, y = S)) +
  geom_tile(aes(fill = log10_BF_bins)) +
  scale_fill_manual(values = brewer.pal(6, 'BuPu'),
                    name = expression(log[10](BF[omega > 1]))) +
  scale_y_continuous(name = expression(S[i]), breaks = seq(0., 10., 1.),
                     expand = expansion(mult = c(0., 0.))) +
  scale_x_continuous(name = expression(N[i]), 
                     breaks = seq(0., 20., 1.),
                     labels = c(0, '', 2, '', 4, '', 6, '', 8, '', 10, '', 
                                12, '', 14, '', 16, '', 18, '', 20),
                     expand = expansion(mult = c(0., 0.))) +
  coord_fixed(xlim = c(-0.5, 20.5)) +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        strip.background = element_rect(fill = 'gray22'),
        strip.text = element_text(colour = 'white', size = 11),
        panel.grid.minor = element_blank())

heatmap_BF_func_SN_Kass_Raftery <- df_BF_Spike_Kass_Raftery %>% 
  ggplot() +
  geom_tile(aes(x = N, y = S, fill = BF_bins_KR)) +
  geom_rect(xmin = 1.5, xmax = 2.5, ymin = 0.5, ymax = 1.5, fill = NA, color = 'black') + 
  geom_rect(xmin = 9.5, xmax = 10.5, ymin = 0.5, ymax = 1.5, fill = NA, color = 'black', linetype = 'dashed') + 
  scale_fill_manual(values = brewer.pal(5, 'BuPu'),
                    name = 'Evidence for\npositive selection\n(BF range)') +
  scale_y_continuous(name = expression(paste('Synonymous mutations count ', S[i])), breaks = seq(0., 10., 1.),
                     expand = expansion(mult = c(0., 0.))) +
  scale_x_continuous(name = expression(paste('Nonsynonymous mutations count ', N[i])), 
                     breaks = seq(0., 20., 1.),
                     labels = c(0, '', 2, '', 4, '', 6, '', 8, '', 10, '', 
                                12, '', 14, '', 16, '', 18, '', 20),
                     expand = expansion(mult = c(0., 0.))) +
  coord_fixed(xlim = c(-0.5, 20.5)) +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        strip.background = element_rect(fill = 'gray22'),
        strip.text = element_text(colour = 'white', size = 11),
        panel.grid.minor = element_blank())

plot(heatmap_BF_func_SN_Kass_Raftery)

pdf('../figure_model_illustration/heatmap_BF_spike_func_N_S.pdf',
    height = 2.5, width = 7.)
plot(heatmap_BF_func_SN)
dev.off()
png('../figure_model_illustration/heatmap_BF_spike_func_N_S.png',
    height = 2.5, width = 7., res = 350, units = 'in')
plot(heatmap_BF_func_SN)
dev.off()

pdf('../figure_model_illustration//heatmap_BF_spike_func_N_S_KS.pdf',
    height = 2.5, width = 8.)
plot(heatmap_BF_func_SN_Kass_Raftery)
dev.off()
png('../figure_model_illustration//heatmap_BF_spike_func_N_S_KS.png',
    height = 2.5, width = 8., res = 350, units = 'in')
plot(heatmap_BF_func_SN_Kass_Raftery)
dev.off()

