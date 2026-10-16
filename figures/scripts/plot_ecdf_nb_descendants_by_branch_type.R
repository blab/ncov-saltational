### This script reproduces figures illustrating the 
### frequency at which saltational branches have descendants
### by branch type

library(tidyverse)
library(scales)
library(ggpubr)
source('config_plot.R')

#############################
### LOADING MODEL RESULTS ###
#############################

# Load name of saltational branches
name_saltational_branches <- 
  read.table('../../results/saltational_branches/name_saltational_branches_CV_1.8_BF_20.tsv', 
             header = F) %>% unlist() %>% as.character()

# Load number of descendants by node
df_descendants <- read_tsv('../../results/bte/df_n_descendants_per_node.tsv')

# Load mutations per branch
df_aa_muts <- read_tsv('../../results/mutations/nodes_aa_muts_summary_coded.tsv') %>% 
  filter(is_coded) %>% # Coded AA mut
  filter(aa_new != aa_initial) %>% # Non-synonymous
  mutate(full_mut = paste0(aa_muts_position, aa_new)) %>% 
  rename(strain = node_id) %>% 
  group_by(strain) %>% summarise(n_nonsyn = n())

# Functions to generate results
get_number_descendants_per_branch_type <- function(df_descendants, df_aa_muts, name_saltational_branches){
  df_desc_per_branch_type <- df_descendants %>%
    mutate(is_saltational = strain %in% name_saltational_branches) %>%
    left_join(df_aa_muts, by = "strain") %>%
    mutate(n_nonsyn = replace_na(n_nonsyn, 0),
           n_terminal_descendants = ifelse(n_terminal_descendants == 1, 0,
                                           n_terminal_descendants),
           n_terminal_descendants_crop = ifelse(n_terminal_descendants > 10, 10,
                                                n_terminal_descendants),
           n_terminal_descendants = ifelse(n_terminal_descendants == 0, 1,
                                           n_terminal_descendants),
           type = case_when(
             is_saltational ~ "Saltational",
             n_nonsyn >= 5 ~ ">= 5 NS muts",
             TRUE ~ "< 5 NS muts"
           )) %>% 
    select(strain, type, n_terminal_descendants) %>% 
    return()
}

df_nb_desc_per_branch <- get_number_descendants_per_branch_type(df_descendants, df_aa_muts, name_saltational_branches)

get_ecdf_bootstrap_band <- function(df_nb_desc_per_branch, n_bootstrap = 200){
  # One row per cluster, with its size and how many of its members fall in
  # the group
  df_for_ecdf <- df_nb_desc_per_branch
  
  vec_size_observed <- c(0, sort(unique(df_nb_desc_per_branch$n_terminal_descendants)))
  
  get_ecdf_over_individuals <- function(vec_row){
    ecdf(df_for_ecdf$n_terminal_descendants[vec_row])
  }
  
  # Run bootstrap
  mat_boot <- Reduce('cbind', lapply(1:n_bootstrap, FUN = function(i_bootstrap){
    vec_ids <- sample(nrow(df_for_ecdf), nrow(df_for_ecdf), replace = TRUE)
    get_ecdf_over_individuals(vec_ids)(vec_size_observed)
  }))
  
  tibble(n_terminal_descendants = vec_size_observed,
         ecdf_estimate = get_ecdf_over_individuals(seq_len(nrow(df_for_ecdf)))(vec_size_observed),
         ecdf_lower = apply(mat_boot, 1, quantile, probs = 0.025),
         ecdf_upper = apply(mat_boot, 1, quantile, probs = 0.975)) %>% 
    return()
}

df_ecdf_boostrap_band <- get_ecdf_bootstrap_band(df_nb_desc_per_branch %>% filter(type == '>= 5 NS muts'), 10000) %>% 
  mutate(log_reverse_ecdf = log(pmax(1 - ecdf_estimate, 1e-8)),
         log_reverse_ecdf_lower = log(pmax(1 - ecdf_upper, 1e-8)),
         log_reverse_ecdf_upper = log(pmax(1 - ecdf_lower, 1e-8)))
df_ecdf_boostrap_band_few_muts <- get_ecdf_bootstrap_band(df_nb_desc_per_branch %>% filter(type == '< 5 NS muts'), 1000) %>% 
  mutate(log_reverse_ecdf = log(pmax(1 - ecdf_estimate, 1e-8)),
         log_reverse_ecdf_lower = log(pmax(1 - ecdf_upper, 1e-8)),
         log_reverse_ecdf_upper = log(pmax(1 - ecdf_lower, 1e-8)))
df_ecdf_boostrap_band_saltational <- get_ecdf_bootstrap_band(df_nb_desc_per_branch %>% filter(type == 'Saltational'), 10000) %>% 
  mutate(log_reverse_ecdf = log(pmax(1 - ecdf_estimate, 1e-8)),
         log_reverse_ecdf_lower = log(pmax(1 - ecdf_upper, 1e-8)),
         log_reverse_ecdf_upper = log(pmax(1 - ecdf_lower, 1e-8)))

# Visualize ECDF of number of descendants by branch type
make_step_ribbon <- function(d) {
  step_boundary <- function(x, y) {
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    
    n <- length(x)
    
    tibble(
      x = c(x[1], rep(x[-1], each = 2)),
      y = c(y[1], as.vector(rbind(y[-n], y[-1])))
    )
  }
  
  d <- d %>%
    arrange(n_terminal_descendants)
  
  upper <- step_boundary(
    d$n_terminal_descendants,
    d$log_reverse_ecdf_upper
  ) %>%
    mutate(side = "upper")
  
  lower <- step_boundary(
    d$n_terminal_descendants,
    d$log_reverse_ecdf_lower
  ) %>%
    mutate(side = "lower")
  
  bind_rows(
    upper,
    lower %>% slice(n():1)
  )
}


plt_ecdf_nb_desc_by_branch_type <- ggplot(mapping = aes(x = n_terminal_descendants)) +
  geom_step(data = df_ecdf_boostrap_band_saltational, 
            aes(y = log_reverse_ecdf, type = 'Saltational', colour = 'Saltational')) +
  geom_polygon(data = make_step_ribbon(df_ecdf_boostrap_band_saltational), 
               aes(x = x, y = y, type = 'Saltational', fill = 'Saltational'),
               alpha = 0.3,
               colour = NA) +
  geom_step(data = df_ecdf_boostrap_band, 
            aes(y = log_reverse_ecdf, type = '>= 5 NS muts', colour = '>= 5 NS muts')) +
  geom_polygon(data = make_step_ribbon(df_ecdf_boostrap_band), 
               aes(x = x, y = y, type = '>= 5 NS muts', fill = '>= 5 NS muts'),
               alpha = 0.3,
               colour = NA) +
  geom_step(data = df_ecdf_boostrap_band_few_muts, 
            aes(y = log_reverse_ecdf, type = '< 5 NS muts', colour = '< 5 NS muts')) +
  geom_polygon(data = make_step_ribbon(df_ecdf_boostrap_band_few_muts), 
               aes(x = x, y = y, type = '< 5 NS muts', fill = '< 5 NS muts'),
               alpha = 0.3,
               colour = NA) +
  scale_colour_manual(
    breaks = c("Saltational", "< 5 NS muts", ">= 5 NS muts"),
    labels = c('Saltational', 
               'Non-saltational with < 5 nonsyn mutations',
               'Non-saltational with at least 5 nonsyn mutations'),
    values = c(Saltational = col_adaptive,
               `< 5 NS muts` = "gray50",
               `>= 5 NS muts` = "gray5"),
    name = "Branch type"
  ) +
  scale_fill_manual(
    breaks = c("Saltational", "< 5 NS muts", ">= 5 NS muts"),
    values = c(Saltational = col_adaptive,
               `< 5 NS muts` = "gray50",
               `>= 5 NS muts` = "gray5"),
    name = "Branch type", guide = "none"
  ) +
  scale_x_continuous(name = 'Number of terminal descendants (n)') +
  scale_y_continuous(name = 'Proportion of branches\nwith ≥ n terminal descendants', breaks = seq(-7, 0), 
                     labels = c(expression(10^{-7}), expression(10^{-6}), expression(10^{-5}),
                                expression(10^{-4}), expression(10^{-3}), expression(10^{-2}),
                                expression(10^{-1}), expression(10^{0}))) +
  coord_cartesian(xlim = c(-1, 500), ylim = c(-7 , 0)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'inside',
        legend.background = element_blank(),
        legend.position.inside = c(0.6, 0.8),
        legend.direction = 'vertical',
        legend.title.position = 'top',
        axis.text = element_text(size = 12))

## Proportion of branches with descendants
plt_prop_branches_with_desc_by_branch_type <- df_nb_desc_per_branch %>% 
  group_by(type) %>% 
  summarise(n_with_desc = sum(n_terminal_descendants > 1),
            n_tot = n()) %>% 
  rowwise() %>%
  mutate(prop = n_with_desc / n_tot,
         ci = list(
           binom.confint(
             x = n_with_desc,
             n = n_tot,
             method = "wilson"
           )[c("lower", "upper")]
         ),
         lower_prop = as.numeric(ci[1]), 
         upper_prop = as.numeric(ci[2])) %>% 
  select(-ci) %>% 
  ggplot(aes(x = type, y = prop, fill = type)) +
  geom_bar(stat = 'identity', alpha = 0.8) +
  geom_linerange(aes(ymin = lower_prop, ymax = upper_prop)) +
  coord_flip() +
  scale_fill_manual(
    breaks = c("Saltational", "< 5 NS muts", ">= 5 NS muts"),
    labels = c('Saltational', 
               'Non-saltational with < 5 nonsyn mutations',
               'Non-saltational with at least 5 nonsyn mutations'),
    values = c(Saltational = col_adaptive,
               `< 5 NS muts` = "gray50",
               `>= 5 NS muts` = "gray5"),
    name = "Branch type", guide = 'none'
  ) +
  scale_x_discrete(name = '', 
                   breaks = c("Saltational", "< 5 NS muts", ">= 5 NS muts"),
                   labels = c('Saltational', 
                              'Non-saltational\nwith < 5 nonsyn\nmutations',
                              'Non-saltational\nwith at least 5 nonsyn\nmutations')) +
  scale_y_continuous(name = 'Proportion of branches with descendants',
                     limits = c(0., NA),
                     expand = expansion(mult = c(0., NA))) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11))
  

panel_desc_count_distrib_by_branch_type <- 
  ggarrange(
    ggarrange(ggplot() + theme_void(),
              plt_prop_branches_with_desc_by_branch_type,
              ggplot() + theme_void(),
              ncol = 3, widths = c(0.1, 1., 0.1)),
    plt_ecdf_nb_desc_by_branch_type, 
    labels = 'AUTO', nrow = 2,
    heights = c(0.5, 1.))



png('../figure_descendant_analysis/panel_dist_nb_desc_by_branch_type.png',
    height = 6., width = 5.5, res = 350, units = 'in')
plot(panel_desc_count_distrib_by_branch_type)
dev.off()
