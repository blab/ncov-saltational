library(tidyverse)
library(pracma)

get_proba_cluster_detected <- function(p_seq, R0){
  proba_det <- 1. + lambertWp(-(1. - p_seq) * R0 * exp(-R0)) / R0
  return(proba_det)
}


expand.grid(log_proba_cluster_starts = seq(log(1e-5), log(1e-1), 0.01),
            log_p_seq = seq(log(1e-5), log(1e0), 0.01),
            R0 = c(0.8)) %>% 
  mutate(p_seq = exp(log_p_seq),
         proba_cluster_starts = exp(log_proba_cluster_starts)) %>% 
  mutate(proba_cluster_detected = get_proba_cluster_detected(p_seq, R0),
         expected_proportion_of_long_branches = proba_cluster_detected * proba_cluster_starts / p_seq) %>% 
  ggplot(aes(x = log_proba_cluster_starts, y = log_p_seq, fill = expected_proportion_of_long_branches)) +
  geom_tile() +
  scale_x_continuous(breaks = c(log(1e-5), log(1e-4), log(1e-3), log(1e-2), log(1e-1), log(1)),
                     labels = c(expression(10^{-5}), 
                                expression(10^{-4}), expression(10^{-3}), expression(10^{-2}),
                                expression(10^{-1}), expression(10^{0})),
                     name = 'Probability for an infection of resulting in a saltational branch') +
  scale_y_continuous(breaks = c(log(1e-5), log(1e-4), log(1e-3), log(1e-2), log(1e-1), log(1)),
                     labels = c(expression(10^{-5}), 
                                expression(10^{-4}), expression(10^{-3}), expression(10^{-2}),
                                expression(10^{-1}), expression(10^{0})),
                     name = 'Probability for an infection of being sequenced') +
  scale_fill_viridis_c(trans = 'log',
                       breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0),
                       labels = c(expression(10^{-4}), expression(10^{-3}), expression(10^{-2}),
                                  expression(10^{-1}), expression(10^{0})),
                       name = 'Expected proportion of\nbranches that are saltational')



get_proba_detect_cluster <- function(p_seq, cluster_size){
  return(
    1. - (1. - p_seq)^cluster_size
  )
}

expand.grid(log_cluster_size = seq(log(1), log(1e3), 0.01),
            log_p_seq = seq(log(1e-6), log(1e0), 0.01)) %>% 
  mutate(proba_detect_cluster = get_proba_detect_cluster(exp(log_p_seq), exp(log_cluster_size))) %>% 
  ggplot(aes(x = log_p_seq, y = log_cluster_size)) +
  geom_tile(aes(fill = proba_detect_cluster)) +
  scale_fill_viridis_c(trans = 'log',
                       breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0),
                       labels = c(expression(10^{-6}), expression(10^{-5}), 
                                  expression(10^{-4}), expression(10^{-3}), expression(10^{-2}),
                                  expression(10^{-1}), expression(10^{0})),
                       option = 'magma',
                       name = 'Probability for a descendant\nof the saltational branch\nof being sequenced') +
  scale_x_continuous(breaks = c(log(1e-6), log(1e-5), log(1e-4), log(1e-3), log(1e-2), log(1e-1), log(1)),
                     labels = c(expression(10^{-6}), expression(10^{-5}), 
                                expression(10^{-4}), expression(10^{-3}), expression(10^{-2}),
                                expression(10^{-1}), expression(10^{0})),
                     name = 'Probability for an infection of being sequenced',
                     expand = expansion(mult = c(0., 0.)),
                     limits = c(NA, log(1))) +
  scale_y_continuous(breaks = c(log(1), log(1e1), log(1e2), log(1e3)),
                     labels = c(expression(10^{0}), expression(10^{1}), 
                                expression(10^{2}), expression(10^{3})),
                     limits = c(NA, log(1e3)),
                     expand = expansion(mult = c(0., 0.)),
                     name = 'Number of descendants following a saltational event') 
  
