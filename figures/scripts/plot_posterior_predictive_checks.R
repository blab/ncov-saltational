library(tidyverse)
library(grid)
source('config_plot.R')


get_summary_from_posterior_simul <- function(gene, CV_prior){
  # Load posterior simulations
  df_posterior_simul <- read_tsv(paste0('../../results/posterior_predictive_checks/posterior_predictive_check_', gene, '_CV_', CV_prior, '.tsv'))
  
  # Get counts of S and N
  df_posterior_count_S <- df_posterior_simul %>% group_by(S_new, i_simul) %>% 
    summarise(n_branches_new = sum(n_branches_new)) %>% 
    ungroup() %>% 
    complete(S_new, i_simul, fill = list(n_branches_new = 0))
  
  df_posterior_count_N <- df_posterior_simul %>% group_by(N_new, i_simul) %>% 
    summarise(n_branches_new = sum(n_branches_new)) %>% 
    ungroup() %>% 
    complete(N_new, i_simul, fill = list(n_branches_new = 0))
  
  df_posterior_simul_count_NS <- df_posterior_simul %>% group_by(S_new, N_new, i_simul) %>% 
    summarise(n_branches_new = sum(n_branches_new)) %>% 
    ungroup() %>% 
    complete(S_new, N_new, i_simul, fill = list(n_branches_new = 0))
  
  # Summary of counts across posterior simulations
  PDF_S <- df_posterior_count_S %>% group_by(S_new) %>% 
    summarise(mean_n_branches_new = mean(n_branches_new),
              median_n_branches_new = median(n_branches_new),
              lower_n_branches_new = quantile(n_branches_new, 0.025),
              upper_n_branches_new = quantile(n_branches_new, 0.975),
              min_n_branches_new = min(n_branches_new),
              max_n_branches_new = max(n_branches_new))
  
  PDF_N <- df_posterior_count_N %>% group_by(N_new) %>% 
    summarise(mean_n_branches_new = mean(n_branches_new),
              median_n_branches_new = median(n_branches_new),
              lower_n_branches_new = quantile(n_branches_new, 0.025),
              upper_n_branches_new = quantile(n_branches_new, 0.975),
              min_n_branches_new = min(n_branches_new),
              max_n_branches_new = max(n_branches_new))
  
  PDF_NS <- df_posterior_simul_count_NS %>% 
    group_by(S_new, N_new) %>% 
    summarise(mean_n_branches_new = mean(n_branches_new),
              median_n_branches_new = median(n_branches_new),
              lower_n_branches_new = quantile(n_branches_new, 0.025),
              upper_n_branches_new = quantile(n_branches_new, 0.975),
              min_n_branches_new = min(n_branches_new),
              max_n_branches_new = max(n_branches_new)) 
  
  # Load actual data
  gene_data <- read_tsv(paste0('../../results/input_model/data_viridian_', gene, '_grouped.tsv'))
  df_true_count_S <- gene_data %>% group_by(synonymous_muts) %>% summarise(n_branches = sum(n_branches)) %>% rename(S = synonymous_muts)
  df_true_count_N <- gene_data %>% group_by(non_synonymous_muts) %>% summarise(n_branches = sum(n_branches)) %>% rename(N = non_synonymous_muts)
  
  df_true_count_SN <- gene_data %>% 
    group_by(synonymous_muts, non_synonymous_muts) %>% summarise(n_branches_true = sum(n_branches)) %>% 
    rename(S = synonymous_muts, N = non_synonymous_muts) %>% 
    group_by(S, N) %>% summarise(n_branches_true = sum(n_branches_true))
  
  
  # Merge PDF
  PDF_N <- PDF_N %>% rename(N = N_new) %>% left_join(df_true_count_N %>% rename(n_branches_true = n_branches)) %>% 
    mutate(n_branches_true = replace_na(n_branches_true, 0),
           mean_n_branches_new = replace_na(mean_n_branches_new, 0),
           median_n_branches_new = replace_na(median_n_branches_new, 0),
           lower_n_branches_new = replace_na(lower_n_branches_new, 0),
           upper_n_branches_new = replace_na(upper_n_branches_new, 0),
           min_n_branches_new = replace_na(min_n_branches_new, 0),
           max_n_branches_new = replace_na(max_n_branches_new, 0))
  PDF_S <- PDF_S %>% rename(S = S_new) %>% left_join(df_true_count_S %>% rename(n_branches_true = n_branches)) %>% 
    mutate(n_branches_true = replace_na(n_branches_true, 0),
           mean_n_branches_new = replace_na(mean_n_branches_new, 0),
           median_n_branches_new = replace_na(median_n_branches_new, 0),
           lower_n_branches_new = replace_na(lower_n_branches_new, 0),
           upper_n_branches_new = replace_na(upper_n_branches_new, 0),
           min_n_branches_new = replace_na(min_n_branches_new, 0),
           max_n_branches_new = replace_na(max_n_branches_new, 0))
  
  PDF_NS <- PDF_NS %>% rename(S = S_new, N = N_new) %>% 
    full_join(df_true_count_SN) %>% 
    mutate(median_n_branches_new = replace_na(median_n_branches_new, 0),
           mean_n_branches_new = replace_na(mean_n_branches_new, 0),
           lower_n_branches_new = replace_na(lower_n_branches_new, 0),
           upper_n_branches_new = replace_na(upper_n_branches_new, 0),
           min_n_branches_new = replace_na(min_n_branches_new, 0),
           max_n_branches_new = replace_na(max_n_branches_new, 0),
           n_branches_true = replace_na(n_branches_true, 0))
  
  return(
    list(
      'PDF_N' = PDF_N %>% mutate(gene = gene, CV = CV_prior),
      'PDF_S' = PDF_S %>% mutate(gene = gene, CV = CV_prior),
      'PDF_NS' = PDF_NS %>% mutate(gene = gene, CV = CV_prior)
    )
  )
}

vec_genes <- c('ORF1a', 'ORF1b', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF8', 'N', 'ORF9b', 'Genome')

list_posterior_simul_by_gene <- lapply(vec_genes, FUN = function(gene){
  get_summary_from_posterior_simul(gene, 1.5) 
})
list_posterior_simul_by_gene_CV_3 <- lapply(vec_genes, FUN = function(gene){
  get_summary_from_posterior_simul(gene, '3.0') 
})

df_PDF_N <- Reduce('bind_rows', lapply(list_posterior_simul_by_gene, `[[`, 'PDF_N')) %>% ungroup()
df_PDF_S <- Reduce('bind_rows', lapply(list_posterior_simul_by_gene, `[[`, 'PDF_S')) %>% ungroup()
df_PDF_NS <- Reduce('bind_rows', lapply(list_posterior_simul_by_gene, `[[`, 'PDF_NS')) %>% ungroup()

df_PDF_N_CV_2 <- Reduce('bind_rows', lapply(list_posterior_simul_by_gene_CV_2, `[[`, 'PDF_N')) %>% ungroup()
df_PDF_S_CV_2 <- Reduce('bind_rows', lapply(list_posterior_simul_by_gene_CV_2, `[[`, 'PDF_S')) %>% ungroup()
df_PDF_NS_CV_2 <- Reduce('bind_rows', lapply(list_posterior_simul_by_gene_CV_2, `[[`, 'PDF_NS')) %>% ungroup()


plot_hist_comp_nb_branches_S_N <- function(gene_to_plot, df_PDF_S, df_PDF_N){
  
  # Function to make a split legend with two colours
  draw_key_diag_split <- function(data, params, size) {
    gTree(children = gList(
      polygonGrob(
        x = c(0, 1, 1),
        y = c(1, 1, 0),
        gp = gpar(fill = col_non_synonymous, col = NA)
      ),
      polygonGrob(
        x = c(0, 0, 1),
        y = c(1, 0, 0),
        gp = gpar(fill = col_synonymous, col = NA)
      ),
      rectGrob(gp = gpar(fill = NA, col = NA, lwd = 0.8))
    ))
  }
  
  plt_hist <- df_PDF_S %>% filter(gene == gene_to_plot) %>% rename(n_aa_mut = S) %>% mutate(type = 'S') %>% 
    bind_rows(df_PDF_N %>% filter(gene == gene_to_plot) %>% rename(n_aa_mut = N) %>% mutate(type = 'N')) %>% 
    mutate(type = ifelse(type == 'S', 'Synonymous', 'Non-synonymous')) %>% 
    mutate(type = factor(type, levels = c('Synonymous', 'Non-synonymous'))) %>% 
    ggplot(aes(x = n_aa_mut)) +
    geom_bar(stat = 'identity', aes(y = n_branches_true, color = type, fill = type, alpha = "Data"), key_glyph = draw_key_diag_split) +
    geom_point(aes(y = mean_n_branches_new, color = 'Posterior draws', alpha = "Posterior draws")) +
    geom_linerange(aes(ymin = lower_n_branches_new, ymax = upper_n_branches_new, 
                       color = 'Posterior draws', alpha = "Posterior draws")) +
    facet_grid(. ~ type, scales = "free_x", space = "free_x") +
    scale_x_continuous(name = 'Number of AA mutations', breaks = 0:11) +
    scale_y_continuous(name = 'Number of branches', breaks = seq(0, 6e6, 1e6),
                       labels = c(0, expression(paste('1\u00B7', 10^6)), expression(paste('2\u00B7', 10^6)), expression(paste('3\u00B7', 10^6)), expression(paste('4\u00B7', 10^6)),
                                  expression(paste('5\u00B7', 10^6)), expression(paste('6 ', 10^6))),
                       expand = expansion(mult = c(0.02, 0.05))) +
    scale_fill_manual(name = '', values = c(col_synonymous, col_non_synonymous), guide = 'none') +
    scale_colour_manual(name = '', values = c(col_synonymous, col_non_synonymous, 'gray22'), guide = 'none') +
    scale_alpha_manual(name = "",
      values = c("Data" = 1, "Posterior draws" = 1),
      breaks = c("Data", "Posterior draws"),
      guide = guide_legend(
        override.aes = list(
          shape    = c(22, 16),
          linetype = c(0, 1),
          size     = c(6, 2)
        )
      )
    ) +
    theme_classic() +
    theme(legend.position = c(0.8, 0.8),
          legend.background = element_blank(),
          strip.background = element_rect(colour = 'white'),
          strip.text = element_text(size = 11))
  
  return(plt_hist)
}

posterior_predictive_checks_S_N_Spike <- plot_hist_comp_nb_branches_S_N('S', df_PDF_S, df_PDF_N)

pdf('../model_illustration/posterior_predictive_checks_S_N_Spike_CV_1.5.pdf', height = 2.8, width = 6.)
plot(posterior_predictive_checks_S_N_Spike)
dev.off()






value_zero <- 5e-2

df_PDF_NS %>% filter(gene != 'Genome') %>% mutate(CV = as.character(CV)) %>% 
  bind_rows(df_PDF_NS_CV_2 %>% filter(gene == 'Genome')) %>% 
  mutate(median_n_branches_new_crop = 
           ifelse(median_n_branches_new == 0, value_zero, median_n_branches_new),
         lower_n_branches_new = ifelse(lower_n_branches_new == 0, value_zero, lower_n_branches_new),
         upper_n_branches_new = ifelse(upper_n_branches_new == 0, value_zero, upper_n_branches_new),
         n_branches_new_crop = ifelse(n_branches_true == 0, value_zero, n_branches_true)
  ) %>% 
  ggplot() +
  geom_abline(intercept = 0., slope = 1., color = 'darkgrey', linetype = 'dashed') +
  geom_point(aes(x = n_branches_new_crop, y = median_n_branches_new_crop, colour = gene)) +
  geom_linerange(aes(x = n_branches_new_crop, ymin = lower_n_branches_new, ymax = upper_n_branches_new, color = gene)) +
  geom_rect(xmin = -Inf, xmax = +Inf,
            ymin = log(value_zero + 0.03), ymax = log(value_zero + 0.3),
            fill = 'white') +
  geom_rect(ymin = -Inf, ymax = +Inf,
            xmin = log(value_zero + 0.03), xmax = log(value_zero + 0.5),
            fill = 'white') +
  #annotate(geom = 'text', x = 5, y = 2e6, , label = paste0('R = ', round(df_cor_by_CV$cor_log_pearson[df_cor_by_CV$CV == 1.5], 3))) +
  scale_x_continuous(trans = 'log', 
                     breaks = c(value_zero, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7),
                     labels = c(0, expression(10^{0}), expression(10^{1}), expression(10^{2}),
                                expression(10^{3}), expression(10^{4}), expression(10^{5}),
                                expression(10^{6}), expression(10^{7})),
                     name = 'True number of branches') +
  scale_y_continuous(trans = 'log', 
                     breaks = c(value_zero, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7),
                     labels = c(0, expression(10^{0}), expression(10^{1}), expression(10^{2}),
                                expression(10^{3}), expression(10^{4}), expression(10^{5}),
                                expression(10^{6}), expression(10^{7})),
                     name = 'Predicted number of branches\n(Median from posterior)') +
  scale_colour_manual(name = 'Genome region', breaks = names(palette_genes), values = palette_genes) +
  coord_fixed() +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        panel.grid.minor = element_blank())
