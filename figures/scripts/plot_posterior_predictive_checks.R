library(tidyverse)
library(grid)
library(ggpubr)
source('config_plot.R')

## Load posterior predictive checks runs
get_summary_from_posterior_simul_PDF_NS <- function(gene, CV_prior){
  # Load posterior simulations
  if(gene == 'Genome'){
    vec_i_simul <- 1:100
    df_posterior_simul <- Reduce('bind_rows', lapply(vec_i_simul, FUN = function(i_simul){
      read_tsv(paste0('/Volumes/homes/viridian-dnds/hpc_run/results/posterior_predictive_checks_genome/posterior_predictive_check_', gene, '_CV_', CV_prior, '_isimul_', i_simul, '.tsv'))
    }))
  } else{
    #df_posterior_simul <- read_tsv(paste0('/Volumes/homes/viridian-dnds/hpc_run/results/posterior_predictive_checks/posterior_predictive_check_', gene, '_CV_', CV_prior, '.tsv'))
    df_posterior_simul <- read_tsv(paste0('../../results/posterior_predictive_checks/posterior_predictive_check_', gene, '_CV_', CV_prior, '.tsv'))
  }
  
  df_posterior_simul_count_NS <- df_posterior_simul %>% group_by(S_new, N_new, i_simul) %>% 
    summarise(n_branches_new = sum(n_branches_new)) %>% 
    ungroup() %>% 
    complete(S_new, N_new, i_simul, fill = list(n_branches_new = 0))
  
  # Summary of counts across posterior simulations
  PDF_NS <- df_posterior_simul_count_NS %>% 
    group_by(S_new, N_new) %>% 
    summarise(mean_n_branches_new = mean(n_branches_new),
              median_n_branches_new = median(n_branches_new),
              lower_n_branches_new = quantile(n_branches_new, 0.025),
              upper_n_branches_new = quantile(n_branches_new, 0.975),
              min_n_branches_new = min(n_branches_new),
              max_n_branches_new = max(n_branches_new)) 
  
  # Load actual data
  if(gene == 'Genome'){
    gene_data <- read_tsv(paste0('../../results/input_model/data_viridian_', gene, '_aggregated.tsv')) %>% 
      group_by(synonymous_muts, non_synonymous_muts) %>% summarise(n_branches = n())
  } else{
    gene_data <- read_tsv(paste0('../../results/input_model/data_viridian_', gene, '_grouped.tsv')) 
  }
  
  df_true_count_SN <- gene_data %>% 
    group_by(synonymous_muts, non_synonymous_muts) %>% summarise(n_branches_true = sum(n_branches)) %>% 
    rename(S = synonymous_muts, N = non_synonymous_muts) %>% 
    group_by(S, N) %>% summarise(n_branches_true = sum(n_branches_true))
  
  
  # Merge PDF
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
    PDF_NS %>% mutate(gene = gene, CV = CV_prior)
  )
}

vec_CV <- sprintf("%.1f", c(1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0))
vec_genes <- c('ORF1a', 'ORF1b', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF9b')

df_PDF_NS <- Reduce('bind_rows', lapply(vec_CV, FUN = function(CV_prior){
  Reduce('bind_rows', lapply(vec_genes, FUN = function(gene){
    get_summary_from_posterior_simul_PDF_NS(gene, CV_prior) 
  }))
}))

CV_maximizing_cor <- df_PDF_NS %>% 
  group_by(CV) %>% 
  summarise(cor = cor(log(n_branches_true + 1), log(median_n_branches_new + 1), method = 'pearson')) %>% 
  filter(cor == max(cor)) %>% select(CV) %>% unlist() %>% as.numeric()

df_PDF_NS_Genome <- Reduce('bind_rows', lapply(CV_maximizing_cor, FUN = function(CV_prior){
  Reduce('bind_rows', lapply('Genome', FUN = function(gene){
    get_summary_from_posterior_simul_PDF_NS(gene, CV_prior) 
  }))
}))

## Illustrate impact of CV on model fit
plot_cor_func_CV <- function(df_PDF_NS){
  plt_cor_func_CV <- df_PDF_NS %>%
    mutate(CV = as.numeric(CV),
           is_genome = (gene == 'Genome')) %>% 
    group_by(CV, is_genome) %>% 
    summarise(cor = cor(log(n_branches_true + 1), 
                        log(median_n_branches_new + 1), 
                        method = 'pearson')) %>% 
    arrange(is_genome) %>% 
    ggplot(aes(x = CV, y = cor, colour = is_genome)) +
    geom_point() +
    scale_x_continuous(name = expression(CV[prior]),
                       breaks = seq(1., 8., 1.)) +
    scale_y_continuous(name = 'Pearson correlation coefficient',
                       #limits = c(0.9, 1.)
                       ) +
    scale_color_manual(name = '', breaks = c(F, T), values = c('gray22', 'darkgrey'), labels = c('Genes', 'Genome')) +
    theme_bw() +
    theme(legend.background = element_blank(),
          legend.position = 'none',
          legend.position.inside = c(0.8, 0.2))
  
  return(plt_cor_func_CV)
}
plot_predicted_number_branches_func_observed_func_genes_CV <- function(df_PDF_NS, CVs_to_plot){
  
  value_zero <- 5e-2
  
  df_for_plot <- df_PDF_NS %>% mutate(CV = as.numeric(CV)) %>% 
    mutate(is_genome = ifelse(gene == 'Genome', 'Genome', 'Genes')) %>% 
    filter(CV %in% CVs_to_plot) %>%
    mutate(CV_label = paste0('CV[prior] == ', CV))
  
  df_cor_for_plot <- df_for_plot %>% 
    group_by(is_genome, CV_label) %>% 
    summarise(cor = cor(log(n_branches_true + 1), log(median_n_branches_new + 1), method = 'pearson'))
  
  panel_plt <- df_for_plot  %>% 
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
    facet_grid(is_genome ~ CV_label, labeller = label_parsed) +
    geom_text(data = df_cor_for_plot, aes(label = paste0('R = ', round(cor, 3))), 
              x = value_zero, y = log(1e6)) +
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
          panel.grid.minor = element_blank(),
          legend.position = 'bottom',
          strip.background = element_rect(fill = 'gray22'),
          strip.text = element_text(colour = 'white', size = 11))
  
  return(panel_plt)
}
get_CV_maximizing_cor <- function(df_PDF_NS){
  df_PDF_NS %>%
    mutate(CV = as.numeric(CV), is_genome = (gene == 'Genome')) %>% 
    group_by(CV, is_genome) %>% 
    summarise(cor = cor(log(n_branches_true + 1), log(median_n_branches_new + 1), 
                        method = 'pearson')) %>% 
    group_by(is_genome) %>% 
    filter(cor == max(cor)) %>% 
    ungroup() %>% 
    return()
}

panel_across_genes_CV <- plot_predicted_number_branches_func_observed_func_genes_CV(df_PDF_NS, 
                                                                                   CVs_to_plot = c(1.2, 1.5, 2.0))
plt_cor_func_CV <- plot_cor_func_CV(df_PDF_NS)

df_max_cor <- get_CV_maximizing_cor(df_PDF_NS)

panel_impact_CV <- ggarrange(panel_across_genes_CV, 
                             plt_cor_func_CV + theme(aspect.ratio = 0.8), 
                             nrow = 2, labels = 'AUTO', heights = c(1., 0.5))

png('../model_illustration/plt_impact_CV_PPC.png', height = 10., width = 10, res = 350, units = 'in')
plot(panel_impact_CV)
dev.off()

# Load posterior predictive checks
list_posterior_simul_by_gene_S <- lapply(vec_genes, FUN = function(gene){
  get_summary_from_posterior_simul(gene, '1.7') 
})

df_PDF_N_CV_1.9 <- Reduce('bind_rows', lapply(list_posterior_simul_by_gene_S, `[[`, 'PDF_N')) %>% ungroup()
df_PDF_S_CV_1.9 <- Reduce('bind_rows', lapply(list_posterior_simul_by_gene_S, `[[`, 'PDF_S')) %>% ungroup()

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
    scale_x_continuous(name = 'Number of mutations', breaks = 0:11) +
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

posterior_predictive_checks_S_N_Spike <- plot_hist_comp_nb_branches_S_N('S', df_PDF_S_CV_1.9, df_PDF_N_CV_1.9)

pdf('../model_illustration/posterior_predictive_checks_S_N_Spike_CV_1.9.pdf', 
    height = 2.8, width = 6.)
plot(posterior_predictive_checks_S_N_Spike)
dev.off()

## Make joint plot
plot_aggregated_predicted_number_branches_func_observed <- function(df_PDF_NS, CV_genome_to_plot, CV_genes_to_plot){
  value_zero <- 5e-2
  
  
  df_plot <- df_PDF_NS %>% 
    filter(((gene == 'Genome') & (CV == CV_genome_to_plot)) | ((gene != 'Genome') & (CV == CV_genes_to_plot))) %>% 
    arrange(gene != 'Genome') %>% 
    mutate(median_n_branches_new_crop = 
             ifelse(median_n_branches_new == 0, value_zero, median_n_branches_new),
           lower_n_branches_new = ifelse(lower_n_branches_new == 0, value_zero, lower_n_branches_new),
           upper_n_branches_new = ifelse(upper_n_branches_new == 0, value_zero, upper_n_branches_new),
           n_branches_new_crop = ifelse(n_branches_true == 0, value_zero, n_branches_true)
    )
    
  df_cor <- df_plot %>% 
    ungroup() %>% 
    summarise(cor = cor(log(n_branches_true + 1), log(median_n_branches_new + 1), 
                        method = 'pearson'))
  
  scatter_plot_branch_counts <- df_plot %>% 
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
    annotate(geom = 'text', x = 5, y = 2e6, , label = paste0('R = ', round(df_cor$cor, 3))) +
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
  
  return(scatter_plot_branch_counts)
}

CV_genome_to_plot <- '6.8'
CV_genes_to_plot <- '1.9'

scatter_plot_comp <- plot_aggregated_predicted_number_branches_func_observed(df_PDF_NS, CV_genome_to_plot, CV_genes_to_plot)

pdf('../model_illustration/scatter_plot_posterior_predictive_checks_main.pdf', 
height = 3.8, width = 4.8)
plot(scatter_plot_comp)
dev.off()
