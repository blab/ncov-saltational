library(tidyverse)

# Load Nextclade QC metrics on terminal saltational branches
nextclade <- read_tsv('../../results/nextclade/nextclade_saltational_branches_CV_1.8_BF_20.tsv')|> 
  select(seqName, totalSubstitutions, totalDeletions, 
         qc.overallScore, qc.overallStatus, 
         qc.missingData.score, qc.missingData.status,
         qc.mixedSites.score, qc.mixedSites.status,
         qc.privateMutations.score, qc.privateMutations.status,
         qc.snpClusters.score, qc.snpClusters.status,
         qc.frameShifts.score, qc.frameShifts.status,
         qc.stopCodons.score, qc.stopCodons.status) |> 
  mutate(qc.overall_score_no_private_mutations = qc.overallScore - ((qc.privateMutations.score)^2) / 100) |> 
  mutate(prop_qc_score_attributable_to_private_mut = ((qc.privateMutations.score)^2) / 100 / qc.overallScore)

# Display results
col_good <- 'forestgreen'
col_mediocre <- 'orange2'
col_bad <- 'firebrick3'


panel_terminal_salt_branches_quality_control_nextclade <- 
  nextclade |> 
  group_by(qc.privateMutations.status) |> summarise(n_branches = n()) |> 
  mutate(status = factor(qc.privateMutations.status, levels = c('good', 'mediocre', 'bad'))) |> 
  select(status, n_branches) |> complete(status, fill = list(n_branches = 0)) |> 
  mutate(qc_type = 'Private mutations') |> 
  bind_rows(nextclade |> 
              group_by(qc.mixedSites.status) |> summarise(n_branches = n()) |> 
              mutate(status = factor(qc.mixedSites.status, levels = c('good', 'mediocre', 'bad'))) |> 
              select(status, n_branches) |> complete(status, fill = list(n_branches = 0))|> 
              mutate(qc_type = 'Mixed sites')) |> 
  bind_rows(nextclade |> 
              group_by(qc.missingData.status) |> summarise(n_branches = n()) |> 
              mutate(status = factor(qc.missingData.status, levels = c('good', 'mediocre', 'bad'))) |> 
              select(status, n_branches) |> complete(status, fill = list(n_branches = 0)) |> 
              mutate(qc_type = 'Missing data')) |> 
  bind_rows(nextclade |> 
              group_by(qc.stopCodons.status) |> summarise(n_branches = n()) |> 
              mutate(status = factor(qc.stopCodons.status, levels = c('good', 'mediocre', 'bad'))) |> 
              select(status, n_branches) |> complete(status, fill = list(n_branches = 0)) |> 
              mutate(qc_type = 'Stop codons')) |> 
  bind_rows(nextclade |> 
              group_by(qc.frameShifts.status) |> summarise(n_branches = n()) |> 
              mutate(status = factor(qc.frameShifts.status, levels = c('good', 'mediocre', 'bad'))) |> 
              select(status, n_branches) |> complete(status, fill = list(n_branches = 0)) |> 
              mutate(qc_type = 'Frameshifts')) |> 
  bind_rows(nextclade |> 
              group_by(qc.snpClusters.status) |> summarise(n_branches = n()) |> 
              mutate(status = factor(qc.snpClusters.status, levels = c('good', 'mediocre', 'bad'))) |> 
              select(status, n_branches) |> complete(status, fill = list(n_branches = 0)) |> 
              mutate(qc_type = 'Mutation clusters')) |> 
  mutate(qc_type = factor(qc_type, levels = c('Private mutations', 'Missing data',
                                              'Mixed sites', 'Stop codons',
                                              'Frameshifts', 'Mutation clusters'))) |> 
  ggplot(aes(x = status, y = n_branches)) +
  geom_bar(stat = 'identity', aes(fill = status)) +
  geom_text(aes(label = n_branches), vjust = - 0.3) +
  facet_wrap(. ~ qc_type) +
  scale_fill_manual(breaks = c('good', 'mediocre', 'bad'), 
                    values = c(col_good, col_mediocre, col_bad),
                    guide = 'none') +
  scale_colour_manual(breaks = c('good', 'mediocre', 'bad'), 
                      values = c(col_good, col_mediocre, col_bad),
                      guide = 'none') +
  scale_x_discrete(name = 'QC status', 
                   breaks = c('good', 'mediocre', 'bad'), 
                   labels = c('Good', 'Mediocre', 'Bad')) +
  scale_y_continuous(name = 'Number of saltational branches',
                     expand = expansion(mult = c(0., 0.1)),
                     breaks = seq(0, 1500, 500)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = 'gray22'),
        strip.text = element_text(colour = 'white', size = 11))

png('../figure_mutation_enrichment/panel_quality_control_saltational_branches.png',
    height = 4.5, width = 6., res = 350, units = 'in')
plot(panel_terminal_salt_branches_quality_control_nextclade)
dev.off()

# Compare branches classified as saltational when removing all mutations flagged by nextclade as 
# mutation clusters, frameshifts, after stop codons or mixed sites
vec_saltational_branches_new <- unlist(read.table('../../results/qc_check/name_saltational_branches_qc_filtered_CV_1.8_BF_20.tsv'))
length(vec_saltational_branches_new)
nrow(nextclade)
