library(RColorBrewer)

# Palette definitions
## Genome regions
palette_genes <- c(
  'ORF1a' = brewer.pal(12, 'Paired')[1],
  'ORF1b' = brewer.pal(12, 'Paired')[2],
  'S' =  'gray22',
  'ORF3a' = brewer.pal(12, 'Paired')[3],
  'E' = brewer.pal(12, 'Paired')[4],
  'M' = brewer.pal(12, 'Paired')[5],
  'ORF6' = brewer.pal(12, 'Paired')[6],
  'ORF7a' = brewer.pal(12, 'Paired')[7],
  'ORF7b' = brewer.pal(12, 'Paired')[8],
  'ORF8' = brewer.pal(12, 'Paired')[9],
  'N' = brewer.pal(12, 'Paired')[10],
  'ORF9b' = brewer.pal(12, 'Paired')[12],
  'Genome' = 'grey80')

## Spike regions
palette_spike <- c(
  'NTD' = 'deeppink3',
  'RBD' = 'cyan3',
  'Other S1' = 'darkslateblue',
  'S2' = 'orange2',
  'Other S' = 'black'
)
