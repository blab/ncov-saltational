library(RColorBrewer)
library(paletteer)

# Script containing color used for the visualization

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

## Palette for NSPs
palette_NSPs <- c(brewer.pal(12, 'Set3'), 'black', 'gray55', 'deeppink3')
names(palette_NSPs) <- paste0('NSP', c(1:10, 12:16))

## Colour for synonymous and non-synonymous
col_synonymous <- 'darkcyan'
col_non_synonymous <- 'firebrick2'

## Colour for posterior and prior
col_prior <- 'hotpink3'
col_posterior <-'black'

## Colour for adaptive and non-adaptive branches
col_adaptive <- 'orange2'
col_non_adaptive <- 'darkgrey'

## Colour for recurrent mutations
col_persistent <- '#fa4d56'
col_minks <- '#6fdc8c'
col_deers <- '#33b1ff'
col_others <- 'gray80'

# Colour for geographic regions
palette_africa <- paletteer_dynamic(palette = "cartography::blue.pal", n = 20)
palette_asia <- paletteer_dynamic(palette = "cartography::red.pal", n = 20)
palette_north_america <- paletteer_dynamic(palette = "cartography::sand.pal", n = 20)
palette_south_america <- paletteer_dynamic(palette = "cartography::brown.pal", n = 20)
palette_europe <- paletteer_dynamic(palette = "cartography::green.pal", n = 20)
palette_oceania <- paletteer_dynamic(palette = "cartography::purple.pal", n = 20)