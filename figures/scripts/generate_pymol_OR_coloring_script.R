# Generate a PyMOL script that colors Spike by residue-level odds ratio
# after filtering the input data to gene == "S" and
# threshold_for_inclusion == 2.

# ---- Input and output files ----
input_tsv <- "../../results/mutation_enrichment/df_AA_pos_enrichment_saltational_CV_1.8_BF_20.tsv"
output_pml <- "plot_7KRQ_OR_surface.pml"

# ---- Read and filter the data ----
filtered <- read_tsv(input_tsv) %>% 
  select(OR, gene, position, threshold_for_inclusion) %>% 
  filter(gene == 'S', threshold_for_inclusion == 1) %>% 
  arrange(position) %>% 
  mutate(log_OR = log(OR))


# ---- Build the PyMOL script ----
pml <- c(
  "set assembly, 1",
  "fetch 7KRQ, async=0",
  "hide everything, 7KRQ",
  "show surface, 7KRQ",
  "color gray80, 7KRQ",
  "set seq_view, 1",
  "bg_color white",
  "set surface_quality, 2",
  "set antialias, 2",
  "set two_sided_lighting, on",
  "set ambient, 0.45",
  "set direct, 0.35",
  "set specular, 0.15",
  "set shininess, 10",
  "alter 7KRQ, b=-999.0"
)

# Assign each residue's log(OR) to all atoms in that residue.
assignment_lines <- sprintf(
  "alter 7KRQ and resi %d, b=%.8f",
  filtered$position,
  filtered$log_OR
)

pml <- c(
  pml,
  assignment_lines,
  "rebuild",
  "spectrum b, blue_white_red, 7KRQ and b>-900, minimum=-3, maximum=3",
  "orient 7KRQ",
  "rotate z, -90",
  "rotate y, 30",
  "ray 2400,2400",
  "png 7KRQ_OR_surface.png, dpi=300"
)
writeLines(pml, output_pml)

message(
  "Wrote ", output_pml,
  " with ", nrow(filtered),
  " unique Spike residue position(s)."
)

