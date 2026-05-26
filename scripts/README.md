# Scripts

Custom scripts used by the Snakemake workflow. Scripts are organized by analysis stage and called from the rule files in `workflow/rules/`.

---

## `mat-manipulation/`

Scripts for processing the Viridian UShER/MAT phylogeny (protobuf format) and extracting per-branch information.

| Script | Description |
|---|---|
| `get_metadata_from_taxonium.py` | Extracts sample metadata (date, country, clade) from a Taxonium JSONL file. |
| `get_nuc_mutation_across_genome.py` | Enumerates nucleotide mutations on every branch of the MAT tree. |
| `get_clades_per_node.py` | Assigns Nextstrain/Pango clade labels to each node in the tree. |
| `get_child_parent_relationship.py` | Extracts parent–child relationships for all nodes in the tree. |
| `get_n_descendants_per_node.py` | Counts total and terminal descendants per node, with earliest/latest tip dates. |
| `get_mutation_opportunity.py` | Computes the number of synonymous and non-synonymous sites available per branch per gene (mutation opportunity). |
| `count_sites_per_branch.py` | Counts synonymous and non-synonymous sites per branch using codon-level logic. |
| `get_stop_codons_across_tree.py` | Identifies branches carrying premature stop codons across the tree. |
| `reconstruct_ancestral_location.py` | Infers ancestral geographic location (country/continent) by parsimony on the tree. |
| `get_branches_leading_to_recombinant_clades.R` | Identifies branches leading to known recombinant clades to exclude from analysis. |
| `get_non_coded_positions.R` | Identifies nucleotide positions that fall outside annotated coding sequences. |
| `tidy_node_translations.R` | Reformats per-node amino acid translation output from matUtils into a tidy TSV. |
| `reconstruct_ancestral_seq.R` | Downloads ancestral sequences for internal nodes from the Viridian alignment. |
| `download_viridian_alignment.R` | Downloads per-sample sequences from the Viridian alignment for a set of sample IDs. |
| `get_metadata_for_taxonium_export.R` | Prepares metadata for export to Taxonium for interactive visualization. |
| `get_terminal_descendants_of_saltational_branches.py` | Lists all terminal (tip) descendants of each saltational branch. |
| `get_cum_aa_changes_from_root.py` | Computes the cumulative number of amino acid changes from the root to each node. |

---

## `bayesian-model/`

Scripts implementing the Bayesian hierarchical model for detecting adaptive evolution (elevated dN/dS) per branch per gene.

| Script | Description |
|---|---|
| `prepare_data_bayesian_model.R` | Joins mutation counts and site counts into the input format expected by the model. |
| `get_prior_parametrization.R` | Estimates Gamma prior parameters (mean dN, mean dS) from tree-wide statistics. |
| `run_bayesian_model_by_gene.R` | Runs the Bayesian model per gene, returning posterior estimates of dN and dS per branch. |
| `run_bayesian_model_genome.R` | Runs the genome-wide version of the Bayesian model. |
| `aggregate_genome_data.R` | Aggregates per-gene mutation and site data into a genome-wide summary. |
| `compute_BF_across_branches.R` | Computes Bayes Factors (BF) for each branch from the posterior dN/dS estimates. |
| `run_posterior_predictive_checks.R` | Runs posterior predictive checks for a single gene across a range of CV prior values. |
| `run_posterior_predictive_checks_genome.R` | Runs posterior predictive checks for the genome-wide model. |
| `utils_model.R` | Utility functions shared by the Bayesian model scripts (Gamma parametrization, model object). |

---

## `branch-analysis/`

| Script | Description |
|---|---|
| `get_saltational_branches.R` | Classifies adaptive branches (BF > threshold) as saltational if their non-synonymous mutation count exceeds a percentile threshold across all branches. |

---

## `mutation-enrichment/`

Scripts for computing enrichment of amino acid mutations and positions in saltational branches relative to background branches.

| Script | Description |
|---|---|
| `calculate_AA_mut_enrichment.R` | Computes odds ratios for each AA mutation being enriched in saltational vs. background branches. |
| `calculate_AA_pos_enrichment.R` | Computes odds ratios for each AA position being enriched in saltational vs. background branches. |
| `calculate_AA_mut_enrichment_account_background.R` | Computes mutation enrichment odds ratios with correction for background mutation opportunity. |
| `get_AA_reference_sequence.py` | Extracts the amino acid reference sequence for each gene from the SARS-CoV-2 reference. |
| `get_AA_frequency_nextstrain_tree.py` | Computes the frequency of each amino acid state at each position across the Nextstrain global tree. |

---

## `sra-metadata/`

Scripts for downloading and curating NCBI SRA/BioSample metadata for all samples in the Viridian tree.

| Script | Description |
|---|---|
| `extract_biosample_ids.py` | Maps SRA run accessions to BioSample IDs using the NCBI SRA_Accessions.tab file. |
| `parse_concatenated_biosample.py` | Parses batch-downloaded BioSample XML records into a flat TSV. |
| `tidy_final_metadata.R` | Merges and cleans all metadata sources; filters to human samples with valid dates. |
| `remove_incorrect_dates.R` | Removes samples with implausible collection dates (before pandemic onset, etc.). |
| `prepare_metadata_countries.R` | Standardizes country names and adds continent/region annotations. |

---

## `tree-proximity/`

Scripts for constructing a background subsample of sequences in phylogenetic proximity to saltational branches.

| Script | Description |
|---|---|
| `convert_tree_units_to_mutations_per_site.py` | Converts MAT tree branch lengths to mutations-per-site units for use with Augur. |
| `get_metadata_proximity_samples.R` | Compiles metadata for background samples selected based on phylogenetic proximity to saltational branches. |
