"""
Viridian tree processing, temporal inference, and mutation analysis.

This module handles:
- Converting UShER protobuf tree to Newick format
- Extracting metadata from Taxonium files
- Running Chronumental for temporal inference
- Extracting nucleotide mutations across the tree
- Getting amino acid translations per node
- Identifying stop codons and non-coded regions
- Characterizing tree topology (parent-child relationships, descendants)
- Counting synonymous and non-synonymous sites per branch

Rules:
    - convert_tree_newick: Convert tree format for Chronumental
    - get_metadata_from_taxonium: Extract and prepare metadata
    - chronumental: Infer dates using molecular clock
    - char_nuc_across_tree: Extract nucleotide mutations
    - get_translations_per_node: Get AA translations with matUtils
    - tidy_node_translations: Reformat translation data
    - get_stop_codons: Identify stop codons across tree
    - add_coded_status_aa_mut: Mark coded vs non-coded mutations
    - get_parent_child_relationship: Extract tree topology
    - get_n_descendants_per_node: Count descendants per node
    - count_sites_per_branch: Count syn/nonsyn sites per branch
"""

rule subset_viridian_tree:
    """
    Subset the Viridian tree to only include strains from the filtered metadata.

    This removes non-human samples (cell cultures, wastewater, etc.) identified
    by the tidy_final_metadata rule.
    """
    input:
        tree = rules.params.input.viridian_tree,
        strains_to_keep = rules.tidy_final_metadata.output.strains_to_keep
    output:
        subset_tree = "results/viridian/viridian_tree_subset.pb.gz"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        matUtils extract \
            --input-mat {input.tree} \
            --samples {input.strains_to_keep} \
            --write-mat {output.subset_tree}
        """


rule convert_tree_newick:
    """Convert UShER tree to newick format - This will be used to run Chronumental"""
    input:
        usher_tree = rules.subset_viridian_tree.output.subset_tree
    output:
        tree_newick = "results/viridian/viridian_tree.nwk"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        matUtils extract \
            --input-mat {input.usher_tree} \
            --write-tree {output.tree_newick}
        """

rule chronumental:
    """Run chronumental to estimate time tree"""
    input:
        tree = rules.convert_tree_newick.output.tree_newick,
        metadata = rules.clean_dates.output.metadata
    output:
        inferred_dates = "results/chronumental/chronumental_inferred_dates.tsv",
        timetree = "results/chronumental/chronumental_time_tree.pb.gz"
    params:
        reference_node = "SRR11597132"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        chronumental \
            --tree {input.tree} \
            --dates {input.metadata} \
            --dates_out {output.inferred_dates} \
            --tree_out {output.timetree} \
            --reference_node {params.reference_node} \
            --steps 2000
        """


rule char_nuc_across_tree:
    """List of nucleotide mutations across the UShER tree"""
    input:
        tree = rules.subset_viridian_tree.output.subset_tree,
    output:
        df_nuc_mut = "results/mutations/df_nuc_mut.tsv"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/get_nuc_mutation_across_genome.py \
            --tree {input.tree} \
            --output_df_nuc_mutations {output.df_nuc_mut}
        """


rule get_translations_per_node:
    """Get mutation information across tree"""
    input:
        mat = rules.subset_viridian_tree.output.subset_tree,
        reference_seq_gtf = rules.params.input.usher_reference_gtf,
        reference_seq_fasta = rules.params.input.usher_reference_fasta,
    output:
        node_translations = "results/mutations/nodes_translations.tsv"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        matUtils summary \
            --input-mat {input.mat} \
            --input-gtf {input.reference_seq_gtf} \
            --input-fasta {input.reference_seq_fasta} \
            --translate {output.node_translations}
        """

rule tidy_node_translations:
    """Tidy the translations from matUtils summary to a more easily workable format"""
    input:
        node_translations = rules.get_translations_per_node.output.node_translations
    output:
        aa_muts = "results/mutations/nodes_aa_muts_summary.tsv",
        nt_muts = "results/mutations/nodes_nt_muts_summary.tsv"
    shell:
        """
        Rscript ./scripts/mat-manipulation/tidy_node_translations.R \
            --nodes_translations {input.node_translations} \
            --output_nt_mutations {output.nt_muts} \
            --output_aa_mutations {output.aa_muts}
        """


rule get_stop_codons:
    """Get stop codons across tree for each node"""
    input:
        tree = rules.subset_viridian_tree.output.subset_tree,
        df_aa_muts = rules.tidy_node_translations.output.aa_muts,
    output:
        df_stop_codons = "results/bte/df_stop_codons.tsv"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/get_stop_codons_across_tree.py \
            --tree {input.tree} \
            --df_aa_muts {input.df_aa_muts} \
            --output_df_stop_codons {output.df_stop_codons}
        """

rule add_coded_status_aa_mut:
    """Add whether aa mutation is coded or not"""
    input:
        df_aa_muts = rules.tidy_node_translations.output.aa_muts,
        df_stop_codons = rules.get_stop_codons.output.df_stop_codons,
    output:
        df_aa_muts_coded = "results/mutations/nodes_aa_muts_summary_coded.tsv",
    shell:
        """
        Rscript ./scripts/mat-manipulation/get_non_coded_positions.R {input.df_aa_muts} {input.df_stop_codons} {output.df_aa_muts_coded}
        """


rule get_parent_child_relationship:
    """Extract parent-child relationships from tree topology"""
    input:
        tree = rules.subset_viridian_tree.output.subset_tree,
    output:
        df_parent_relationships = "results/bte/df_parent_child_relationship.tsv"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/get_child_parent_relationship.py \
            --tree {input.tree} \
            --output_parent_relationship {output.df_parent_relationships}
        """


rule get_n_descendants_per_node:
    """Get number of descendants per node"""
    input:
        tree = rules.subset_viridian_tree.output.subset_tree,
    output:
        df_n_descendants = "results/bte/df_n_descendants_per_node.tsv"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/get_n_descendants_per_node.py \
            --tree {input.tree} \
            --output_df_n_descendants {output.df_n_descendants}
        """

rule count_sites_per_branch:
    """Count the number of possible syn and nonsyn mutations available per gene on each branch of the MAT tree"""
    input:
        tree = rules.subset_viridian_tree.output.subset_tree,
        reference_seq_gtf = rules.params.input.usher_reference_gtf,
        reference_seq_fasta = rules.params.input.usher_reference_fasta,
    output:
        df_sites = "results/bte/df_sites_syn_non_syn_truncated_at_stop.tsv"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/count_sites_per_branch.py \
            --ref-fasta {input.reference_seq_fasta} \
            --ref-gtf {input.reference_seq_gtf} \
            --tree-path {input.tree} \
            --results-file {output.df_sites} \
            --truncate-at-stop
        """

rule reconstruct_ancestral_country:
    "Infer country of internal nodes in the viridian tree"
    input:
        tree = rules.subset_viridian_tree.output.subset_tree,
        metadata = rules.prepare_metadata_with_locations.output.metadata,
        country_mapping = rules.params.input.country_mapping
    output:
        inferred_countries = "results/bte/inferred_countries.tsv"
    params:
        location_col = "country"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/reconstruct_ancestral_location.py \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --country_mapping {input.country_mapping} \
            --location_col {params.location_col} \
            --output_locations {output.inferred_countries}
        """

rule reconstruct_ancestral_subregion:
    "Infer subregion of internal nodes in the viridian tree"
    input:
        tree = rules.subset_viridian_tree.output.subset_tree,
        metadata = rules.prepare_metadata_with_locations.output.metadata,
        country_mapping = rules.params.input.country_mapping
    output:
        inferred_subregions = "results/bte/inferred_subregions.tsv"
    params:
        location_col = "subregion"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/reconstruct_ancestral_location.py \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --country_mapping {input.country_mapping} \
            --location_col {params.location_col} \
            --output_locations {output.inferred_subregions}
        """

rule reconstruct_ancestral_continent:
    "Infer continent of internal nodes in the viridian tree"
    input:
        tree = rules.subset_viridian_tree.output.subset_tree,
        metadata = rules.prepare_metadata_with_locations.output.metadata,
        country_mapping = rules.params.input.country_mapping
    output:
        inferred_continents = "results/bte/inferred_continents.tsv"
    params:
        location_col = "continent"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/reconstruct_ancestral_location.py \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --country_mapping {input.country_mapping} \
            --location_col {params.location_col} \
            --output_locations {output.inferred_continents}
        """