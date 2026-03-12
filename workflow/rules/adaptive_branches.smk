"""
Adaptive branch identification and sequence preparation.

This module handles:
- Identifying branches with significant Bayes Factors
- Concatenating adaptive branches across all genes
- Separating internal nodes vs tip sequences
- Reconstructing ancestral sequences for internal nodes
- Downloading tip sequences from Viridian
- Preparing final alignments and metadata for adaptive branches

Rules:
    - get_name_adaptive_branches: Identify adaptive branches per gene
    - concat_name_adaptive_branches: Combine adaptive branches across genes
    - get_name_internal_nodes_adaptive: Separate internal from tip nodes
    - ancestral_sequence_reconstruction: Reconstruct sequences for internal nodes
    - get_viridian_alignment_adaptive_branches: Download tip sequences
    - get_alignment_all_adaptive_branches: Combine all adaptive sequences and get alignment
    - get_metadata_all_adaptive_branches: Prepare metadata with dates from adaptive branches
"""


rule get_name_adaptive_branches:
    """Get the names of the adaptive paths across the tree for each gene and each CV prior"""
    input:
        data_viridian = rules.prepare_data_model.output.data_viridian,
        data_viridian_grouped = rules.prepare_data_model.output.data_viridian_grouped,
        df_BF = rules.get_BF_across_model_runs.output.df_BF
    output:
        name_adaptive_branches = "results/adaptive_branches/name_adaptive_branches_{gene}_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        tsv-join -H \
            --filter-file <(
                tsv-filter -H --gt BF_omega_greater_1:{wildcards.BF_threshold} {input.df_BF} |
                tsv-join -H --filter-file {input.data_viridian_grouped} --key-fields i_row
            ) \
            --key-fields synonymous_muts,non_synonymous_muts,syn_sites,nonsyn_sites {input.data_viridian} |
            tsv-select -H -f strain |
            sed 1d > {output.name_adaptive_branches}
        """


rule concat_name_adaptive_branches:
    """Concatenate adaptive branch names across all genes"""
    input:
        name_adaptive_branches_across_genes = lambda wildcards: expand("results/adaptive_branches/name_adaptive_branches_{gene}_CV_{curr_cv_prior}_BF_{curr_BF_threshold}.tsv", gene=GENES, curr_cv_prior = wildcards.cv_prior, curr_BF_threshold = wildcards.BF_threshold),
    output:
        concat_name_adaptive_branches = "results/adaptive_branches/concat_name_adaptive_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        tsv-append -H {input.name_adaptive_branches_across_genes} | tsv-uniq >> {output.concat_name_adaptive_branches}
        """


rule get_name_internal_nodes_adaptive:
    """Identify which adaptive branches are internal nodes (vs tips)"""
    input:
        name_adaptive_branches = rules.concat_name_adaptive_branches.output.concat_name_adaptive_branches,
    output:
        name_internal_nodes_adaptive = "results/adaptive_branches/name_adaptive_internal_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        tsv-filter -H --regex 'strain:^node_' <(printf "strain\n"; cat {input.name_adaptive_branches}) | sed 1d > {output.name_internal_nodes_adaptive}
        """


rule ancestral_sequence_reconstruction:
    """Perform ancestral sequence reconstruction for internal adaptive branches"""
    input:
        df_child_parent_relationships = rules.get_parent_child_relationship.output.df_parent_relationships,
        name_internal_nodes_adaptive = rules.get_name_internal_nodes_adaptive.output.name_internal_nodes_adaptive,
        df_nuc_mut = rules.char_nuc_across_tree.output.df_nuc_mut,
        ref_fasta = rules.params.input.usher_reference_fasta,
    output:
        output_fasta = "results/alignment_adaptive_branches/alignment_adaptive_internal_nodes_CV_{cv_prior}_BF_{BF_threshold}.fasta"
    shell:
        """
        Rscript ./scripts/mat-manipulation/reconstruct_ancestral_seq.R \
            --df_nuc_mut {input.df_nuc_mut} \
            --ref_fasta {input.ref_fasta} \
            --nodes_to_reconstruct {input.name_internal_nodes_adaptive} \
            --df_child_parent {input.df_child_parent_relationships} \
            --output_ancestral_fasta {output.output_fasta}
        """


rule get_viridian_alignment_adaptive_branches:
    """Download sequences for adaptive branch tips from Viridian"""
    input:
        name_adaptive_branches = rules.concat_name_adaptive_branches.output.concat_name_adaptive_branches,
    output:
        alignment_adaptive_final_nodes = "results/alignment_adaptive_branches/alignment_adaptive_final_nodes_CV_{cv_prior}_BF_{BF_threshold}.fasta",
    params:
        url_viridian_1 = "https://api.figshare.com/v2/file/download/45969777",
        url_viridian_2 = "https://api.figshare.com/v2/file/download/49692480",
        file_path_viridian_1 = "data/tmp/viridian_alignment_1.tar",
        file_path_viridian_2 = "data/tmp/viridian_alignment_2.tar",
        dir_unpacked_viridian = "data/tmp/",
    shell:
        """
        Rscript ./scripts/mat-manipulation/download_viridian_alignment.R \
            {params.url_viridian_1} {params.url_viridian_2} {params.file_path_viridian_1} {params.file_path_viridian_2} \
            {params.dir_unpacked_viridian} {params.dir_unpacked_viridian} \
            {input.name_adaptive_branches} {output.alignment_adaptive_final_nodes}
        """

rule get_alignment_all_adaptive_branches:
    """Combine internal node and tip sequences for all adaptive branches"""
    input:
        alignment_adaptive_final_nodes = rules.get_viridian_alignment_adaptive_branches.output.alignment_adaptive_final_nodes,
        alignment_adaptive_internal_nodes = rules.ancestral_sequence_reconstruction.output.output_fasta,
    output:
        alignment_all_adaptive_branches = "results/alignment_adaptive_branches/alignment_adaptive_branches_CV_{cv_prior}_BF_{BF_threshold}.fasta"
    shell:
        """
        cat {input.alignment_adaptive_final_nodes} {input.alignment_adaptive_internal_nodes} > {output.alignment_all_adaptive_branches}
        """

rule get_metadata_all_adaptive_branches:
    """Prepare metadata for all adaptive branches with inferred dates"""
    input:
        chronumental_inferred_dates = rules.chronumental.output.inferred_dates,
        name_adaptive_branches = rules.concat_name_adaptive_branches.output.concat_name_adaptive_branches,
        metadata = rules.get_metadata_from_taxonium.output.metadata,
    output:
        metadata = "results/alignment_adaptive_branches/metadata_all_adaptive_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        Rscript ./scripts/mat-manipulation/get_dates_adaptive_branches.R \
            --chronumental {input.chronumental_inferred_dates} \
            --names_adapt_branches {input.name_adaptive_branches} \
            --metadata {input.metadata} \
            --output_metadata {output.metadata}
        """
