"""
Tree building rules for adaptive branches with phylogenetic and temporal context.

Steps:
1. Align sequences to reference
2. Build maximum-likelihood tree with IQ-TREE
3. Refine tree with temporal information
4. Prepare node data for Auspice (branch_type annotations)
5. Export tree in Auspice JSON format
"""
rule get_name_internal_nodes_adaptive:
    """Identify which adaptive branches are internal nodes (vs tips)"""
    input:
        name_adaptive_branches = rules.concat_name_adaptive_branches.output.concat_name_adaptive_branches,
    output:
        name_internal_nodes_adaptive = "results/tree_building/name_adaptive_internal_branches_{BF_threshold}.tsv"
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
        output_fasta = "results/tree_building/alignment_adaptive_internal_nodes_BF_{BF_threshold}.fasta"
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
        alignment_adaptive_final_nodes = "results/tree_building/alignment_adaptive_final_nodes_BF_{BF_threshold}.fasta",
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
        alignment_all_adaptive_branches = "results/alignment_adaptive_branches/alignment_adaptive_branches_BF_{BF_threshold}.fasta"
    shell:
        """
        cat {input.alignment_adaptive_final_nodes} {input.alignment_adaptive_internal_nodes} > {output.alignment_all_adaptive_branches}
        """









rule get_metadata_all_adaptive_branches:
    """Prepare metadata for all adaptive branches with inferred dates"""
    input:
        chronumental_inferred_dates = rules.chronumental.output.inferred_dates,
        name_adaptive_branches = rules.concat_name_adaptive_branches.output.concat_name_adaptive_branches,
        metadata = rules.tidy_final_metadata.output.subset_metadata,
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

rule rename_internal_nodes_fasta:
    """Rename internal nodes in FASTA to avoid conflicts when adding to tree (node_XXX -> adaptive_node_XXX)"""
    input:
        internal_nodes = rules.ancestral_sequence_reconstruction.output.output_fasta,
    output:
        renamed_fasta = "results/proximity/internal_nodes_renamed_CV_{cv_prior}_BF_{BF_threshold}.fasta"
    shell:
        """
        # Rename sequence headers: node_XXX -> adaptive_node_XXX
        sed 's/>node_/>adaptive_node_/g' {input.internal_nodes} > {output.renamed_fasta}
        """

rule rename_internal_nodes_list:
    """Create renamed list of internal adaptive nodes for extraction"""
    input:
        adaptive_branches = rules.concat_name_adaptive_branches.output.concat_name_adaptive_branches,
    output:
        renamed_nodes = "results/proximity/internal_nodes_renamed_list_CV_{cv_prior}_BF_{BF_threshold}.txt",
        tips_only = "results/proximity/adaptive_branches_tips_only_CV_{cv_prior}_BF_{BF_threshold}.txt"
    shell:
        """
        # Extract internal nodes and rename them
        grep "^node_" {input.adaptive_branches} | sed 's/^node_/adaptive_node_/' > {output.renamed_nodes}

        # Extract tips only (no internal nodes)
        grep -v "^node_" {input.adaptive_branches} > {output.tips_only}
        """

rule convert_internal_nodes_to_vcf:
    """Convert renamed internal nodes from FASTA to VCF format for UShER placement"""
    input:
        reference = rules.params.input.usher_reference_fasta,
        internal_nodes = rules.rename_internal_nodes_fasta.output.renamed_fasta,
    output:
        vcf = "results/proximity/internal_nodes_renamed_CV_{cv_prior}_BF_{BF_threshold}.vcf"
    conda:
        "../envs/tree-proximity.yaml"
    shell:
        """
        # Concatenate reference + renamed internal nodes
        cat {input.reference} {input.internal_nodes} > {output.vcf}.temp.fasta

        # Convert to VCF
        faToVcf {output.vcf}.temp.fasta {output.vcf}

        # Clean up
        rm {output.vcf}.temp.fasta
        """

rule add_internal_nodes_to_tree:
    """Add renamed internal nodes to Viridian tree using UShER"""
    input:
        tree = rules.params.input.viridian_tree,
        vcf = rules.convert_internal_nodes_to_vcf.output.vcf,
    output:
        tree_updated = "results/proximity/viridian_tree_with_internal_nodes_CV_{cv_prior}_BF_{BF_threshold}.pb"
    threads: 4
    conda:
        "../envs/tree-proximity.yaml"
    shell:
        """
        usher \
            --load-mutation-annotated-tree {input.tree} \
            --vcf {input.vcf} \
            --save-mutation-annotated-tree {output.tree_updated} \
            --threads {threads}
        """

rule combine_all_adaptive_samples:
    """Combine tips + renamed internal nodes for extraction"""
    input:
        tips = rules.rename_internal_nodes_list.output.tips_only,
        renamed_nodes = rules.rename_internal_nodes_list.output.renamed_nodes,
    output:
        all_samples = "results/proximity/all_adaptive_samples_CV_{cv_prior}_BF_{BF_threshold}.txt"
    shell:
        """
        cat {input.tips} {input.renamed_nodes} > {output.all_samples}
        """

rule get_proximal_sequences_from_tree:
    """Find 5 nearest phylogenetic neighbors for ALL adaptive branches (tips + internal) using matUtils extract"""
    input:
        tree = rules.add_internal_nodes_to_tree.output.tree_updated,
        all_samples = rules.combine_all_adaptive_samples.output.all_samples,
    output:
        proximal_samples_all = "results/proximity/proximal_samples_all_CV_{cv_prior}_BF_{BF_threshold}.txt"
    conda:
        "../envs/tree-proximity.yaml"
    shell:
        """
        matUtils extract \
            --input-mat {input.tree} \
            --samples {input.all_samples} \
            --select-nearest 5 \
            --used-samples {output.proximal_samples_all}
        """


rule download_viridian_proximal_alignment:
    """Download Viridian sequences for the 5 nearest phylogenetic neighbors (includes adaptive branches)"""
    input:
        strain_names = rules.get_proximal_sequences_from_tree.output.proximal_samples_all,
    output:
        alignment_proximal = "results/proximity/alignment_proximal_CV_{cv_prior}_BF_{BF_threshold}.fasta",
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
            {input.strain_names} {output.alignment_proximal}
        """


rule combine_adaptive_with_proximal_context:
    """Combine adaptive branches with their 5 nearest phylogenetic neighbors"""
    input:
        adaptive = rules.get_alignment_all_adaptive_branches.output.alignment_all_adaptive_branches,
        proximal = rules.download_viridian_proximal_alignment.output.alignment_proximal,
    output:
        combined = "results/proximity/adaptive_with_phylo_context_CV_{cv_prior}_BF_{BF_threshold}.fasta"
    shell:
        """
        cat {input.adaptive} {input.proximal} > {output.combined}
        """

