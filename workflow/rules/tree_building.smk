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
        name_internal_nodes_adaptive = "results/tree_building/name_adaptive_internal_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv"
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
        output_fasta = "results/tree_building/alignment_adaptive_internal_nodes_CV_{cv_prior}_BF_{BF_threshold}.fasta"
    threads: 4
    shell:
        """
        Rscript ./scripts/mat-manipulation/reconstruct_ancestral_seq.R \
            --df_nuc_mut {input.df_nuc_mut} \
            --ref_fasta {input.ref_fasta} \
            --nodes_to_reconstruct {input.name_internal_nodes_adaptive} \
            --df_child_parent {input.df_child_parent_relationships} \
            --output_ancestral_fasta {output.output_fasta} \
            --n_cores {threads}
        """

rule get_viridian_alignment_adaptive_branches:
    """Download sequences for adaptive branch tips from Viridian"""
    input:
        name_adaptive_branches = rules.concat_name_adaptive_branches.output.concat_name_adaptive_branches,
    output:
        alignment_adaptive_final_nodes = "results/tree_building/alignment_adaptive_final_nodes_CV_{cv_prior}_BF_{BF_threshold}.fasta",
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

rule rename_internal_nodes_fasta:
    """Rename internal nodes in FASTA to avoid conflicts when adding to tree (node_XXX -> adaptive_node_XXX)"""
    input:
        internal_nodes = rules.ancestral_sequence_reconstruction.output.output_fasta,
    output:
        renamed_fasta = temp("results/tree_building/internal_nodes_renamed_CV_{cv_prior}_BF_{BF_threshold}.fasta")
    shell:
        """
        # Rename sequence headers: node_XXX -> adaptive_node_XXX
        sed 's/>node_/>adaptive_node_/g' {input.internal_nodes} > {output.renamed_fasta}
        """


rule get_alignment_all_adaptive_branches:
    """Combine internal node and tip sequences for all adaptive branches"""
    input:
        alignment_adaptive_final_nodes = rules.get_viridian_alignment_adaptive_branches.output.alignment_adaptive_final_nodes,
        alignment_adaptive_internal_nodes = rules.rename_internal_nodes_fasta.output.renamed_fasta,
    output:
        alignment_all_adaptive_branches = "results/alignment_adaptive_branches/alignment_adaptive_branches_CV_{cv_prior}_BF_{BF_threshold}.fasta"
    shell:
        """
        cat {input.alignment_adaptive_final_nodes} {input.alignment_adaptive_internal_nodes} > {output.alignment_all_adaptive_branches}
        """

rule convert_internal_nodes_to_vcf:
    """Convert renamed internal nodes from FASTA to VCF format for UShER placement"""
    input:
        reference = rules.params.input.usher_reference_fasta,
        internal_nodes = rules.rename_internal_nodes_fasta.output.renamed_fasta,
    output:
        vcf = temp("results/tree_building/internal_nodes_renamed_CV_{cv_prior}_BF_{BF_threshold}.vcf")
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
    """Add renamed internal nodes to Viridian tree using UShER - this tree will then be used for proximity-based subsampling of background sequences for adaptive branches"""
    input:
        tree = rules.subset_viridian_tree.output.subset_tree,
        vcf = rules.convert_internal_nodes_to_vcf.output.vcf,
    output:
        tree_updated = "results/tree_building/viridian_tree_with_internal_nodes_CV_{cv_prior}_BF_{BF_threshold}.pb"
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

rule get_adaptive_nodes_list_for_proximity:
    """Create renamed list of internal adaptive nodes and combine with tips"""
    input:
        adaptive_branches = rules.concat_name_adaptive_branches.output.concat_name_adaptive_branches,
    output:
        renamed_nodes = "results/tree_building/internal_nodes_renamed_list_CV_{cv_prior}_BF_{BF_threshold}.txt",
        tips_only = "results/tree_building/adaptive_branches_tips_only_CV_{cv_prior}_BF_{BF_threshold}.txt",
        combined_adaptive_branches = "results/tree_building/adaptive_branches_combined_CV_{cv_prior}_BF_{BF_threshold}.txt"
    shell:
        """
        # Extract internal nodes and rename them
        grep "^node_" {input.adaptive_branches} | sed 's/^node_/adaptive_node_/' > {output.renamed_nodes}

        # Extract tips only (no internal nodes)
        grep -v "^node_" {input.adaptive_branches} > {output.tips_only}

        # Combine renamed internal nodes and tips
        cat {output.renamed_nodes} {output.tips_only} > {output.combined_adaptive_branches}
        """


rule sample_proximity_from_tree:
    """Sample background sequences in proximity to adaptive branches from updated tree with internal nodes"""
    input:
        tree = rules.add_internal_nodes_to_tree.output.tree_updated,
        name_adaptive_branches = rules.get_adaptive_nodes_list_for_proximity.output.combined_adaptive_branches,
    output:
        mat = "results/tree_building/proximity_tree_CV_{cv_prior}_BF_{BF_threshold}.pb.gz",
        tree = "results/tree_building/proximity_tree_CV_{cv_prior}_BF_{BF_threshold}.nwk",
        taxodium = "results/tree_building/proximity_tree_CV_{cv_prior}_BF_{BF_threshold}.jsonl.gz",
        proximity_samples = "results/tree_building/proximity_samples_CV_{cv_prior}_BF_{BF_threshold}.txt"
    params:
        proximity_param = 5,
        n_random_samples = 6000,
        reference_gtf = rules.params.input.usher_reference_gtf,
        reference_fasta = rules.params.input.usher_reference_fasta
    conda:
        "../envs/tree-proximity.yaml"
    shell:
        """
        matUtils extract \
            --input-mat {input.tree} \
            --samples {input.name_adaptive_branches} \
            --select-nearest {params.proximity_param} \
            --add-random {params.n_random_samples} \
            --used-samples {output.proximity_samples} \
            --write-mat {output.mat} \
            --write-tree {output.tree} \
            --write-taxodium {output.taxodium} \
            --title "Phylogenetic tree with adaptive branches" \
            --include-nt \
            --input-gtf {params.reference_gtf} \
            --input-fasta {params.reference_fasta}

        """

rule normalize_tree_branch_lengths:
    """Normalize branch lengths from mutation counts to substitutions per site for TreeTime compatibility"""
    input:
        tree = rules.sample_proximity_from_tree.output.mat,
    output:
        tree_normalized = "results/tree_building/proximity_tree_normalized_CV_{cv_prior}_BF_{BF_threshold}.nwk"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python scripts/tree-proximity/convert_tree_units_to_mutations_per_site.py \
            --tree {input.tree} \
            --genome_length 29903 \
            --output_newick {output.tree_normalized}
        """

rule get_proximity_background_samples_only:
    """Filter proximity samples to exclude ALL adaptive branches (both internal and terminal) - keep only background/context samples"""
    input:
        proximity_samples = rules.sample_proximity_from_tree.output.proximity_samples,
        adaptive_branches = rules.get_adaptive_nodes_list_for_proximity.output.combined_adaptive_branches,
    output:
        background_samples = "results/tree_building/proximity_background_only_CV_{cv_prior}_BF_{BF_threshold}.txt"
    shell:
        """
        # Remove adaptive_node_* internal nodes, then remove all adaptive branches
        comm -23 \
            <(grep -v "^adaptive_node_" {input.proximity_samples} | sort) \
            <(cat {input.adaptive_branches} | sort) \
            > {output.background_samples}
        """

rule download_viridian_sequences_proximity_background:
    """Download sequences for background/context samples from Viridian (excludes adaptive branches)"""
    input:
        background_samples = rules.get_proximity_background_samples_only.output.background_samples,
    output:
        alignment_background = "results/tree_building/alignment_proximity_background_CV_{cv_prior}_BF_{BF_threshold}.fasta"
    params:
        url_viridian_1 = "https://api.figshare.com/v2/file/download/45969777",
        url_viridian_2 = "https://api.figshare.com/v2/file/download/49692480",
        file_path_viridian_1 = "data/tmp/viridian_alignment_1_background_CV_{cv_prior}_BF_{BF_threshold}.tar",
        file_path_viridian_2 = "data/tmp/viridian_alignment_2_background_CV_{cv_prior}_BF_{BF_threshold}.tar",
        dir_unpacked_viridian = "data/tmp/viridian_background_CV_{cv_prior}_BF_{BF_threshold}/",
    shell:
        """
        Rscript ./scripts/mat-manipulation/download_viridian_alignment.R \
            {params.url_viridian_1} {params.url_viridian_2} \
            {params.file_path_viridian_1} {params.file_path_viridian_2} \
            {params.dir_unpacked_viridian} {params.dir_unpacked_viridian} \
            {input.background_samples} {output.alignment_background}
        """

rule combine_proximity_alignment:
    """Combine adaptive branches alignment with background sequences for full proximity tree alignment"""
    input:
        adaptive_alignment = rules.get_alignment_all_adaptive_branches.output.alignment_all_adaptive_branches,
        background_alignment = rules.download_viridian_sequences_proximity_background.output.alignment_background,
    output:
        combined_alignment = "results/tree_building/alignment_proximity_all_CV_{cv_prior}_BF_{BF_threshold}.fasta"
    shell:
        """
        cat {input.adaptive_alignment} {input.background_alignment} > {output.combined_alignment}
        """

rule augur_align_proximity_sequences:
    """Align proximity sequences (adaptive + background) to reference genome"""
    input:
        sequences = rules.combine_proximity_alignment.output.combined_alignment,
        reference = rules.params.input.usher_reference_fasta,
    output:
        alignment = "results/tree_building/aligned_proximity_sequences_CV_{cv_prior}_BF_{BF_threshold}.fasta"
    threads: 4
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --nthreads {threads}
        """

rule get_metadata_proximity_samples:
    """Aggregate metadata (dates) for all proximity samples (adaptive + context + random)"""
    input:
        proximity_samples = rules.sample_proximity_from_tree.output.proximity_samples,
        chronumental = rules.chronumental.output.inferred_dates,
        metadata = rules.tidy_final_metadata.output.subset_metadata,
    output:
        metadata = "results/tree_building/metadata_proximity_samples_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        Rscript ./scripts/tree-proximity/get_metadata_proximity_samples.R \
            --proximity_samples {input.proximity_samples} \
            --chronumental {input.chronumental} \
            --metadata {input.metadata} \
            --output_metadata {output.metadata}
        """

rule proximity_time_tree:
    """Refine proximity tree with temporal information using augur refine"""
    input:
        tree = rules.normalize_tree_branch_lengths.output.tree_normalized,
        alignment = rules.augur_align_proximity_sequences.output.alignment,
        metadata = rules.get_metadata_proximity_samples.output.metadata,
    output:
        tree = "results/tree_building/proximity_time_tree_CV_{cv_prior}_BF_{BF_threshold}.nwk",
        node_data = "results/tree_building/proximity_branch_lengths_CV_{cv_prior}_BF_{BF_threshold}.json"
    threads: 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent opt \
            --divergence-units mutations \
            --keep-root \
            --year-bounds 2019 2024
        """

rule export_time_tree:
    """Export proximity time tree as Auspice JSON for visualization"""
    input:
        tree = rules.proximity_time_tree.output.tree,
        node_data = rules.proximity_time_tree.output.node_data,
    output:
        auspice_json = "results/tree_building/auspice_tree_CV_{cv_prior}_BF_{BF_threshold}.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --node-data {input.node_data} \
            --output {output.auspice_json}
        """