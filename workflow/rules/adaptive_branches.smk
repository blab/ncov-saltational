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
"""


rule get_name_adaptive_branches:
    """Get the names of the adaptive paths across the tree for each gene and each CV prior"""
    input:
        data_viridian = "results/input_model/data_viridian_{gene}.tsv",
        data_viridian_grouped = rules.prepare_data_model.output.data_viridian_grouped,
        df_BF = rules.get_BF_across_model_runs.output.df_BF
    output:
        name_adaptive_branches = "results/adaptive_branches/name_adaptive_branches_{gene}_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    params:
        join_keys = lambda wildcards: (
            "syn_sites_E,syn_sites_M,syn_sites_N,syn_sites_ORF1a,syn_sites_ORF1b,syn_sites_ORF3a,syn_sites_ORF6,syn_sites_ORF7a,syn_sites_ORF7b,syn_sites_ORF8,syn_sites_ORF9b,syn_sites_S,nonsyn_sites_E,nonsyn_sites_M,nonsyn_sites_N,nonsyn_sites_ORF1a,nonsyn_sites_ORF1b,nonsyn_sites_ORF3a,nonsyn_sites_ORF6,nonsyn_sites_ORF7a,nonsyn_sites_ORF7b,nonsyn_sites_ORF8,nonsyn_sites_ORF9b,nonsyn_sites_S,synonymous_muts_E,synonymous_muts_M,synonymous_muts_N,synonymous_muts_ORF1a,synonymous_muts_ORF1b,synonymous_muts_ORF3a,synonymous_muts_ORF6,synonymous_muts_ORF7a,synonymous_muts_ORF7b,synonymous_muts_ORF8,synonymous_muts_ORF9b,synonymous_muts_S,non_synonymous_muts_E,non_synonymous_muts_M,non_synonymous_muts_N,non_synonymous_muts_ORF1a,non_synonymous_muts_ORF1b,non_synonymous_muts_ORF3a,non_synonymous_muts_ORF6,non_synonymous_muts_ORF7a,non_synonymous_muts_ORF7b,non_synonymous_muts_ORF8,non_synonymous_muts_ORF9b,non_synonymous_muts_S"
            if wildcards.gene == "Genome"
            else "synonymous_muts,non_synonymous_muts,syn_sites,nonsyn_sites"
        )
    shell:
        """
        tsv-join -H \
            --filter-file <(
                tsv-filter -H --gt BF_omega_greater_1:{wildcards.BF_threshold} {input.df_BF} |
                tsv-select -H -f i_row |
                tsv-join -H --filter-file {input.data_viridian_grouped} --key-fields i_row --append-fields {params.join_keys}
            ) \
            --key-fields {params.join_keys} {input.data_viridian} |
            tsv-select -H -f strain |
            sed 1d > {output.name_adaptive_branches}
        """


rule concat_name_adaptive_branches:
    """Concatenate adaptive branch names across all genes"""
    input:
        name_adaptive_branches_across_genes = lambda wildcards: expand("results/adaptive_branches/name_adaptive_branches_{gene}_CV_{cv_prior}_BF_{BF_threshold}.tsv",
                   gene= GENES, cv_prior=wildcards.cv_prior, BF_threshold=wildcards.BF_threshold),
    output:
        concat_name_adaptive_branches = "results/adaptive_branches/concat_name_adaptive_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        tsv-append -H {input.name_adaptive_branches_across_genes} | tsv-uniq >> {output.concat_name_adaptive_branches}
        """

rule name_saltational_branches:
    """Get the names of the saltational paths across the tree for each gene and each CV prior"""
    input:
        name_adaptive_branches = rules.concat_name_adaptive_branches.output.concat_name_adaptive_branches,
        df_aa_muts = rules.add_coded_status_aa_mut.output.df_aa_muts_coded
    output:
        name_saltational_branches = "results/saltational_branches/name_saltational_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv",
    params:
        percentile_threshold_nonsyn = 0.995
    shell:
        """
        Rscript ./scripts/branch-analysis/get_saltational_branches.R \
            --name_adaptive_branches {input.name_adaptive_branches} \
            --df_aa_muts_coded {input.df_aa_muts} \
            --percentile_threshold_nonsyn {params.percentile_threshold_nonsyn} \
            --output_name_saltational_branches {output.name_saltational_branches}
        """

rule name_saltational_branches_by_genome_region:
    """Get the names of the saltational paths across the tree for each gene and each CV prior"""
    input:
        name_adaptive_branches = rules.get_name_adaptive_branches.output.name_adaptive_branches,
        df_aa_muts = rules.add_coded_status_aa_mut.output.df_aa_muts_coded
    output:
        name_saltational_branches = "results/saltational_branches/name_saltational_branches_{gene}_CV_{cv_prior}_BF_{BF_threshold}.tsv",
    params:
        percentile_threshold_nonsyn = 0.995
    shell:
        """
        Rscript ./scripts/branch-analysis/get_saltational_branches.R \
            --name_adaptive_branches {input.name_adaptive_branches} \
            --df_aa_muts_coded {input.df_aa_muts} \
            --percentile_threshold_nonsyn {params.percentile_threshold_nonsyn} \
            --output_name_saltational_branches {output.name_saltational_branches}
        """

rule get_name_internal_nodes_saltational_branches:
    """Identify which saltational branches are internal nodes (vs tips)"""
    input:
        name_saltational_branches = rules.name_saltational_branches.output.name_saltational_branches,
    output:
        name_internal_nodes_saltational = "results/saltational_branches/name_saltational_internal_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        tsv-filter -H --regex 'strain:^node_' <(printf "strain\n"; cat {input.name_saltational_branches}) | sed 1d > {output.name_internal_nodes_saltational}
        """

rule ancestral_sequence_reconstruction_saltational:
    """Perform ancestral sequence reconstruction for internal saltational branches"""
    input:
        df_child_parent_relationships = rules.get_parent_child_relationship.output.df_parent_relationships,
        name_internal_nodes_saltational = rules.get_name_internal_nodes_saltational_branches.output.name_internal_nodes_saltational,
        df_nuc_mut = rules.char_nuc_across_tree.output.df_nuc_mut,
        ref_fasta = rules.params.input.usher_reference_fasta,
    output:
        output_fasta = "results/saltational_branches/sequences_saltational_internal_CV_{cv_prior}_BF_{BF_threshold}.fasta"
    threads: 4
    shell:
        """
        Rscript ./scripts/mat-manipulation/reconstruct_ancestral_seq.R \
            --df_nuc_mut {input.df_nuc_mut} \
            --ref_fasta {input.ref_fasta} \
            --nodes_to_reconstruct {input.name_internal_nodes_saltational} \
            --df_child_parent {input.df_child_parent_relationships} \
            --output_ancestral_fasta {output.output_fasta} \
            --n_cores {threads}
        """

rule rename_internal_nodes_saltational_fasta:
    """Rename internal nodes in FASTA to avoid conflicts when adding to tree (node_XXX -> adaptive_node_XXX)"""
    input:
        internal_nodes = rules.ancestral_sequence_reconstruction_saltational.output.output_fasta,
    output:
        renamed_fasta = temp("results/saltational_branches/sequences_saltational_internal_renamed_CV_{cv_prior}_BF_{BF_threshold}.fasta")
    shell:
        """
        # Rename sequence headers: node_XXX -> adaptive_node_XXX
        sed 's/>node_/>adaptive_node_/g' {input.internal_nodes} > {output.renamed_fasta}
        """

rule convert_internal_nodes_saltational_to_vcf:
    """Convert renamed internal nodes from FASTA to VCF format for UShER placement"""
    input:
        reference = rules.params.input.usher_reference_fasta,
        internal_nodes = rules.rename_internal_nodes_saltational_fasta.output.renamed_fasta,
    output:
        vcf = temp("results/saltational_branches/sequences_saltational_internal_renamed_CV_{cv_prior}_BF_{BF_threshold}.vcf")
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

rule add_internal_nodes_saltational_to_tree:
    """Add renamed internal nodes to Viridian tree using UShER - this tree will then be used for proximity-based subsampling of background sequences for adaptive branches"""
    input:
        tree = rules.subset_viridian_tree.output.subset_tree,
        vcf = rules.convert_internal_nodes_saltational_to_vcf.output.vcf,
    output:
        tree_updated = "results/saltational_branches/viridian_tree_with_internal_nodes_CV_{cv_prior}_BF_{BF_threshold}.pb.gz"
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

rule get_saltational_nodes_list_for_proximity:
    """Create renamed list of internal saltational nodes and combine with tips"""
    input:
        saltational_branches = rules.name_saltational_branches.output.name_saltational_branches,
    output:
        renamed_nodes = "results/saltational_branches/internal_nodes_renamed_list_CV_{cv_prior}_BF_{BF_threshold}.txt",
        tips_only = "results/saltational_branches/saltational_branches_tips_only_CV_{cv_prior}_BF_{BF_threshold}.txt",
        combined_saltational_branches = "results/saltational_branches/saltational_branches_combined_CV_{cv_prior}_BF_{BF_threshold}.txt"
    shell:
        """
        # Extract internal nodes and rename them
        grep "^node_" {input.saltational_branches} | sed 's/^node_/adaptive_node_/' > {output.renamed_nodes}

        # Extract tips only (no internal nodes)
        grep -v "^node_" {input.saltational_branches} > {output.tips_only}

        # Combine renamed internal nodes and tips
        cat {output.renamed_nodes} {output.tips_only} > {output.combined_saltational_branches}
        """

rule sample_proximity_to_saltational_from_tree:
    """Sample background sequences in proximity to saltational branches from updated tree with internal nodes"""
    input:
        tree = rules.add_internal_nodes_saltational_to_tree.output.tree_updated,
        name_saltational_branches = rules.get_saltational_nodes_list_for_proximity.output.combined_saltational_branches,
    output:
        mat = "results/tree_with_saltational/proximity_tree_CV_{cv_prior}_BF_{BF_threshold}.pb.gz",
        tree = "results/tree_with_saltational/proximity_tree_CV_{cv_prior}_BF_{BF_threshold}.nwk",
        taxodium = "results/tree_with_saltational/proximity_tree_CV_{cv_prior}_BF_{BF_threshold}.jsonl.gz",
        proximity_samples = "results/tree_with_saltational/proximity_samples_CV_{cv_prior}_BF_{BF_threshold}.txt"
    params:
        proximity_param = 10,
        n_random_samples = 3000,
        reference_gtf = rules.params.input.usher_reference_gtf,
        reference_fasta = rules.params.input.usher_reference_fasta
    conda:
        "../envs/tree-proximity.yaml"
    shell:
        """
        matUtils extract \
            --input-mat {input.tree} \
            --samples {input.name_saltational_branches} \
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

rule normalize_saltational_tree_branch_lengths:
    """Normalize branch lengths from mutation counts to substitutions per site for TreeTime compatibility"""
    input:
        tree = rules.sample_proximity_to_saltational_from_tree.output.mat,
    output:
        tree_normalized = "results/tree_with_saltational/proximity_tree_normalized_CV_{cv_prior}_BF_{BF_threshold}.nwk"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python scripts/tree-proximity/convert_tree_units_to_mutations_per_site.py \
            --tree {input.tree} \
            --genome_length 29903 \
            --output_newick {output.tree_normalized}
        """