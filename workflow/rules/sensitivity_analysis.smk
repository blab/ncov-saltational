rule subsample_tips_viridian_tree_per_year:
    """
    Subsample the Viridian tree to have the same number of tips per year.
    """
    input:
        metadata_dates = rules.clean_dates.output.metadata
    output:
        subsample_tips = "results/sensitivity_analysis/name_subsample_tips.txt"
    shell:
        """
        Rscript ./scripts/sensitivity-analysis/subsample_tips_same_every_year.R \
            --sra_metadata_dates {input.metadata_dates} \
            --seed 23047 \
            --name_subsamples {output.subsample_tips}
        """

rule downsample_viridian_tree_per_year:
    """
    Downsample the Viridian tree to have the same number of tips per year.
    """
    input:
        tree = rules.subset_viridian_tree.output.subset_tree,
        subsample_tips = rules.subsample_tips_viridian_tree_per_year.output.subsample_tips
    output:
        downsampled_tree = "results/sensitivity_analysis/viridian_tree_downsampled_per_year.pb.gz"
    shell:
        """
        matUtils extract \
            --input-mat {input.tree} \
            --samples {input.subsample_tips} \
            --write-mat {output.downsampled_tree}
        """

rule count_sites_per_branch_downsampled_tree:
    """Count the number of possible syn and nonsyn mutations available per gene on each branch of the MAT tree"""
    input:
        tree = rules.downsample_viridian_tree_per_year.output.downsampled_tree,
        reference_seq_gtf = rules.params.input.usher_reference_gtf,
        reference_seq_fasta = rules.params.input.usher_reference_fasta,
    output:
        df_sites = "results/sensitivity_analysis/df_sites_syn_non_syn_truncated_at_stop.tsv"
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

rule get_translations_per_node_downsampled_tree:
    """Get mutation information across tree"""
    input:
        mat = rules.downsample_viridian_tree_per_year.output.downsampled_tree,
        reference_seq_gtf = rules.params.input.usher_reference_gtf,
        reference_seq_fasta = rules.params.input.usher_reference_fasta,
    output:
        node_translations = "results/sensitivity_analysis/nodes_translations.tsv"
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

rule tidy_node_translations_downsampled_tree:
    """Tidy the translations from matUtils summary to a more easily workable format"""
    input:
        node_translations = rules.get_translations_per_node_downsampled_tree.output.node_translations
    output:
        aa_muts = temp("results/sensitivity_analysis/nodes_aa_muts_summary.tsv"),
        nt_muts = temp("results/sensitivity_analysis/nodes_nt_muts_summary.tsv")
    shell:
        """
        Rscript ./scripts/mat-manipulation/tidy_node_translations.R \
            --nodes_translations {input.node_translations} \
            --output_nt_mutations {output.nt_muts} \
            --output_aa_mutations {output.aa_muts}
        """

rule get_stop_codons_downsampled_tree:
    """Get stop codons across tree for each node"""
    input:
        tree = rules.downsample_viridian_tree_per_year.output.downsampled_tree,
        df_aa_muts = rules.tidy_node_translations_downsampled_tree.output.aa_muts,
    output:
        df_stop_codons = temp("results/sensitivity_analysis/df_stop_codons.tsv")
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/get_stop_codons_across_tree.py \
            --tree {input.tree} \
            --df_aa_muts {input.df_aa_muts} \
            --output_df_stop_codons {output.df_stop_codons}
        """

rule add_coded_status_aa_mut_downsampled_tree:
    """Add whether aa mutation is coded or not"""
    input:
        df_aa_muts = rules.tidy_node_translations_downsampled_tree.output.aa_muts,
        df_stop_codons = rules.get_stop_codons_downsampled_tree.output.df_stop_codons,
    output:
        df_aa_muts_coded = "results/sensitivity_analysis/nodes_aa_muts_summary_coded.tsv",
    shell:
        """
        Rscript ./scripts/mat-manipulation/get_non_coded_positions.R {input.df_aa_muts} {input.df_stop_codons} {output.df_aa_muts_coded}
        """

rule get_clades_per_node_downsampled_tree:
    """Get clade for each node"""
    input:
        tree = rules.downsample_viridian_tree_per_year.output.downsampled_tree,
    output:
        df_clades = "results/sensitivity_analysis/df_clade_per_node.tsv"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/get_clades_per_node.py \
            --tree {input.tree} \
            --output_df_clades {output.df_clades}
        """

rule get_parent_child_relationship_downsampled_tree:
    """Extract parent-child relationships from tree topology"""
    input:
        tree = rules.downsample_viridian_tree_per_year.output.downsampled_tree,
    output:
        df_parent_relationships = "results/sensitivity_analysis/df_parent_child_relationship.tsv"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/get_child_parent_relationship.py \
            --tree {input.tree} \
            --output_parent_relationship {output.df_parent_relationships}
        """

rule get_branches_leading_to_recombinant_clades_downsampled_tree:
    """Get branches leading to recombinant clades"""
    input:
        df_clades = rules.get_clades_per_node_downsampled_tree.output.df_clades,
        df_parent_relationships = rules.get_parent_child_relationship_downsampled_tree.output.df_parent_relationships
    output:
        name_strains_leading_to_recombinant_clades = "results/sensitivity_analysis/name_strains_leading_to_recombinant_clades.txt"
    shell:
        """
        Rscript scripts/mat-manipulation/get_branches_leading_to_recombinant_clades.R \
            --df_clades {input.df_clades} \
            --df_child_parent {input.df_parent_relationships} \
            --name_branches_leading_to_recombinant {output.name_strains_leading_to_recombinant_clades}
        """

rule prepare_data_model_downsampled_tree:
    """Prepare data for model run"""
    input:
        df_sites = rules.count_sites_per_branch_downsampled_tree.output.df_sites,
        df_aa_muts_coded = rules.add_coded_status_aa_mut_downsampled_tree.output.df_aa_muts_coded,
        name_recombinants = rules.get_branches_leading_to_recombinant_clades_downsampled_tree.output.name_strains_leading_to_recombinant_clades,
    output:
        data_viridian = "results/sensitivity_analysis/input_model/data_viridian_{gene}.tsv",
        data_viridian_grouped = "results/sensitivity_analysis/input_model/data_viridian_{gene}_grouped.tsv",
    shell:
        """
        Rscript ./scripts/bayesian-model/prepare_data_bayesian_model.R \
            --df_sites {input.df_sites} \
            --df_aa_muts_coded {input.df_aa_muts_coded} \
            --name_recombinants_to_remove {input.name_recombinants} \
            --output_data_viridian {output.data_viridian} \
            --output_data_viridian_grouped {output.data_viridian_grouped} \
            --genome_region {wildcards.gene}
        """

rule aggregate_genome_data_downsampled_tree:
    """Convert genome-wide gene-level data to aggregated format for downstream compatibility"""
    input:
        data_viridian_wide = "results/sensitivity_analysis/input_model/data_viridian_Genome.tsv"
    output:
        data_viridian_aggregated = "results/sensitivity_analysis/input_model/data_viridian_Genome_aggregated.tsv"
    shell:
        """
        Rscript ./scripts/bayesian-model/aggregate_genome_data.R \
            --data_viridian_wide {input.data_viridian_wide} \
            --output {output.data_viridian_aggregated}
        """

rule get_prior_parametrization_downsampled_tree:
    """Get prior parametrization for model run"""
    input:
        df_sites = rules.count_sites_per_branch_downsampled_tree.output.df_sites,
        df_aa_muts_coded = rules.add_coded_status_aa_mut_downsampled_tree.output.df_aa_muts_coded,
        df_parent_child_relationships = rules.get_parent_child_relationship_downsampled_tree.output.df_parent_relationships,
        name_recombinants = rules.get_branches_leading_to_recombinant_clades_downsampled_tree.output.name_strains_leading_to_recombinant_clades,
    output:
        df_prior_parametrization = "results/sensitivity_analysis/input_model/df_prior_parametrization.tsv"
    shell:
        """
        Rscript ./scripts/bayesian-model/get_prior_parametrization.R \
            --name_recombinants_to_remove {input.name_recombinants} \
            --df_sites {input.df_sites} \
            --df_aa_muts_coded {input.df_aa_muts_coded} \
            --df_parent_child {input.df_parent_child_relationships} \
            --df_prior_param {output.df_prior_parametrization}
        """

rule run_bayesian_model_downsampled_tree:
    """Run the Bayesian model for each gene and each CV prior"""
    input:
        data_viridian_grouped = rules.prepare_data_model_downsampled_tree.output.data_viridian_grouped,
        df_prior_parametrization = rules.get_prior_parametrization_downsampled_tree.output.df_prior_parametrization,
    output:
        df_param_estim = "results/sensitivity_analysis/output_model/df_param_estim_{gene}_CV_{cv_prior}.tsv"
    threads: 8
    params:
        script = lambda wildcards: (
            "./scripts/bayesian-model/run_bayesian_model_genome.R"
            if wildcards.gene == "Genome"
            else "./scripts/bayesian-model/run_bayesian_model_by_gene.R"
        ),
        n_cores_arg = lambda wildcards, threads: (
            f"--n_cores {threads}" if wildcards.gene == "Genome" else ""
        )
    shell:
        """
        Rscript {params.script} \
            --data_viridian_grouped {input.data_viridian_grouped} \
            --df_prior_param {input.df_prior_parametrization} \
            --CV_prior {wildcards.cv_prior} \
            --output_param_estim {output.df_param_estim} \
            --genome_region {wildcards.gene} \
            {params.n_cores_arg}
        """

rule get_BF_across_model_runs_downsampled_tree:
    """Get Bayes Factors across model runs for each gene and each CV prior"""
    input:
        df_param_estim = rules.run_bayesian_model_downsampled_tree.output.df_param_estim,
    output:
        df_BF = "results/sensitivity_analysis/output_model/df_BF_{gene}_CV_{cv_prior}.tsv"
    shell:
        """
        Rscript ./scripts/bayesian-model/compute_BF_across_branches.R \
            --df_param_estim {input.df_param_estim} \
            --df_BF {output.df_BF}
        """

rule get_name_adaptive_branches_downsampled_tree:
    """Get the names of the adaptive paths across the tree for each gene and each CV prior"""
    input:
        data_viridian = "results/sensitivity_analysis/input_model/data_viridian_{gene}.tsv",
        data_viridian_grouped = rules.prepare_data_model_downsampled_tree.output.data_viridian_grouped,
        df_BF = rules.get_BF_across_model_runs_downsampled_tree.output.df_BF
    output:
        name_adaptive_branches = "results/sensitivity_analysis/adaptive_branches/name_adaptive_branches_{gene}_CV_{cv_prior}_BF_{BF_threshold}.tsv"
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

rule concat_name_adaptive_branches_downsampled_tree:
    """Concatenate adaptive branch names across all genes"""
    input:
        name_adaptive_branches_across_genes = lambda wildcards: expand("results/sensitivity_analysis/adaptive_branches/name_adaptive_branches_{gene}_CV_{cv_prior}_BF_{BF_threshold}.tsv",
                   gene= GENES, cv_prior=wildcards.cv_prior, BF_threshold=wildcards.BF_threshold),
    output:
        concat_name_adaptive_branches = "results/sensitivity_analysis/adaptive_branches/concat_name_adaptive_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        tsv-append -H {input.name_adaptive_branches_across_genes} | tsv-uniq >> {output.concat_name_adaptive_branches}
        """

rule name_saltational_branches_downsampled_tree:
    """Get the names of the saltational paths across the tree for each gene and each CV prior"""
    input:
        name_adaptive_branches = rules.concat_name_adaptive_branches_downsampled_tree.output.concat_name_adaptive_branches,
        df_aa_muts = rules.add_coded_status_aa_mut_downsampled_tree.output.df_aa_muts_coded
    output:
        name_saltational_branches = "results/sensitivity_analysis/saltational_branches/name_saltational_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv",
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

rule convert_tree_newick_downsampled_tree:
    """Convert UShER tree to newick format - This will be used to run Chronumental"""
    input:
        usher_tree = rules.downsample_viridian_tree_per_year.output.downsampled_tree
    output:
        tree_newick = "results/sensitivity_analysis/viridian_tree.nwk"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        matUtils extract \
            --input-mat {input.usher_tree} \
            --write-tree {output.tree_newick}
        """

rule chronumental_downsampled_tree:
    """Run chronumental to estimate time tree"""
    input:
        tree = rules.convert_tree_newick_downsampled_tree.output.tree_newick,
        metadata = rules.clean_dates.output.metadata
    output:
        inferred_dates = "results/sensitivity_analysis/chronumental/chronumental_inferred_dates.tsv",
        timetree = "results/sensitivity_analysis/chronumental/chronumental_time_tree.pb.gz"
    params:
        reference_node = "ERR4424384"
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