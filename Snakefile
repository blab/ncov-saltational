GENES = ['ORF1a', 'ORF1b', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF9b', 'Genome']
CV_PRIORS = [0.95, 1.5, 2.0]
BF_THRESHOLDS = [31.62277660]

rule all:
    input:
        df_nuc_mut = "results/mutations/df_nuc_mut.tsv",
        df_aa_muts_coded = "results/mutations/nodes_aa_muts_summary_coded.tsv",
        df_sites = "results/bte/df_sites_syn_non_syn_truncated_at_stop.tsv",

        df_parent_relationships = "results/bte/df_parent_child_relationship.tsv",
        df_n_descendants = "results/bte/df_n_descendants_per_node.tsv",

        data_viridian = expand("results/input_model/data_viridian_{gene}.tsv", gene=GENES),
        data_viridian_grouped = expand("results/input_model/data_viridian_{gene}_grouped.tsv", gene=GENES),
        df_prior_parametrization = "results/input_model/df_prior_parametrization.tsv",
        #posterior_simul = expand("results/posterior_predictive_checks/posterior_predictive_check_{gene}_CV_{cv_prior}.tsv", gene=GENES, cv_prior = CV_PRIORS),

        df_param_estim = expand("results/output_model/df_param_estim_{gene}_CV_{cv_prior}.tsv", gene=GENES, cv_prior = CV_PRIORS),
        df_BF = expand("results/output_model/df_BF_{gene}_CV_{cv_prior}.tsv", gene=GENES, cv_prior = CV_PRIORS),
        name_adaptive_branches = expand("results/adaptive_branches/name_adaptive_branches_{gene}_CV_{cv_prior}_BF_{BF_threshold}.tsv", gene=GENES, cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),
        concat_name_adaptive_branches = expand("results/adaptive_branches/concat_name_adaptive_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),

        ancestral_fasta = expand("results/alignment_adaptive_branches/alignment_adaptive_internal_nodes_CV_{cv_prior}_BF_{BF_threshold}.fasta", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),
        alignment_adaptive_final_nodes = expand("results/alignment_adaptive_branches/alignment_adaptive_final_nodes_CV_{cv_prior}_BF_{BF_threshold}.fasta", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),
        alignment_all_adaptive_branches = expand("results/alignment_adaptive_branches/alignment_adaptive_branches_CV_{cv_prior}_BF_{BF_threshold}.fasta", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),

rule params:
    input:
        usher_reference_gtf = "params/reference/reference_seq_all.gtf",
        usher_reference_fasta = "params/reference/reference_seq.fasta",
        viridian_tree = "data/viridian/tree.all_viridian.202409.pb.gz",

rule char_nuc_across_tree:
    "List of nucleotide mutations across the UShER tree"
    input:
        tree = rules.params.input.viridian_tree,
    output:
        df_nuc_mut = "results/mutations/df_nuc_mut.tsv"
    conda:
        "envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/get_nuc_mutation_across_genome.py \
            --tree {input.tree} \
            --output_df_nuc_mutations {output.df_nuc_mut}
        """

rule get_translations_per_node:
    "Get mutation information across tree"
    input:
        mat = rules.params.input.viridian_tree,
        reference_seq_gtf = rules.params.input.usher_reference_gtf,
        reference_seq_fasta = rules.params.input.usher_reference_fasta,
    output:
        node_translations = "results/mutations/nodes_translations.tsv"
    shell:
        """
        matUtils summary \
            --input-mat {input.mat} \
            --input-gtf {input.reference_seq_gtf} \
            --input-fasta {input.reference_seq_fasta} \
            --translate {output.node_translations}
        """

rule tidy_node_translations:
    "Tidy the translations from matUtils summary to a more easily workable format"
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
    "Get stop codons across tree for each node"
    input:
        tree = rules.params.input.viridian_tree,
        df_aa_muts = rules.tidy_node_translations.output.aa_muts,
    output:
        df_stop_codons = "results/bte/df_stop_codons.tsv"
    conda:
        "envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/get_stop_codons_across_tree.py \
            --tree {input.tree} \
            --df_aa_muts {input.df_aa_muts} \
            --output_df_stop_codons {output.df_stop_codons}
        """

rule add_coded_status_aa_mut:
    "Add whether aa mutation is coded or not"
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
    input:
        tree = rules.params.input.viridian_tree,
    output:
        df_parent_relationships = "results/bte/df_parent_child_relationship.tsv"
    conda:
        "envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/get_child_parent_relationship.py \
            --tree {input.tree} \
            --output_parent_relationship {output.df_parent_relationships}
        """

rule get_n_descendants_per_node:
    "Get number of descendants per node"
    input:
        tree = rules.params.input.viridian_tree,
    output:
        df_n_descendants = "results/bte/df_n_descendants_per_node.tsv",
        df_n_children = "results/bte/df_n_children_per_node.tsv"
    conda:
        "envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/get_n_descendants_per_node.py \
            --tree {input.tree} \
            --output_df_n_descendants {output.df_n_descendants} \
            --output_df_n_children {output.df_n_children}
        """

rule count_sites_per_branch:
    "Count the number of possible syn and nonsyn mutations available per gene on each branch of the MAT tree"
    input:
        tree = rules.params.input.viridian_tree,
        reference_seq_gtf = rules.params.input.usher_reference_gtf,
        reference_seq_fasta = rules.params.input.usher_reference_fasta,
    output:
        df_sites = "results/bte/df_sites_syn_non_syn_truncated_at_stop.tsv"
    conda:
        "envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/count_sites_per_branch.py \
            --ref-fasta {input.reference_seq_fasta} \
            --ref-gtf {input.reference_seq_gtf} \
            --tree-path {input.tree} \
            --results-file {output.df_sites} \
            --truncate-at-stop
        """

rule prepare_data_model:
    "Prepare data for model run"
    input:
        df_sites = rules.count_sites_per_branch.output.df_sites,
        df_aa_muts_coded = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
    output:
        data_viridian = "results/input_model/data_viridian_{gene}.tsv", 
        data_viridian_grouped = "results/input_model/data_viridian_{gene}_grouped.tsv",
    shell:
        """
        Rscript ./scripts/bayesian-model/prepare_data_bayesian_model.R \
            --df_sites {input.df_sites} \
            --df_aa_muts_coded {input.df_aa_muts_coded} \
            --output_data_viridian {output.data_viridian} \
            --output_data_viridian_grouped {output.data_viridian_grouped} \
            --genome_region {wildcards.gene}
        """

rule get_prior_parametrization:
    "Get prior parametrization for model run"
    input:
        df_sites = rules.count_sites_per_branch.output.df_sites,
        df_aa_muts_coded = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
        df_parent_child_relationships = rules.get_parent_child_relationship.output.df_parent_relationships,
    output:
        df_prior_parametrization = "results/input_model/df_prior_parametrization.tsv"
    shell:
        """
        Rscript ./scripts/bayesian-model/get_prior_parametrization.R \
            --df_sites {input.df_sites} \
            --df_aa_muts_coded {input.df_aa_muts_coded} \
            --df_parent_child {input.df_parent_child_relationships} \
            --df_prior_param {output.df_prior_parametrization}
        """

rule run_posterior_predictive_checks:
    "Run posterior predictive checks for each gene and each CV prior"
    input:
        data_viridian_grouped = rules.prepare_data_model.output.data_viridian_grouped,
        df_prior_parametrization = rules.get_prior_parametrization.output.df_prior_parametrization,
    output:
        posterior_simul = "results/posterior_predictive_checks/posterior_predictive_check_{gene}_CV_{cv_prior}.tsv"
    shell:
        """
        Rscript ./scripts/bayesian-model/run_posterior_predictive_checks.R \
            --df_prior_param {input.df_prior_parametrization} \
            --genome_region {wildcards.gene} \
            --CV_prior {wildcards.cv_prior} \
            --data_viridian_grouped {input.data_viridian_grouped} \
            --output_posterior_simul {output.posterior_simul}
        """

rule run_bayesian_model:
    "Run the Bayesian model for each gene and each CV prior"
    input:
        data_viridian_grouped = rules.prepare_data_model.output.data_viridian_grouped,
        df_prior_parametrization = rules.get_prior_parametrization.output.df_prior_parametrization,
    output:
        df_param_estim = "results/output_model/df_param_estim_{gene}_CV_{cv_prior}.tsv"
    shell:
        """
        Rscript ./scripts/bayesian-model/run_bayesian_model_by_gene.R \
            --data_viridian_grouped {input.data_viridian_grouped} \
            --df_prior_param {input.df_prior_parametrization} \
            --CV_prior {wildcards.cv_prior} \
            --output_param_estim {output.df_param_estim} \
            --genome_region {wildcards.gene}
        """

rule get_BF_across_model_runs:
    "Get Bayes Factors across model runs for each gene and each CV prior"
    input:
        df_param_estim = rules.run_bayesian_model.output.df_param_estim,
    output:
        df_BF = "results/output_model/df_BF_{gene}_CV_{cv_prior}.tsv"
    shell:
        """
        Rscript ./scripts/bayesian-model/compute_BF_across_branches.R \
            --df_param_estim {input.df_param_estim} \
            --df_BF {output.df_BF}
        """

rule get_name_adaptive_branches:
    "Get the names of the adaptive paths across the tree for each gene and each CV prior"
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
    input:
        name_adaptive_branches_across_genes = lambda wildcards: expand("results/adaptive_branches/name_adaptive_branches_{gene}_CV_{curr_cv_prior}_BF_{curr_BF_threshold}.tsv", gene=GENES, curr_cv_prior = wildcards.cv_prior, curr_BF_threshold = wildcards.BF_threshold),
    output:
        concat_name_adaptive_branches = "results/adaptive_branches/concat_name_adaptive_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        tsv-append -H {input.name_adaptive_branches_across_genes} | tsv-uniq >> {output.concat_name_adaptive_branches}
        """


rule get_name_internal_nodes_adaptive:
    input:
        name_adaptive_branches = rules.concat_name_adaptive_branches.output.concat_name_adaptive_branches,
    output:
        name_internal_nodes_adaptive = "results/adaptive_branches/name_adaptive_internal_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        tsv-filter -H --regex 'strain:^node_' <(printf "strain\n"; cat {input.name_adaptive_branches}) | sed 1d > {output.name_internal_nodes_adaptive}
        """

rule ancestral_sequence_reconstruction:
    "Perform ancestral sequence reconstruction for internal adaptive branches"
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
    input:
        alignment_adaptive_final_nodes = rules.get_viridian_alignment_adaptive_branches.output.alignment_adaptive_final_nodes,
        alignment_adaptive_internal_nodes = rules.ancestral_sequence_reconstruction.output.output_fasta,
    output:
        alignment_all_adaptive_branches = "results/alignment_adaptive_branches/alignment_adaptive_branches_CV_{cv_prior}_BF_{BF_threshold}.fasta"
    shell:
        """
        cat {input.alignment_adaptive_final_nodes} {input.alignment_adaptive_internal_nodes} > {output.alignment_all_adaptive_branches}
        """