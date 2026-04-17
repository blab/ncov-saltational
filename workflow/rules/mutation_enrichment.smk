rule calculate_AA_pos_enrichment:
    """Calculate enrichment of AA mutations at each position within adaptive branches"""
    input:
        df_aa_muts_coded = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
        name_adaptive_branches = rules.concat_name_adaptive_branches.output.concat_name_adaptive_branches,
    output:
        df_AA_pos_enrichment = "results/mutation_enrichment/df_AA_pos_enrichment_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        Rscript ./scripts/mutation-enrichment/calculate_AA_pos_enrichment.R \
            --df_aa_muts {input.df_aa_muts_coded} \
            --name_adaptive_branches {input.name_adaptive_branches} \
            --output_df_enrichment {output.df_AA_pos_enrichment}
        """

rule calculate_AA_mut_enrichment:
    """Calculate enrichment of AA mutations at each position within adaptive branches"""
    input:
        df_aa_muts_coded = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
        name_adaptive_branches = rules.concat_name_adaptive_branches.output.concat_name_adaptive_branches,
    output:
        df_AA_mut_enrichment = "results/mutation_enrichment/df_AA_mut_enrichment_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        Rscript ./scripts/mutation-enrichment/calculate_AA_mut_enrichment.R \
            --df_aa_muts {input.df_aa_muts_coded} \
            --name_adaptive_branches {input.name_adaptive_branches} \
            --output_df_enrichment {output.df_AA_mut_enrichment}
        """

rule calculate_AA_pos_enrichment_saltational:
    """Calculate enrichment of AA mutations at each position within saltational branches"""
    input:
        df_aa_muts_coded = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
        name_saltational_branches = rules.name_saltational_branches.output.name_saltational_branches,
    output:
        df_AA_pos_enrichment = "results/mutation_enrichment/df_AA_pos_enrichment_saltational_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        Rscript ./scripts/mutation-enrichment/calculate_AA_pos_enrichment.R \
            --df_aa_muts {input.df_aa_muts_coded} \
            --name_adaptive_branches {input.name_saltational_branches} \
            --output_df_enrichment {output.df_AA_pos_enrichment}
        """

rule calculate_AA_mut_enrichment_saltational:
    """Calculate enrichment of AA mutations at each position within saltational branches"""
    input:
        df_aa_muts_coded = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
        name_saltational_branches = rules.name_saltational_branches.output.name_saltational_branches,
    output:
        df_AA_mut_enrichment = "results/mutation_enrichment/df_AA_mut_enrichment_saltational_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        Rscript ./scripts/mutation-enrichment/calculate_AA_mut_enrichment.R \
            --df_aa_muts {input.df_aa_muts_coded} \
            --name_adaptive_branches {input.name_saltational_branches} \
            --output_df_enrichment {output.df_AA_mut_enrichment}
        """

rule get_mutation_opportunity_saltational:
    """Get mutation opportunity counts for each gene and mutation"""
    input:
        df_aa_muts_coded = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
        name_saltational_branches = rules.name_saltational_branches.output.name_saltational_branches,
        tree = rules.subset_viridian_tree.output.subset_tree,
        df_stop_codons = rules.get_stop_codons.output.df_stop_codons,
        reference_seq_gtf = rules.params.input.usher_reference_gtf,
        reference_seq_fasta = rules.params.input.usher_reference_fasta,
    output:
        df_opportunity = "results/mutation_enrichment/df_mutation_opportunity_saltational_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    conda:
        "../envs/mutation-opportunity.yaml"
    shell:
        """
        python scripts/mat-manipulation/get_mutation_opportunity.py \
            --df_aa_muts {input.df_aa_muts_coded} \
            --saltational_branches {input.name_saltational_branches} \
            --tree {input.tree} \
            --df_stop_codons {input.df_stop_codons} \
            --reference_gtf {input.reference_seq_gtf} \
            --reference_fasta {input.reference_seq_fasta} \
            --output_df_opportunity {output.df_opportunity}
        """
