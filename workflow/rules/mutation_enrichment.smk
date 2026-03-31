rule calculate_AA_pos_enrichment:
    """Calculate enrichment of AA mutations at each position within adaptive branches"""
    input:
        df_aa_muts_coded = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
        name_adaptive_branches = rules.concat_name_adaptive_branches.output.concat_name_adaptive_branches,
    output:
        df_AA_pos_enrichment = "results/mutation_enrichment/df_AA_pos_enrichment_BF_{BF_threshold}.tsv"
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
        df_AA_mut_enrichment = "results/mutation_enrichment/df_AA_mut_enrichment_BF_{BF_threshold}.tsv"
    shell:
        """
        Rscript ./scripts/mutation-enrichment/calculate_AA_mut_enrichment.R \
            --df_aa_muts {input.df_aa_muts_coded} \
            --name_adaptive_branches {input.name_adaptive_branches} \
            --output_df_enrichment {output.df_AA_mut_enrichment}
        """

