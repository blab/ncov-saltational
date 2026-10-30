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

rule calculate_AA_mut_enrichment_saltational_account_opportunity:
    """Calculate enrichment of AA mutations at each position within saltational branches"""
    input:
        df_aa_muts_coded = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
        name_saltational_branches = rules.name_saltational_branches.output.name_saltational_branches,
        df_opportunity = rules.get_mutation_opportunity_saltational.output.df_opportunity,
    output:
        df_AA_mut_enrichment = "results/mutation_enrichment/df_AA_mut_enrichment_saltational_opportunity_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        Rscript ./scripts/mutation-enrichment/calculate_AA_mut_enrichment_account_background.R \
            --df_aa_muts {input.df_aa_muts_coded} \
            --df_mutation_opportunity {input.df_opportunity} \
            --name_adaptive_branches {input.name_saltational_branches} \
            --output_df_enrichment {output.df_AA_mut_enrichment}
        """

# Enrichment restricted to saltational branches of a given branch type.
# Saltational events occur mostly on terminal branches. A virus on a terminal branch left no
# sampled descendant and may have adapted to within-host replication rather than to
# transmission, whereas an internal branch transmitted onward by definition. These rules
# recompute the mutation enrichment separately for saltational events on internal and on
# terminal branches, keeping all non-saltational branches as the background so that each
# stratified odds ratio stays comparable to the one computed over all saltational branches.

rule get_mutation_opportunity_by_branch_type:
    """
    Get mutation opportunity counts, splitting the saltational branches into internal and
    terminal ones. The background stays pooled across all non-saltational branches.
    """
    input:
        df_aa_muts_coded = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
        name_saltational_branches = rules.name_saltational_branches.output.name_saltational_branches,
        tree = rules.subset_viridian_tree.output.subset_tree,
        df_stop_codons = rules.get_stop_codons.output.df_stop_codons,
        reference_seq_gtf = rules.params.input.usher_reference_gtf,
        reference_seq_fasta = rules.params.input.usher_reference_fasta,
    output:
        df_opportunity = "results/mutation_enrichment/df_mutation_opportunity_by_branch_type_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    log:
        "logs/mutation_enrichment/get_mutation_opportunity_by_branch_type_CV_{cv_prior}_BF_{BF_threshold}.log"
    conda:
        "../envs/mutation-opportunity.yaml"
    shell:
        """
        python scripts/mat-manipulation/get_mutation_opportunity_by_branch_type.py \
            --df_aa_muts {input.df_aa_muts_coded} \
            --saltational_branches {input.name_saltational_branches} \
            --tree {input.tree} \
            --df_stop_codons {input.df_stop_codons} \
            --reference_gtf {input.reference_seq_gtf} \
            --reference_fasta {input.reference_seq_fasta} \
            --output_df_opportunity {output.df_opportunity} > {log} 2>&1
        """

rule calculate_AA_mut_enrichment_by_branch_type:
    """
    Calculate enrichment of AA mutations on saltational branches, separately for saltational
    events occurring on internal and on terminal branches.
    """
    input:
        df_aa_muts_coded = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
        name_saltational_branches = rules.name_saltational_branches.output.name_saltational_branches,
        df_opportunity = rules.get_mutation_opportunity_by_branch_type.output.df_opportunity,
    output:
        df_AA_mut_enrichment = "results/mutation_enrichment/df_AA_mut_enrichment_by_branch_type_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    log:
        "logs/mutation_enrichment/calculate_AA_mut_enrichment_by_branch_type_CV_{cv_prior}_BF_{BF_threshold}.log"
    shell:
        """
        Rscript ./scripts/mutation-enrichment/calculate_AA_mut_enrichment_by_branch_type.R \
            --df_aa_muts {input.df_aa_muts_coded} \
            --name_adaptive_branches {input.name_saltational_branches} \
            --df_mutation_opportunity {input.df_opportunity} \
            --output_df_enrichment {output.df_AA_mut_enrichment} > {log} 2>&1
        """

rule get_AA_reference_sequence:
    """Extract amino acid sequence from reference genome for each gene"""
    input:
        reference_seq_gtf = rules.params.input.usher_reference_gtf,
        reference_seq_fasta = rules.params.input.usher_reference_fasta,
    output:
        df_aa_reference_sequence = "results/mutation_enrichment/df_aa_reference_sequence.tsv"
    conda:
        "../envs/mutation-opportunity.yaml"
    shell:
        """
        python ./scripts/mutation-enrichment/get_AA_reference_sequence.py \
            --reference_gtf {input.reference_seq_gtf} \
            --reference_fasta {input.reference_seq_fasta} \
            --output_df_aa_reference_sequence {output.df_aa_reference_sequence}
        """

rule get_AA_frequency_nextstrain_clades:
    input:
        nextstrain_json = rules.params.input.nextstrain_json,
        df_aa_muts = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
        name_saltational_branches = rules.name_saltational_branches.output.name_saltational_branches,
    output:
        df_aa_freq_combined = "results/mutation_enrichment/df_aa_freq_combined_CV_{cv_prior}_BF_{BF_threshold}.tsv",
    shell:
        """
        python ./scripts/mutation-enrichment/get_AA_frequency_nextstrain_tree.py \
            --nextstrain_global_tree {input.nextstrain_json} \
            --df_aa_muts {input.df_aa_muts} \
            --saltational_branches {input.name_saltational_branches} \
            --output_df_combined {output.df_aa_freq_combined}
        """