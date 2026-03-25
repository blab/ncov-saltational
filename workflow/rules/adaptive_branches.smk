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
        data_viridian = lambda wildcards: (
            "results/input_model/data_viridian_Genome_aggregated.tsv"
            if wildcards.gene == "Genome"
            else f"results/input_model/data_viridian_{wildcards.gene}.tsv"
        ),
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
        name_adaptive_branches_across_genes = lambda wildcards: (
            expand("results/adaptive_branches/name_adaptive_branches_{gene}_CV_{cv_prior}_BF_{BF_threshold}.tsv",
                   gene=[g for g in GENES if g != 'Genome'],
                   cv_prior=CV_prior_gene,
                   BF_threshold=wildcards.BF_threshold) +
            expand("results/adaptive_branches/name_adaptive_branches_{gene}_CV_{cv_prior}_BF_{BF_threshold}.tsv",
                   gene=['Genome'],
                   cv_prior=CV_prior_genome,
                   BF_threshold=wildcards.BF_threshold)
        ),
    output:
        concat_name_adaptive_branches = "results/adaptive_branches/concat_name_adaptive_branches_BF_{BF_threshold}.tsv"
    shell:
        """
        tsv-append -H {input.name_adaptive_branches_across_genes} | tsv-uniq >> {output.concat_name_adaptive_branches}
        """

