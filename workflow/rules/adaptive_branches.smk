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

