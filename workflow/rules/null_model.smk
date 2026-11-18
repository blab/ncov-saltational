"""
Null model for saltational evolution.

Each branch keeps its observed number of coding mutations in each gene; only the
synonymous / non-synonymous split is redrawn, with dN/dS fixed to each gene's
tree-wide average. The adaptive and saltational criteria are then applied exactly as
they are to the observed tree, giving the number of saltational branches expected
when no branch is adaptive.

Rules:
    - fixed_omega_null_saltational: Expected number of adaptive and saltational branches
"""


rule fixed_omega_null_saltational:
    """
    Expected number of adaptive and saltational branches when each branch keeps its
    number of mutations in each gene and only the synonymous / non-synonymous split is
    redrawn, with dN/dS fixed to each gene's tree-wide average.
    Uses the observed prior parametrization and the same criteria as the observed tree.
    """
    input:
        data_viridian = expand("results/input_model/data_viridian_{gene}.tsv", gene=GENES_PPC),
        df_prior_parametrization = rules.get_prior_parametrization.output.df_prior_parametrization,
    output:
        replicates = "results/null_model/null_replicates_CV_{cv_prior}_BF_{BF_threshold}.tsv",
        summary = "results/null_model/null_summary_CV_{cv_prior}_BF_{BF_threshold}.tsv",
        dir_branch_counts = directory("results/null_model/branch_counts_CV_{cv_prior}_BF_{BF_threshold}")
    params:
        dir_input_model = "results/input_model",
        min_nonsyn = 5,
        n_replicates = N_REPLICATES_FIXED_OMEGA_NULL,
        seed = NULL_BASE_SEED
    threads: 6
    log:
        "logs/null_model/fixed_omega_null_saltational_CV_{cv_prior}_BF_{BF_threshold}.log"
    shell:
        """
        Rscript ./scripts/null-model/fixed_omega_null_saltational.R \
            --dir_input_model {params.dir_input_model} \
            --df_prior_param {input.df_prior_parametrization} \
            --CV_prior {wildcards.cv_prior} \
            --BF_threshold {wildcards.BF_threshold} \
            --min_nonsyn {params.min_nonsyn} \
            --n_replicates {params.n_replicates} \
            --seed {params.seed} \
            --n_cores {threads} \
            --output_replicates {output.replicates} \
            --output_summary {output.summary} \
            --dir_output_branch_counts {output.dir_branch_counts} > {log} 2>&1
        """
