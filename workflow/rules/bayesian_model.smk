"""
Bayesian model for detecting adaptive evolution.

This module handles:
- Preparing data for Bayesian model input
- Parametrizing priors based on tree-wide statistics
- Running posterior predictive checks
- Executing the Bayesian model to estimate selection parameters
- Computing Bayes Factors to identify adaptive branches

Rules:
    - prepare_data_model: Format mutation data for model input
    - get_prior_parametrization: Calculate prior parameters from data
    - run_posterior_predictive_checks: Validate model
    - run_bayesian_model: Estimate dN,dS and selection parameters
    - get_BF_across_model_runs: Compute Bayes Factors for branch selection (dN/dS > 1)
"""


rule prepare_data_model:
    """Prepare data for model run"""
    input:
        df_sites = rules.count_sites_per_branch.output.df_sites,
        df_aa_muts_coded = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
        name_recombinants = rules.get_branches_leading_to_recombinant_clades.output.name_strains_leading_to_recombinant_clades,
    output:
        data_viridian = "results/input_model/data_viridian_{gene}.tsv",
        data_viridian_grouped = "results/input_model/data_viridian_{gene}_grouped.tsv",
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


rule aggregate_genome_data:
    """Convert genome-wide gene-level data to aggregated format for downstream compatibility"""
    input:
        data_viridian_wide = "results/input_model/data_viridian_Genome.tsv"
    output:
        data_viridian_aggregated = "results/input_model/data_viridian_Genome_aggregated.tsv"
    shell:
        """
        Rscript ./scripts/bayesian-model/aggregate_genome_data.R \
            --data_viridian_wide {input.data_viridian_wide} \
            --output {output.data_viridian_aggregated}
        """


rule get_prior_parametrization:
    """Get prior parametrization for model run"""
    input:
        df_sites = rules.count_sites_per_branch.output.df_sites,
        df_aa_muts_coded = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
        df_parent_child_relationships = rules.get_parent_child_relationship.output.df_parent_relationships,
        name_recombinants = rules.get_branches_leading_to_recombinant_clades.output.name_strains_leading_to_recombinant_clades,
    output:
        df_prior_parametrization = "results/input_model/df_prior_parametrization.tsv"
    shell:
        """
        Rscript ./scripts/bayesian-model/get_prior_parametrization.R \
            --name_recombinants_to_remove {input.name_recombinants} \
            --df_sites {input.df_sites} \
            --df_aa_muts_coded {input.df_aa_muts_coded} \
            --df_parent_child {input.df_parent_child_relationships} \
            --df_prior_param {output.df_prior_parametrization}
        """


rule run_posterior_predictive_checks:
    """Run posterior predictive checks for each gene and each CV prior"""
    input:
        data_viridian_grouped = rules.prepare_data_model.output.data_viridian_grouped,
        df_prior_parametrization = rules.get_prior_parametrization.output.df_prior_parametrization,
    output:
        posterior_simul = "results/posterior_predictive_checks/posterior_predictive_check_{gene}_CV_{cv_prior}.tsv"
    log:
        "logs/posterior_predictive_checks/{gene}_CV_{cv_prior}.log"
    benchmark:
        "benchmarks/posterior_predictive_checks/{gene}_CV_{cv_prior}.tsv"
    threads: lambda wildcards: 8 if wildcards.gene == "Genome" else 1
    params:
        script = lambda wildcards: (
            "./scripts/bayesian-model/run_posterior_predictive_checks_genome.R"
            if wildcards.gene == "Genome"
            else "./scripts/bayesian-model/run_posterior_predictive_checks.R"
        ),
        n_cores_arg = lambda wildcards, threads: (
            f"--n_cores {threads}" if wildcards.gene == "Genome" else ""
        )
    shell:
        """
        Rscript {params.script} \
            --df_prior_param {input.df_prior_parametrization} \
            --genome_region {wildcards.gene} \
            --CV_prior {wildcards.cv_prior} \
            --data_viridian_grouped {input.data_viridian_grouped} \
            --output_posterior_simul {output.posterior_simul} \
            {params.n_cores_arg} > {log} 2>&1
        """


rule run_bayesian_model:
    """Run the Bayesian model for each gene and each CV prior"""
    input:
        data_viridian_grouped = rules.prepare_data_model.output.data_viridian_grouped,
        df_prior_parametrization = rules.get_prior_parametrization.output.df_prior_parametrization,
    output:
        df_param_estim = "results/output_model/df_param_estim_{gene}_CV_{cv_prior}.tsv"
    log:
        "logs/bayesian_model/run_bayesian_model_{gene}_CV_{cv_prior}.log"
    benchmark:
        "benchmarks/bayesian_model/run_bayesian_model_{gene}_CV_{cv_prior}.tsv"
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
            {params.n_cores_arg} > {log} 2>&1
        """


rule get_BF_across_model_runs:
    """Get Bayes Factors across model runs for each gene and each CV prior"""
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
