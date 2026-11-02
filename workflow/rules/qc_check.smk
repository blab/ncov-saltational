"""
Quality control check on terminal saltational branches.
"""


rule get_qc_regions_from_nextclade:
    """
    Get the regions of each saltational branch where an amino acid mutation call
    may be an artefact rather than a genuine mutation based on Nextclade annotations.
    (Conservative)
    """
    input:
        nextclade = rules.run_nextclade_saltational_branches.output.nextclade_saltational_branches,
        gene_coordinates = "figures/param/proteins_aa_coordinates.csv",
    output:
        qc_regions = "results/qc_check/qc_regions_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    log:
        "logs/qc_check/get_qc_regions_from_nextclade_CV_{cv_prior}_BF_{BF_threshold}.log"
    shell:
        """
        Rscript ./scripts/qc-check/get_qc_regions_from_nextclade.R \
            --nextclade {input.nextclade} \
            --gene_coordinates {input.gene_coordinates} \
            --output_qc_regions {output.qc_regions} > {log} 2>&1
        """


rule flag_qc_affected_mutations:
    """
    Identify the amino acid mutations of the saltational branches that fall in a flagged
    region, and the mutation counts that remain once they are removed.
    """
    input:
        qc_regions = rules.get_qc_regions_from_nextclade.output.qc_regions,
        df_aa_muts_coded = rules.add_coded_status_aa_mut.output.df_aa_muts_coded,
        name_saltational_branches = rules.name_saltational_branches.output.name_saltational_branches,
        data_viridian = expand("results/input_model/data_viridian_{gene}.tsv", gene=GENES_PPC),
    output:
        flagged_mutations = "results/qc_check/qc_flagged_mutations_CV_{cv_prior}_BF_{BF_threshold}.tsv",
        counts = "results/qc_check/df_counts_qc_filtered_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    params:
        dir_input_model = "results/input_model"
    log:
        "logs/qc_check/flag_qc_affected_mutations_CV_{cv_prior}_BF_{BF_threshold}.log"
    shell:
        """
        Rscript ./scripts/qc-check/flag_qc_affected_mutations.R \
            --qc_regions {input.qc_regions} \
            --df_aa_muts {input.df_aa_muts_coded} \
            --name_saltational_branches {input.name_saltational_branches} \
            --dir_input_model {params.dir_input_model} \
            --output_flagged_mutations {output.flagged_mutations} \
            --output_counts {output.counts} > {log} 2>&1
        """


rule recompute_saltational_after_qc:
    """
    Recompute the Bayes Factors and the saltational status of the affected branches using
    the filtered mutation counts. Unaffected branches keep their published status.
    """
    input:
        counts = rules.flag_qc_affected_mutations.output.counts,
        df_prior_parametrization = rules.get_prior_parametrization.output.df_prior_parametrization,
        name_saltational_branches = rules.name_saltational_branches.output.name_saltational_branches,
    output:
        name_saltational_branches = "results/qc_check/name_saltational_branches_qc_filtered_CV_{cv_prior}_BF_{BF_threshold}.tsv",
        summary = "results/qc_check/qc_check_summary_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    params:
        min_nonsyn = 5,
        n_draws_genome = N_DRAWS_GENOME,
        seed = 1
    threads: 6
    log:
        "logs/qc_check/recompute_saltational_after_qc_CV_{cv_prior}_BF_{BF_threshold}.log"
    shell:
        """
        Rscript ./scripts/qc-check/recompute_saltational_after_qc.R \
            --df_counts {input.counts} \
            --df_prior_param {input.df_prior_parametrization} \
            --name_saltational_branches {input.name_saltational_branches} \
            --CV_prior {wildcards.cv_prior} \
            --BF_threshold {wildcards.BF_threshold} \
            --min_nonsyn {params.min_nonsyn} \
            --n_draws_genome {params.n_draws_genome} \
            --seed {params.seed} \
            --n_cores {threads} \
            --output_name_saltational_branches {output.name_saltational_branches} \
            --output_summary {output.summary} > {log} 2>&1
        """
