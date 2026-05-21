rule run_nextclade_adaptive_branches:
    """Run Nextclade on the sequences corresponding to the adaptive branches"""
    input:
        sequences_adaptive_branches = rules.get_alignment_all_adaptive_branches.output.alignment_all_adaptive_branches
    params:
        nextclade_folder = "data/nextclade"
    output:
        nextclade_adaptive_branches = "results/nextclade/nextclade_adaptive_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        nextclade run \
            --input-dataset={params.nextclade_folder} \
            --output-tsv={output.nextclade_adaptive_branches} \
            {input.sequences_adaptive_branches}
        """
