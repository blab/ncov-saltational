rule get_viridian_alignment_saltational_branches:
    """Download the Viridian sequences of the saltational branches (terminal only: internal nodes have no Viridian sequence)"""
    input:
        name_saltational_branches = rules.name_saltational_branches.output.name_saltational_branches
    output:
        alignment_saltational_final_nodes = "results/nextclade/alignment_saltational_final_nodes_CV_{cv_prior}_BF_{BF_threshold}.fasta"
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
            {input.name_saltational_branches} {output.alignment_saltational_final_nodes}
        """

rule run_nextclade_saltational_branches:
    """Run Nextclade on the sequences corresponding to the saltational branches"""
    input:
        sequences_saltational_branches = rules.get_viridian_alignment_saltational_branches.output.alignment_saltational_final_nodes
    params:
        nextclade_folder = "data/nextclade"
    output:
        nextclade_saltational_branches = "results/nextclade/nextclade_saltational_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv"
    shell:
        """
        nextclade run \
            --input-dataset={params.nextclade_folder} \
            --output-tsv={output.nextclade_saltational_branches} \
            {input.sequences_saltational_branches}
        """
