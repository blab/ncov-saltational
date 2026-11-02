rule get_viridian_alignment_saltational_branches:
    """Keep the Viridian sequences of the saltational branches (terminal only: internal nodes have no Viridian sequence)"""
    input:
        alignment_adaptive_final_nodes = rules.get_viridian_alignment_adaptive_branches.output.alignment_adaptive_final_nodes,
        name_saltational_branches = rules.name_saltational_branches.output.name_saltational_branches
    output:
        alignment_saltational_final_nodes = "results/nextclade/alignment_saltational_final_nodes_CV_{cv_prior}_BF_{BF_threshold}.fasta"
    shell:
        """
        awk 'FNR == NR {{ salt[$1] = 1; next }}
             /^>/ {{ keep = (substr($1, 2) in salt) }}
             keep' {input.name_saltational_branches} {input.alignment_adaptive_final_nodes} \
            > {output.alignment_saltational_final_nodes}
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
