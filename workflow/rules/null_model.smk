"""
Null model for saltational evolution.

The null keeps the observed tree and the observed branch lengths (number of
substitutions per branch) fixed, and redraws only which substitutions occurred.

This module handles:
- Estimating the mutational spectrum used to draw the simulated substitutions,
  from 4-fold degenerate sites where all three alternative bases are synonymous
  and the estimate is therefore not filtered by amino-acid-level selection

Rules:
    - get_mutational_spectrum: Estimate the 12-way mutational spectrum
"""


rule get_mutational_spectrum:
    """
    Estimate the 12-way mutational spectrum from 4-fold degenerate sites.

    Degeneracy is determined per branch from the ancestral codon rather than from
    the reference codon, since mutations at the first two codon positions change
    whether the third position is 4-fold degenerate.
    """
    input:
        tree = rules.subset_viridian_tree.output.subset_tree,
        reference_seq_fasta = rules.params.input.usher_reference_fasta,
        reference_seq_gtf = rules.params.input.usher_reference_gtf,
    output:
        spectrum = "results/null_model/mutational_spectrum.tsv"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python scripts/mat-manipulation/get_mutational_spectrum.py \
            --tree {input.tree} \
            --reference_fasta {input.reference_seq_fasta} \
            --reference_gtf {input.reference_seq_gtf} \
            --output_spectrum {output.spectrum}
        """
