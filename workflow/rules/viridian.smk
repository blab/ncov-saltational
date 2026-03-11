"""
Viridian tree processing and temporal inference.

This module handles:
- Converting UShER protobuf tree to Newick format
- Extracting metadata from Taxonium files
- Running Chronumental for temporal inference

Rules:
    - convert_tree_newick: Convert tree format for Chronumental
    - get_metadata_from_taxonium: Extract and prepare metadata
    - chronumental: Infer dates using molecular clock
"""


rule convert_tree_newick:
    """Convert UShER tree to newick format - This will be used to run Chronumental"""
    input:
        usher_tree = rules.params.input.viridian_tree
    output:
        tree_newick = "results/viridian/viridian_tree.nwk"
    conda:
        "workflow/envs/bte.yaml"
    shell:
        """
        matUtils extract \
            --input-mat {input.usher_tree} \
            --write-tree {output.tree_newick}
        """


rule get_metadata_from_taxonium:
    """Get metadata stored in viridian tree and add reference sequence information"""
    input:
        taxonium = rules.params.input.viridian_taxonium
    output:
        metadata = "results/viridian/metadata.tsv"
    params:
        name_reference = "MN908947",
        date_reference = "2019-12-26", # hospital admission date
        country_reference = "China",
        pango_reference = "B"
    conda:
        "workflow/envs/bte.yaml"
    shell:
        """
        python3 scripts/mat-manipulation/get_metadata_from_taxonium.py \
            --input_taxonium {input.taxonium} \
            --output_metadata {output.metadata} \
            --name_reference {params.name_reference} \
            --date_reference {params.date_reference} \
            --country_reference {params.country_reference} \
            --pango_reference {params.pango_reference}
        """


rule chronumental:
    """Run chronumental to estimate time tree"""
    input:
        tree = rules.convert_tree_newick.output.tree_newick,
        metadata = rules.get_metadata_from_taxonium.output.metadata
    output:
        inferred_dates = "results/chronumental/chronumental_inferred_dates.tsv",
        timetree = "results/chronumental/chronumental_time_tree.pb.gz"
    params:
        reference_node = "MN908947"
    conda:
        "workflow/envs/bte.yaml"
    shell:
        """
        chronumental \
            --tree {input.tree} \
            --dates {input.metadata} \
            --dates_out {output.inferred_dates} \
            --tree_out {output.timetree} \
            --steps 2000
        """
