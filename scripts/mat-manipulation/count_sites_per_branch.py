#!/usr/bin/env python3
"""
Count the number of possible syn and nonsyn mutations available per gene
on each branch of the MAT tree
Author: Katie Kistler
"""

import argparse
import json
import random
import time
from itertools import product

import bte
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def readin_reference(fasta_path, gtf_path):
    """ """
    # read in reference sequence and gene locations
    with open(fasta_path, "r") as handle:
        for seq_record in SeqIO.parse(handle, "fasta"):
            reference_genome_seq = seq_record.seq

    # read in gene annotation GTF file, and get start and end nts for each gene
    gene_coords = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            start, end = int(fields[3]), int(fields[4])
            gene = fields[8].split('gene_id "')[1].split('"')[0]
            gene_coords[gene] = {"start": start, "end": end}

    return reference_genome_seq, gene_coords


def count_syn_nonsyn_mutations(codon):
    """
    Of the 9 possible one step mutations, return the number that are syn or nonsyn
    """
    original_aa = Seq(codon).translate()
    synonymous = nonsynonymous = 0

    for i in range(3):
        for base in "ACGT":
            if base != codon[i]:
                mutated = codon[:i] + base + codon[i + 1 :]
                if Seq(mutated).translate() == original_aa:
                    synonymous += 1
                else:
                    nonsynonymous += 1

    return nonsynonymous, synonymous


def compute_baseline_counts_per_gene(reference_seq, gene_coords, mut_counts_mapper):
    """
    Compute baseline syn/nonsyn mutation counts for each gene in the reference sequence.

    Return a dict: {gene_name: {'syn': count, 'nonsyn': count}}
    """
    baseline_counts = {}

    for gene, coords in gene_coords.items():
        # Extract CDS sequence (convert 1-based to 0-based indexing)
        start_idx = coords["start"] - 1
        end_idx = coords["end"]  # end is inclusive in GTF, so this works for slicing
        gene_seq = str(reference_seq[start_idx:end_idx])

        # Split into codons
        num_codons = len(gene_seq) // 3

        total_syn = 0
        total_nonsyn = 0

        for i in range(num_codons):
            codon = gene_seq[i * 3 : (i + 1) * 3]
            if len(codon) == 3:  # Ensure complete codon
                nonsyn, syn = mut_counts_mapper[codon]
                total_syn += syn
                total_nonsyn += nonsyn
            else:
                print("error divying codons")

        baseline_counts[gene] = {"syn": total_syn, "nonsyn": total_nonsyn}

    return baseline_counts


def update_counts_for_mutation(
    gene_name, position, new_base, current_seq, gene_coords, mut_counts_mapper
):
    """
    Calculate the change in syn/nonsyn counts when a mutation occurs in a gene.
    """
    coords = gene_coords[gene_name]
    gene_start = coords["start"]
    gene_end = coords["end"]

    # Check if mutation is within this gene
    if position < gene_start or position > gene_end:
        print("error, position not in gene")

    # Position within the gene (0-based)
    pos_in_gene = position - gene_start

    # Which codon is affected (0-based codon index)
    codon_idx = pos_in_gene // 3

    # Position within the codon (0, 1, or 2)
    pos_in_codon = pos_in_gene % 3

    # Extract the gene sequence
    gene_seq = str(current_seq[gene_start - 1 : gene_end])

    # Get the old codon
    old_codon = gene_seq[codon_idx * 3 : (codon_idx + 1) * 3]

    if len(old_codon) != 3:
        print("codon counting error")

    # Create new codon
    new_codon = old_codon[:pos_in_codon] + new_base + old_codon[pos_in_codon + 1 :]

    # Get old and new counts
    old_nonsyn, old_syn = mut_counts_mapper[old_codon]
    new_nonsyn, new_syn = mut_counts_mapper[new_codon]

    # Return the delta
    return {"syn": new_syn - old_syn, "nonsyn": new_nonsyn - old_nonsyn}


def find_premature_stop(gene_seq, num_codons):
    """
    Scan gene for premature stop codon
    Returns the 0-based codon index of the first premature stop, or None.
    """
    stop_codons = {"TAA", "TAG", "TGA"}
    for i in range(num_codons - 1):  # exclude last codon (the normal stop)
        codon = gene_seq[i * 3 : (i + 1) * 3]
        if codon in stop_codons:
            return i
    return None


def recompute_counts_truncated(gene_seq, stop_codon_idx, mut_counts_mapper):
    """
    Compute syn/nonsyn counts for codons 0 to stop_codon_idx - 1 (before the stop).
    """
    total_syn = 0
    total_nonsyn = 0
    for i in range(stop_codon_idx):
        codon = gene_seq[i * 3 : (i + 1) * 3]
        if len(codon) == 3:
            nonsyn, syn = mut_counts_mapper[codon]
            total_syn += syn
            total_nonsyn += nonsyn
    return {"syn": total_syn, "nonsyn": total_nonsyn}


def get_affected_genes(position, gene_coords):
    """
    Find which gene(s) contain the mutated position (handles overlaps).
    """
    affected_gene = []
    for gene, coords in gene_coords.items():
        if coords["start"] <= position <= coords["end"]:
            affected_gene.append(gene)
    return affected_gene


def parse_mutation(mut_string):
    """
    Parse mutation string like 'T28144C' into (position, orig_nt, new_nt)
    """
    mut_from = mut_string[0]
    mut_to = mut_string[-1]
    position = int(mut_string[1:-1])

    return position, mut_from, mut_to


def read_tree(tree_file):
    tree = bte.MATree(tree_file)
    return tree


def compute_node_counts(
    tree,
    all_nodes,
    reference_seq,
    gene_coords,
    baseline_counts,
    mut_counts_mapper,
    results_file_name,
    truncate_at_stop=False,
):
    """
    Compute available syn/nonsyn mutations per gene at each node in the tree.
    Save results as json
    """
    # node_counts = {}
    rows_for_df = []

    for node in all_nodes:
        # Start with a copy of the reference sequence
        current_seq = str(reference_seq)

        # Start with baseline counts (deep copy)
        current_counts = {
            gene: counts.copy() for gene, counts in baseline_counts.items()
        }

        # Get path from root to this node
        node_path = tree.rsearch(node.id, include_self=True)

        # Apply mutations along the path (root to tip order)
        for ancestor in reversed(node_path):
            for mut_string in ancestor.mutations:
                # Parse the mutation
                position, ref_base, new_base = parse_mutation(mut_string)
                # Find which gene(s) are affected
                affected_genes = get_affected_genes(position, gene_coords)

                # Skip if mutation is outside all coding regions
                if len(affected_genes) > 0:
                    # Update counts for each affected gene
                    for gene in affected_genes:
                        delta = update_counts_for_mutation(
                            gene,
                            position,
                            new_base,
                            current_seq,
                            gene_coords,
                            mut_counts_mapper,
                        )
                        current_counts[gene]["syn"] += delta["syn"]
                        current_counts[gene]["nonsyn"] += delta["nonsyn"]

                    # Apply mutation to sequence
                    # this will correctly account for multiple mutations at same codon along path
                    current_seq = (
                        current_seq[: position - 1] + new_base + current_seq[position:]
                    )

        # If we're looking for premature stops, check each branch for them and then truncate counts
        early_stop_genes = []
        if truncate_at_stop:
            for gene, coords in gene_coords.items():
                start_idx = coords["start"] - 1
                end_idx = coords["end"]
                gene_seq = str(current_seq[start_idx:end_idx])
                num_codons = len(gene_seq) // 3
                # will return an index if it finds premature stop, otherwise None
                stop_idx = find_premature_stop(gene_seq, num_codons)
                if stop_idx is not None:
                    early_stop_genes.append(gene)
                    current_counts[gene] = recompute_counts_truncated(
                        gene_seq, stop_idx, mut_counts_mapper
                    )
        # genes where we found premature stops
        early_stop_str = ",".join(early_stop_genes)

        # Store the final counts for this node
        for g, cs in current_counts.items():
            row = {
                "node_id": node.id,
                "gene": g,
                "syn_sites": cs["syn"],
                "nonsyn_sites": cs["nonsyn"],
            }
            if truncate_at_stop:
                row["early_stop"] = early_stop_str
            rows_for_df.append(row)

    # put results in tidy dataframe
    df = pd.DataFrame(rows_for_df)

    # save as tsv
    with open(results_file_name, "w") as f:
        df.to_csv(f, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Count possible syn and nonsyn mutations available on each branch of tree"
    )
    parser.add_argument(
        "--ref-fasta",
        help="Path to reference fasta sequence",
    )
    parser.add_argument(
        "--ref-gtf",
        help="Path to reference GTF annotations",
    )
    parser.add_argument(
        "--tree-path",
        help="Path to the viridian MAT tree",
    )
    parser.add_argument(
        "--results-file",
        help="Path to save the results TSV",
    )
    parser.add_argument(
        "--truncate-at-stop",
        action="store_true",
        default=False,
        help="If set, truncate gene site counts at premature stop codons",
    )

    args = parser.parse_args()

    # get reference sequence and gene coordinates
    reference_genome_seq, gene_coords = readin_reference(args.ref_fasta, args.ref_gtf)

    # make a dictionary mapper saying how many syn and nonsyn available via single mutation from each starting codon
    all_codons = ["".join(c) for c in product("ACGT", repeat=3)]

    mut_counts_codon_mapper = {c: count_syn_nonsyn_mutations(c) for c in all_codons}

    # get the total number of syn and nonsyn for the reference sequence (also the tree root)
    # do this per gene
    baseline_counts = compute_baseline_counts_per_gene(
        reference_genome_seq, gene_coords, mut_counts_codon_mapper
    )

    # load the tree
    tree = read_tree(args.tree_path)

    # get all node all node names
    all_nodes = tree.depth_first_expansion()
    # if you want to test with just a small number of branches
    # test_nodes = random.sample(all_nodes, 50)

    start_time = time.time()
    # count available muts for all nodes
    compute_node_counts(
        tree,
        all_nodes,
        reference_genome_seq,
        gene_coords,
        baseline_counts,
        mut_counts_codon_mapper,
        args.results_file,
        args.truncate_at_stop,
    )
    end_time = time.time()

    print(f"Took {end_time - start_time:.2f} seconds")
