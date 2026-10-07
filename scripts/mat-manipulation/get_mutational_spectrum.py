import bte
import argparse
import pandas as pd
import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from itertools import product
from collections import Counter

# Note: read_tree, readin_reference, parse_mutation and get_affected_genes are
# duplicated from count_sites_per_branch.py. That file is not an importable
# package and cross-script imports would depend on the working directory, so the
# functions are copied here (as read_tree already is across several scripts).


def read_tree(tree_file):
    tree = bte.MATree(tree_file)
    return tree


def readin_reference(fasta_path, gtf_path):
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


def parse_mutation(mut_string):
    """
    Parse mutation string like 'T28144C' into (position, orig_nt, new_nt)
    """
    mut_from = mut_string[0]
    mut_to = mut_string[-1]
    position = int(mut_string[1:-1])

    return position, mut_from, mut_to


def is_fourfold(codon):
    """Third position is 4-fold degenerate if all 4 bases there give the same AA."""
    aa = str(Seq(codon).translate())
    if aa == '*':
        return False
    return len({str(Seq(codon[:2] + b).translate()) for b in 'ACGT'}) == 1


def get_fourfold_mapper():
    all_codons = ["".join(c) for c in product("ACGT", repeat=3)]
    fourfold_mapper = {c: is_fourfold(c) for c in all_codons}
    return fourfold_mapper


def get_dict_mut_from_ref(tree):
    """This returns a dictionary storing the nucleotide site that are mutated compared with the reference"""
    dict_mut_from_ref = {}

    for node in tqdm.tqdm(tree.depth_first_expansion()): # Tree traversal from root to tips
        if node.parent is None:
            dict_mut_from_ref[node.id] = {}
        else:
            dict_mut_from_ref[node.id] = dict(dict_mut_from_ref[node.parent.id])   # Copy parent dictionary
        for m in node.mutations:
            pos, a, b = parse_mutation(m)
            dict_mut_from_ref[node.id][pos] = b

    return dict_mut_from_ref


def get_affected_genes(position, gene_coords):
    """
    Find which gene(s) contain the mutated position (handles overlaps).
    """
    affected_gene = []
    for gene, coords in gene_coords.items():
        if coords["start"] <= position <= coords["end"]:
            affected_gene.append(gene)
    return affected_gene


def get_countable_codon_position(position, gene_coords):
    """Computes whether a susbtitution at `position` is a candidate for the 4-fold spectrum:
    (i) it must fall in exactly one gene (no overlapping frames)
    (ii) it must be the third position of a codon. It returns (gene, codon_start) if countable, else None."""

    # Get genes containing the position
    affected_genes = get_affected_genes(position, gene_coords)
    if len(affected_genes) != 1:
        return None
    else:
        gene = affected_genes[0]
        gene_start = gene_coords[gene]["start"]
        if (position - gene_start) % 3 != 2:
            return None
        else:
            return gene, position - 2


def get_codon_corresponding_to_site(codon_start, node_id, dict_mut_from_ref, reference_genome_seq):
    """
    Return the 3-base codon beginning at `codon_start` (1-based genomic position),
    as it exists on `node_id` - i.e. the reference sequence with that node's
    accumulated mutations applied.
    """
    # Get substitutions that occurred at node
    mut_at_node_from_ref = dict_mut_from_ref[node_id]
    #
    return "".join(
        mut_at_node_from_ref.get(p, reference_genome_seq[p - 1])
        for p in (codon_start, codon_start + 1, codon_start + 2)
    )


def get_spectrum_from_fourfold_degenerate_sites(tree, dict_mut_from_ref, reference_genome_seq, gene_coords):
    # Initialize spectrum counter
    spectrum = Counter()

    # Get reference sequence
    ref = str(reference_genome_seq)
    # Get four-fold mapper
    fourfold_mapper = get_fourfold_mapper()

    # Tree traversal and count relevant substitutions across tree
    for node in tqdm.tqdm(tree.depth_first_expansion()):
        if node.parent is None:
            continue # The root is not a branch
        for m in node.mutations:
            pos, a, b = parse_mutation(m)
            hit = get_countable_codon_position(pos, gene_coords)
            if hit is None:
                continue
            gene, codon_start = hit
            codon = get_codon_corresponding_to_site(codon_start, node.parent.id, dict_mut_from_ref, ref)

            if fourfold_mapper.get(codon, False):
                spectrum[(a, b)] += 1

    return spectrum


def get_df_spectrum(spectrum):
    """Convert the spectrum Counter to a tidy dataframe.

    raw_proportion is the type-first SAMPLING WEIGHT (share of counts), not a
    per-site mutation rate: raw(a->b) is proportional to rate(a->b) * n_sites(a),
    so it already accounts for how many sites carry each starting base.
    """
    total = sum(spectrum.values())
    df_spectrum = pd.DataFrame(
        [{"from_base": a, "to_base": b, "count": n, "raw_proportion": n / total}
         for (a, b), n in sorted(spectrum.items(), key=lambda kv: -kv[1])]
    )

    return df_spectrum


if __name__ == "__main__":

    # Load arguments
    parser = argparse.ArgumentParser(
        description="Get the mutational spectrum from 4-fold degenerate sites across the tree",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help='UShER tree in protobuf format')
    parser.add_argument('--reference_fasta', type=str, required=True, help='Reference sequence in FASTA format')
    parser.add_argument('--reference_gtf', type=str, required=True, help='Reference GTF annotations')
    parser.add_argument('--output_spectrum', type=str, required=True, help='TSV output file with the mutational spectrum')
    args = parser.parse_args()

    tree = read_tree(args.tree)

    reference_genome_seq, gene_coords = readin_reference(args.reference_fasta, args.reference_gtf)

    # Get the mutations carried by each node relative to the reference
    dict_mut_from_ref = get_dict_mut_from_ref(tree)

    # Count substitutions at sites that are 4-fold degenerate on their own branch
    spectrum = get_spectrum_from_fourfold_degenerate_sites(
        tree, dict_mut_from_ref, reference_genome_seq, gene_coords)

    df_spectrum = get_df_spectrum(spectrum)

    # Write results to output file
    df_spectrum.to_csv(args.output_spectrum, sep = "\t", index = False)
