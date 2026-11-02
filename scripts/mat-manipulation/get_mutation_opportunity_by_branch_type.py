import argparse
import os
import sys
import tqdm
import pandas as pd

# Reuse the helpers of get_mutation_opportunity.py unchanged. Python puts the running
# script's directory on sys.path, but the Snakemake rule is run from the repository root,
# so the directory is added explicitly.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from get_mutation_opportunity import (read_tree, get_df_aa_muts, get_set_saltational_branches,
                                      get_stop_codons_index, build_translations_index,
                                      dict_AA_by_gene, POSITION_BATCH_SIZE)

# Mutation opportunity counts, with the saltational branches split into internal and
# terminal ones. The background stays pooled (all non-saltational branches), so that the
# stratified odds ratios keep the same denominator as the published analysis.

BRANCH_TYPES = ["saltational_internal", "saltational_terminal", "non_saltational"]


def get_branch_type(node_id, saltational_set):
    """Internal nodes are named node_XXX, terminal ones carry a sample name"""
    if node_id not in saltational_set:
        return "non_saltational"
    return "saltational_internal" if node_id.startswith("node_") else "saltational_terminal"


def count_opportunity_per_gene_by_branch_type(tree_expansion_from_root, translations_index, dict_gene_fasta,
                                              saltational_set, stop_codon_index, ns_muts_per_branch, gene, triples):
    """Same traversal as count_opportunity_per_gene in get_mutation_opportunity.py, with the
    saltational counts split by branch type.

    For each branch where opportunity exists (parent AA != aa_new AND pos < parent_min_stop),
    we add the branch's total NS mutation count (all genes) to the denominator of its group.
    """
    pos_to_aa_new = {}
    for pos, aa_new in triples:
        pos_to_aa_new.setdefault(pos, []).append(aa_new)

    counts = {t: {b: 0 for b in BRANCH_TYPES} for t in triples}
    positions = list(pos_to_aa_new.keys())
    translate = dict_gene_fasta[gene]['translate']
    n_batches = (len(positions) + POSITION_BATCH_SIZE - 1) // POSITION_BATCH_SIZE

    for batch_idx in range(n_batches):
        batch_pos = positions[batch_idx * POSITION_BATCH_SIZE : (batch_idx + 1) * POSITION_BATCH_SIZE]
        batch_triples = [(pos, aa_new) for pos in batch_pos for aa_new in pos_to_aa_new[pos]]
        ref_aa = {pos: translate[pos - 1] for pos in batch_pos}
        node_aa = {}

        desc = gene if n_batches == 1 else f"{gene} batch {batch_idx + 1}/{n_batches}"
        for node in tqdm.tqdm(tree_expansion_from_root, desc=desc):
            node_muts = translations_index.get(node.id, {}).get(gene, {})

            if node.parent is None:
                node_aa[node.id] = {pos: ref_aa[pos] for pos in batch_pos}
            else:
                parent_state = node_aa[node.parent.id]
                node_aa[node.id] = {pos: node_muts.get(pos, parent_state[pos]) for pos in batch_pos}

                node_min_stop = stop_codon_index.get((node.id, gene), float("inf"))
                branch_type = get_branch_type(node.id, saltational_set)
                n_ns = ns_muts_per_branch.get(node.id, 0)
                for pos, aa_new in batch_triples:
                    if pos <= node_min_stop and parent_state[pos] != aa_new:
                        counts[(pos, aa_new)][branch_type] += n_ns

    return counts


def get_df_opportunity_for_internal_terminal_split(
        df_aa_muts_coded, tree_expansion_from_root, translations_index,
        dict_gene_fasta, saltational_set, stop_codon_index, ns_muts_per_branch
    ):
    """Opportunity counts per mutation, split by branch type.

    n_saltational_opportunity is the sum of the internal and terminal columns, so the output
    is a superset of get_mutation_opportunity.py's and can be checked against it.
    """
    results = []

    for gene, gene_df in df_aa_muts_coded.groupby("gene"):
        triples = list(gene_df[["aa_muts_position", "aa_new"]].itertuples(index=False, name=None))
        counts = count_opportunity_per_gene_by_branch_type(
            tree_expansion_from_root, translations_index, dict_gene_fasta,
            saltational_set, stop_codon_index, ns_muts_per_branch, gene, triples
        )
        for pos, aa_new in triples:
            curr_counts = counts[(pos, aa_new)]
            results.append({
                "gene": gene,
                "aa_muts_position": pos,
                "aa_new": aa_new,
                "mutation": f"{gene}:{pos}{aa_new}",
                "n_saltational_opportunity_internal": curr_counts["saltational_internal"],
                "n_saltational_opportunity_terminal": curr_counts["saltational_terminal"],
                "n_saltational_opportunity": curr_counts["saltational_internal"] + curr_counts["saltational_terminal"],
                "n_non_saltational_opportunity": curr_counts["non_saltational"],
            })

    df_opportunity = pd.DataFrame(results)
    return df_opportunity


if __name__ == "__main__":

    # Load arguments
    parser = argparse.ArgumentParser(
        description="Get mutation opportunity counts, splitting saltational branches into internal and terminal ones",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help='UShER tree in protobuf format')
    parser.add_argument('--df_aa_muts', type=str, required=True, help='TSV file with amino acid mutations')
    parser.add_argument('--saltational_branches', type=str, required=True, help='TSV file with saltational branches names')
    parser.add_argument('--df_stop_codons', type=str, required=True, help='TSV file with stop codons')
    parser.add_argument('--reference_gtf', type=str, required=True, help='Reference GTF')
    parser.add_argument('--reference_fasta', type=str, required=True, help='Reference FASTA')
    parser.add_argument('--output_df_opportunity', type=str, required=True, help='Output TSV file with mutation opportunity counts')

    args = parser.parse_args()

    # Get viridian tree
    tree = read_tree(args.tree)

    # Get saltational branches
    saltational_set = get_set_saltational_branches(args.saltational_branches)

    # Get AA mutations and NS mut counts per branch
    df_aa_muts_coded, ns_muts_per_branch = get_df_aa_muts(args.df_aa_muts, saltational_set)

    # Get stop codons information
    stop_codon_index = get_stop_codons_index(args.df_stop_codons)

    # Get all the nodes in the tree (starting from the root)
    tree_expansion_from_root = tree.breadth_first_expansion()

    # Get AA muts across the tree
    translations = tree.translate(args.reference_gtf, args.reference_fasta)
    # Pre-index translations: {node_id: {gene: {pos: new_aa}}}
    translations_index = build_translations_index(translations)

    # Get AA translations for reference sequence
    dict_gene_fasta = dict_AA_by_gene(args.reference_gtf, args.reference_fasta)

    # Get mutation opportunity counts for each gene and mutation, split by branch type
    df_opportunity = get_df_opportunity_for_internal_terminal_split(
        df_aa_muts_coded, tree_expansion_from_root, translations_index,
        dict_gene_fasta, saltational_set, stop_codon_index, ns_muts_per_branch
    )

    # Write results to output file
    df_opportunity.to_csv(args.output_df_opportunity, sep = "\t", index = False)
