import bte
import argparse
import sys
import os
import tqdm
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from gtfparse import read_gtf

def read_tree(tree_file):
    tree = bte.MATree(tree_file)
    return tree


def get_df_aa_muts(aa_muts_file, saltational_set):
    df_aa_muts = pd.read_csv(aa_muts_file, sep = "\t")

    df_aa_muts_coded = df_aa_muts[
        (df_aa_muts["is_coded"]) &
        (df_aa_muts["aa_initial"] != df_aa_muts["aa_new"]) &
        (df_aa_muts["node_id"].isin(saltational_set))
    ][["gene", "aa_muts_position", "aa_new"]].drop_duplicates().copy()

    ns_muts_per_branch = (
        df_aa_muts[df_aa_muts["is_coded"] & (df_aa_muts["aa_initial"] != df_aa_muts["aa_new"])]
        .groupby("node_id")
        .size()
        .to_dict()
    )

    return df_aa_muts_coded, ns_muts_per_branch

def get_set_saltational_branches(salt_branches_file):
    strains = []

    with open(salt_branches_file, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0]:
                strains.append(parts[0])
    return set(strains)

def get_stop_codons_index(stop_codons_file):
    df_stop_codons = pd.read_csv(stop_codons_file, sep="\t")

    stop_codon_index = (
        df_stop_codons
        .groupby(["node_id", "gene"])["position"]
        .min()
        .to_dict()
    )
    return stop_codon_index

def build_translations_index(translations):
    """Pre-index translations into {node_id: {gene: {position: new_aa}}} """
    index = {}
    for node_id, muts in translations.items():
        gene_dict = {}
        for aa_mut in muts:
            aa_str = aa_mut.aa
            pos = int(aa_str[1:-1])
            new_aa = aa_str[-1]
            gene_dict.setdefault(aa_mut.gene, {})[pos] = new_aa
        index[node_id] = gene_dict
    return index

def get_dict_gene_coordinate_from_gtf(reference_gtf_path):
    # Load GTF
    reference_gtf = read_gtf(reference_gtf_path).to_pandas()
    # Get list of gene names
    list_gene_names = list(reference_gtf['gene_id'].unique())
    
    # Get dictionary of gene coordinates
    dict_gene_coordinates = {}
    for gene in list_gene_names:
        start_start_codon = int(reference_gtf[(reference_gtf['gene_id'] == gene) & (reference_gtf['feature'] == "start_codon")]['start'].iloc[0])
        end_start_codon = int(reference_gtf[(reference_gtf['gene_id'] == gene) & (reference_gtf['feature'] == "start_codon")]['end'].iloc[0])
        start_end_codon = int(reference_gtf[(reference_gtf['gene_id'] == gene) & (reference_gtf['feature'] == "stop_codon")]['start'].iloc[0])
        end_end_codon = int(reference_gtf[(reference_gtf['gene_id'] == gene) & (reference_gtf['feature'] == "stop_codon")]['end'].iloc[0])
        dict_gene_coordinates[gene] = {'start_codon': [start_start_codon, end_start_codon], 'stop_codon': [start_end_codon, end_end_codon]}
    
    return dict_gene_coordinates

def subset_sequence(seq, start_position, stop_position):
    # start_position and stop_position are >= 1
    nuc_string = str(seq.seq)[(start_position - 1):stop_position]
    return Seq(nuc_string)

def dict_AA_by_gene(reference_gtf_path, reference_alignment_path):
    # Read alignment fasta
    reference_seq = SeqIO.read(reference_alignment_path, "fasta")

    # Get gene coordinates
    dict_gene_coordinates = get_dict_gene_coordinate_from_gtf(reference_gtf_path)

    # Get dictionary with gene alignments
    dict_gene_fasta = {}
    for gene, gene_dict in dict_gene_coordinates.items():
        nuc_alignment = subset_sequence(reference_seq, gene_dict['start_codon'][0], gene_dict['stop_codon'][1])
        nuc_translate = nuc_alignment.translate()
        dict_gene_fasta[gene] = {'alignment': nuc_alignment, 'translate': nuc_translate}
        
    return dict_gene_fasta


POSITION_BATCH_SIZE = 200


def count_opportunity_per_gene(tree_expansion_from_root, translations_index, dict_gene_fasta,
                               saltational_set, stop_codon_index, ns_muts_per_branch, gene, triples):
    """Traversal per batch of positions for one gene. Batching caps memory to POSITION_BATCH_SIZE
    positions at a time, keeping node_aa bounded regardless of how many positions the gene has.

    For each branch where opportunity exists (parent AA != aa_new AND pos < parent_min_stop),
    we add the branch's total NS mutation count (all genes) to the denominator. This gives
    "NS mutations on branches with opportunity" — the correct denominator for the Fisher's test.
    """
    pos_to_aa_new = {}
    for pos, aa_new in triples:
        pos_to_aa_new.setdefault(pos, []).append(aa_new)

    counts = {t: {"saltational": 0, "non_saltational": 0} for t in triples}
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

                parent_min_stop = stop_codon_index.get((node.parent.id, gene), float("inf"))
                branch_type = "saltational" if node.id in saltational_set else "non_saltational"
                n_ns = ns_muts_per_branch.get(node.id, 0)
                for pos, aa_new in batch_triples:
                    if pos < parent_min_stop and parent_state[pos] != aa_new:
                        counts[(pos, aa_new)][branch_type] += n_ns

    return counts


def get_df_opportunity_by_gene(
        df_aa_muts_coded, tree_expansion_from_root, translations_index,
        dict_gene_fasta, saltational_set, stop_codon_index, ns_muts_per_branch
    ):
    results = []

    for gene, gene_df in df_aa_muts_coded.groupby("gene"):
        triples = list(gene_df[["aa_muts_position", "aa_new"]].itertuples(index=False, name=None))
        counts = count_opportunity_per_gene(
            tree_expansion_from_root, translations_index, dict_gene_fasta,
            saltational_set, stop_codon_index, ns_muts_per_branch, gene, triples
        )
        for pos, aa_new in triples:
            results.append({
                "gene": gene,
                "aa_muts_position": pos,
                "aa_new": aa_new,
                "mutation": f"{gene}:{pos}{aa_new}",
                "n_saltational_opportunity": counts[(pos, aa_new)]["saltational"],
                "n_non_saltational_opportunity": counts[(pos, aa_new)]["non_saltational"],
            })

    df_opportunity = pd.DataFrame(results)
    return df_opportunity


if __name__ == "__main__":

    # Load arguments
    parser = argparse.ArgumentParser(
        description="Get parent / child node relationships",
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

    # Get mutation opportunity counts for each gene and mutation
    df_opportunity = get_df_opportunity_by_gene(
        df_aa_muts_coded, tree_expansion_from_root, translations_index,
        dict_gene_fasta, saltational_set, stop_codon_index, ns_muts_per_branch
    )

    # Write results to output file
    df_opportunity.to_csv(args.output_df_opportunity, sep = "\t", index = False)
