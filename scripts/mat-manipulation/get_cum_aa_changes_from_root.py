import bte
import argparse
import sys
import re
import pandas as pd
import tqdm

def read_tree(tree_file):
    tree = bte.MATree(tree_file)
    return tree

def read_name_saltational_branches(salt_branches_file):
    strains = []

    with open(salt_branches_file, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0]:
                strains.append(parts[0])
    return strains

def is_non_synonymous(aa_string):
    aa_pattern = re.compile(r"^([A-Z*])\d+([A-Z*])$")
    match = aa_pattern.match(aa_string)
    aa_init = match.groups()[0]
    aa_new = match.groups()[1]
    return aa_init != aa_new

def get_cum_aa_changes_from_root(tree_expansion_from_root, translations, vec_name_saltational_branches):

    dict_aa_changes_from_root = {}
    dict_aa_changes_saltational_from_root = {}

    for node in tqdm.tqdm(tree_expansion_from_root):
        if node.parent == None:
            dict_aa_changes_from_root[node.id] = 0
            dict_aa_changes_saltational_from_root[node.id] = 0
        elif not node.id in translations:
            dict_aa_changes_from_root[node.id] = dict_aa_changes_from_root[node.parent.id]
            dict_aa_changes_saltational_from_root[node.id] = dict_aa_changes_saltational_from_root[node.parent.id]
        else:
            n_aa_changes_on_branch = sum(1 for trans in translations[node.id] if is_non_synonymous(trans.aa))
            dict_aa_changes_from_root[node.id] = n_aa_changes_on_branch + dict_aa_changes_from_root[node.parent.id]
    
            if node.id in vec_name_saltational_branches:
                dict_aa_changes_saltational_from_root[node.id] = n_aa_changes_on_branch + dict_aa_changes_saltational_from_root[node.parent.id]
            else:
                dict_aa_changes_saltational_from_root[node.id] = dict_aa_changes_saltational_from_root[node.parent.id]

    df_cum_aa_changes = pd.DataFrame({
        "cum_aa_changes_from_root": dict_aa_changes_from_root,
        "cum_saltational_aa_changes_from_root": dict_aa_changes_saltational_from_root,
    }).rename_axis("node_id").reset_index()

    return(df_cum_aa_changes)

if __name__ == "__main__":

    # Load arguments
    parser = argparse.ArgumentParser(
        description="Get parent / child node relationships",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help='UShER tree in protobuf format')
    parser.add_argument('--name_saltational_branches', type=str, required=True, help='File containing names of saltational branches')
    parser.add_argument('--reference_gtf', type=str, required=True, help='Reference GTF file with gene coordinates')
    parser.add_argument('--reference_fasta', type=str, required=True, help='Reference FASTA file')
    parser.add_argument('--output_cum_aa_changes', type=str, required=True, help='TSV output file with cumulative aa changes from root')
    args = parser.parse_args()

    # Read tree
    tree = read_tree(args.tree)

    # Read name of saltational branches
    vec_name_saltational_branches = read_name_saltational_branches(args.name_saltational_branches)

    # Get AA muts across the tree
    translations = tree.translate(args.reference_gtf, args.reference_fasta)

    # Tree-traversal starting from the root
    tree_expansion_from_root = tree.breadth_first_expansion()

    df_cum_aa_changes = get_cum_aa_changes_from_root(tree_expansion_from_root, translations, vec_name_saltational_branches)

    # Write results to output file
    df_cum_aa_changes.to_csv(args.output_cum_aa_changes, sep = "\t", index = False)
