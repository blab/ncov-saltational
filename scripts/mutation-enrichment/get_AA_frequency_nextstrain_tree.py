import argparse
import sys
import tqdm
import pandas as pd
import json
from Bio import SeqIO
from augur.utils import json_to_tree,annotate_parents_for_tree
from collections import defaultdict, Counter

def read_json_nextstrain_tree(tree_path):
    with open(tree_path, 'r') as f:
        tree_json = json.load(f)

    # Convert tree in Bio.phylo format
    tree = json_to_tree(tree_json)
    tree = annotate_parents_for_tree(tree)

    return tree

def get_df_aa_muts(aa_muts_file, saltational_set):
    df_aa_muts = pd.read_csv(aa_muts_file, sep = "\t")

    df_aa_positions_mutated_in_salt_branches = df_aa_muts[
        (df_aa_muts["is_coded"]) &
        (df_aa_muts["aa_initial"] != df_aa_muts["aa_new"]) &
        (df_aa_muts["node_id"].isin(saltational_set))
    ][["gene", "aa_muts_position"]].drop_duplicates().copy()

    return df_aa_positions_mutated_in_salt_branches

def get_set_saltational_branches(salt_branches_file):
    strains = []

    with open(salt_branches_file, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0]:
                strains.append(parts[0])
    return set(strains)

def get_clade_defining_nodes(tree):
    clade_roots = {}

    for node in tree.find_clades(order="preorder"):
        clade = node.node_attrs.get("clade_membership", {}).get("value")
        parent = getattr(node, "parent", None)
        parent_clade = None if parent is None else parent.node_attrs.get("clade_membership", {}).get("value")

        if clade is not None and clade != parent_clade:
            clade_roots[clade] = node

    return clade_roots

def get_mutations_leading_to_clade(tree, clade_defining_nodes):
    mutations_leading_to_clade_dict = {}
    for clade, node in clade_defining_nodes.items():
        mutations_leading_to_clade_dict[clade] = node.branch_attrs.get("mutations", {})

    return mutations_leading_to_clade_dict

def get_clades_with_defining_mut_at_pos(clade_defining_mutations, gene_of_interest, position_of_interest):
    dict_clade_mutation_on_branch = {}

    for clade, aa_muts in clade_defining_mutations.items():
        aa_muts_on_gene = aa_muts.get(gene_of_interest)
        if aa_muts_on_gene is not None:
            aa_mut_of_interest = [aa_mut for aa_mut in aa_muts_on_gene if aa_mut[1:-1] == str(position_of_interest)]
            if aa_mut_of_interest:
                dict_clade_mutation_on_branch[clade] = aa_mut_of_interest
    return dict_clade_mutation_on_branch

def build_clade_dict(tree):
    nextstrain_clade_dict = {}
    for node in tree.find_clades():
        nextstrain_clade_dict[node.name] = node.node_attrs['clade_membership']['value']
    return nextstrain_clade_dict

def get_clade_sizes(clade_dict):
    return Counter(clade_dict.values())

def build_aa_by_pos_dict(tree, gene, pos, root_value="root_aa"):
    AA_by_node_dict = {}

    for node in tree.find_clades(order="preorder"):
        # Get value from parent
        aa = root_value if node.parent is None else AA_by_node_dict[node.parent.name]

        # Overwrite if a mutation occurred on the branch at this site
        muts = node.branch_attrs.get("mutations", {}).get(gene, [])
        for mut in muts:
            # e.g. "A1000L"
            if len(mut) >= 3 and mut[1:-1] == str(pos):
                aa = mut[-1]

        AA_by_node_dict[node.name] = aa

    return AA_by_node_dict

def get_aa_count_by_clade(tree, clade_dict, gene, pos):
    aa_by_node = build_aa_by_pos_dict(tree, gene, pos, root_value = "root_aa")
    counts = defaultdict(Counter)

    for node_name, aa in aa_by_node.items():
        clade = clade_dict.get(node_name)
        if clade is not None:
            counts[clade][aa] += 1

    return counts

def get_within_clade_mutations(tree, clade_dict, clade_defining_nodes, gene, pos):
    """
    For a given (gene, pos), find clades where the mutant AA arose on a
    non-defining branch inside the clade (i.e. after the clade split).
    Returns {clade: set_of_AAs_that_arose_within_clade}.
    """
    # Build a set of clade-defining node names for fast lookup
    clade_defining_node_names = {node.name for node in clade_defining_nodes.values()}

    within_clade_muts = defaultdict(set)

    for node in tree.find_clades(order="preorder"):
        if node.parent is None:
            continue

        clade = clade_dict.get(node.name)
        parent_clade = clade_dict.get(node.parent.name)

        # Only consider nodes that are inside a clade (not at the clade root)
        if clade is None or clade != parent_clade:
            continue

        # Skip clade-defining root nodes (already captured by is_lineage_defining)
        if node.name in clade_defining_node_names:
            continue

        # Check if this node carries a mutation at (gene, pos)
        muts = node.branch_attrs.get("mutations", {}).get(gene, [])
        for mut in muts:
            if len(mut) >= 3 and mut[1:-1] == str(pos):
                new_aa = mut[-1]
                within_clade_muts[clade].add(new_aa)

    return within_clade_muts

def build_combined_df(tree, clade_dict, clade_defining_mutations, clade_defining_nodes, df_positions):
    """
    For each (gene, pos) in df_positions, compute per-clade AA counts,
    clade sizes, and annotate whether each (clade, AA) is lineage-defining
    or arose within the clade.
    Returns a single flat DataFrame.
    """
    clade_sizes = get_clade_sizes(clade_dict)
    rows_all = []

    for _, row in tqdm.tqdm(df_positions.iterrows(), total=len(df_positions)):
        gene = row["gene"]
        pos = row["aa_muts_position"]

        # AA counts per clade
        counts_dict = get_aa_count_by_clade(tree, clade_dict, gene, pos)

        # Lineage-defining: {clade: [mut_string, ...]}
        lineage_def = get_clades_with_defining_mut_at_pos(clade_defining_mutations, gene, pos)
        # Map clade -> set of AAs that are lineage-defining at this pos
        lineage_def_aas = {
            clade: {mut_str[-1] for mut_str in mut_list}
            for clade, mut_list in lineage_def.items()
        }

        # Within-clade mutations: {clade: set of AAs}
        within_clade = get_within_clade_mutations(tree, clade_dict, clade_defining_nodes, gene, pos)

        for clade, aa_counts in counts_dict.items():
            for aa, count in aa_counts.items():
                rows_all.append({
                    "gene": gene,
                    "pos": pos,
                    "AA": aa,
                    "clade": clade,
                    "count": count,
                    "clade_size": clade_sizes[clade],
                    "is_lineage_defining": aa in lineage_def_aas.get(clade, set()),
                    "is_within_clade_mut": aa in within_clade.get(clade, set()),
                })

    return pd.DataFrame(rows_all, columns=["gene", "pos", "AA", "clade", "count", "clade_size", "is_lineage_defining", "is_within_clade_mut"])


if __name__ == "__main__":

    # Load arguments
    parser = argparse.ArgumentParser(
        description="Get AA frequency by clade from Nextstrain tree, with lineage-defining and within-clade mutation flags",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--nextstrain_global_tree', type=str, required=True, help='Nextstrain global tree JSON file')
    parser.add_argument('--df_aa_muts', type=str, required=True, help='TSV file with amino acid mutations')
    parser.add_argument('--saltational_branches', type=str, required=True, help='TSV file with saltational branches names')

    parser.add_argument('--output_df_combined', type=str, required=True, help='Output TSV with AA counts, clade sizes, and mutation origin flags per (gene, pos, AA, clade)')

    args = parser.parse_args()

    # Get nextstrain tree
    nextstrain_tree = read_json_nextstrain_tree(args.nextstrain_global_tree)

    # Get saltational branches
    saltational_set = get_set_saltational_branches(args.saltational_branches)

    # Get AA mutations positions found in saltational branches
    df_aa_positions_mutated_in_salt_branches = get_df_aa_muts(args.df_aa_muts, saltational_set)

    # Get clade-defining nodes and their mutations
    clade_defining_nodes = get_clade_defining_nodes(nextstrain_tree)
    clade_defining_mutations = get_mutations_leading_to_clade(nextstrain_tree, clade_defining_nodes)

    # Get clade associated with each node in the Nextstrain tree
    clade_all_nodes_dict = build_clade_dict(nextstrain_tree)

    # Build combined output
    df_combined = build_combined_df(
        nextstrain_tree,
        clade_all_nodes_dict,
        clade_defining_mutations,
        clade_defining_nodes,
        df_aa_positions_mutated_in_salt_branches
    )

    # Write results to output file
    df_combined.to_csv(args.output_df_combined, sep="\t", index=False)
