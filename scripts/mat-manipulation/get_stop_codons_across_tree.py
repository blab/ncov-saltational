import bte
import argparse
import sys
import pandas as pd
import tqdm

def read_tree(tree_file):
    tree = bte.MATree(tree_file)
    return tree

def read_df_aa_muts(aa_muts_file):
    df_aa_muts = pd.read_csv(aa_muts_file, sep = "\t")
    return df_aa_muts

def get_df_stop_codons_across_tree(tree, df_aa_muts):
    # Create dictionaries with newly occuring stop codons and corresponding branch names
    df_stop_codons = df_aa_muts[df_aa_muts["aa_new"] == "*"]
    dict_all_stop_codons = (df_stop_codons.groupby("node_id")[["gene", "aa_muts_position"]]
                            .apply(lambda g: g.values.tolist())
                            .to_dict())
    
    # Create dictionaries with the stop codons that are reversed (AA change away from *) and corresponding branch names
    df_reverse_stop_codons = df_aa_muts[df_aa_muts["aa_initial"] == "*"]
    df_reverse_stop_codons = df_reverse_stop_codons[df_reverse_stop_codons["aa_new"] != "*"]
    dict_all_reverse_stop_codons = (df_reverse_stop_codons.groupby("node_id")[["gene", "aa_muts_position"]]
                                    .apply(lambda g: g.values.tolist())
                                    .to_dict())
    
    # Get all the nodes in the tree (starting from the root)
    tree_expansion = tree.breadth_first_expansion(reverse = False)

    # Traverse tree and keep tract of all stop codons that have occured on each branch
    # This includes propagating the stop codons from the parent if they aren't reversed

    # Initiate dictionary
    dict_stop_codon_across_tree = dict()
    dict_stop_codon_across_tree[tree_expansion[0].id] = []

    # Traverse tree
    for id_tree in tqdm.tqdm(range(1,len(tree_expansion))):
        curr_node = tree_expansion[id_tree]
        curr_node_id = curr_node.id
        curr_node_parent_id = curr_node.parent.id


        parent_set_stop_codons = {tuple(pos) for pos in dict_stop_codon_across_tree.get(curr_node_parent_id, [])}
        new_stop_codon_to_add = {(gene, int(pos)) for gene, pos in dict_all_stop_codons.get(curr_node_id, [])}
        new_reverse_stop_codons_to_remove = {(gene, int(pos)) for gene, pos in dict_all_reverse_stop_codons.get(curr_node_id, [])}

        curr_stop_codon_set = (parent_set_stop_codons | new_stop_codon_to_add) - new_reverse_stop_codons_to_remove
        
        dict_stop_codon_across_tree[curr_node_id] = [list(x) for x in curr_stop_codon_set]

    # Convert to dataframe
    rows = [
        {"node_id": nid, "gene": g, "position": int(p)}
        for nid, pairs in dict_stop_codon_across_tree.items()
        for g, p in pairs
    ]
    df_summary_stop_codon_across_tree = pd.DataFrame(rows, columns=["node_id", "gene", "position"])

    return df_summary_stop_codon_across_tree


if __name__ == "__main__":

    # Load arguments
    parser = argparse.ArgumentParser(
        description="Get stop codons across tree",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help='UShER tree in protobuf format')
    parser.add_argument('--df_aa_muts', type=str, required=True, help='Dataframe of amino acid mutations')
    parser.add_argument('--output_df_stop_codons', type=str, required=True, help='Output file for stop codons dataframe')
    

    args = parser.parse_args()


    tree = read_tree(args.tree)
    df_aa_muts = read_df_aa_muts(args.df_aa_muts)

    df_stop_codons = get_df_stop_codons_across_tree = get_df_stop_codons_across_tree(tree, df_aa_muts)

    # Write results to output file
    df_stop_codons.to_csv(args.output_df_stop_codons, sep = "\t", index = False)
