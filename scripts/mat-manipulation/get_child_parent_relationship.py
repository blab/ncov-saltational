import bte
import argparse
import sys
from collections import defaultdict, Counter
import pandas as pd
import tqdm

def read_tree(tree_file):
    tree = bte.MATree(tree_file)
    return tree

def get_df_nodes_from_tree(tree):
    # Get all the nodes in the tree (starting from the tips)
    tree_expansion = tree.breadth_first_expansion(reverse = True)
    # Get number of nodes in the tree
    n_nodes = len(tree_expansion)

    list_tree_expansion_ids = []
    list_nodes_ids = []
    list_parents = []
    list_branch_lengths = []
    # Create a dataframe with nodes and their parents
    for node_id in range(n_nodes):
        list_tree_expansion_ids.append(node_id)
        list_nodes_ids.append(tree_expansion[node_id].id)
        list_branch_lengths.append(tree_expansion[node_id].branch_length)
        curr_parent = tree_expansion[node_id].parent
        if curr_parent is None:
            list_parents.append("None")
        else:
            list_parents.append(tree_expansion[node_id].parent.id)
        

    # Create a dataframe with node information
    df_nodes = pd.DataFrame({'node_id': list_tree_expansion_ids, 'strain': list_nodes_ids, 'parent_name': list_parents,
                            'branch_length': list_branch_lengths})

    return df_nodes

if __name__ == "__main__":

    # Load arguments
    parser = argparse.ArgumentParser(
        description="Get parent / child node relationships",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help='UShER tree in protobuf format')
    parser.add_argument('--output_parent_relationship', type=str, required=True, help='TSV output file with relationships')
    args = parser.parse_args()


    tree = read_tree(args.tree)

    df_nodes = get_df_nodes_from_tree(tree)

    # Write results to output file
    df_nodes.to_csv(args.output_parent_relationship, sep = "\t", index = False)
