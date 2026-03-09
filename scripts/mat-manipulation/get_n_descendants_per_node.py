import bte
import argparse
import sys
from collections import defaultdict, Counter
import pandas as pd
import tqdm

def read_tree(tree_file):
    tree = bte.MATree(tree_file)
    return tree

def get_df_descendant_counts(tree):
    # Get all the nodes in the tree (starting from the tips)
    tree_expansion = tree.breadth_first_expansion(reverse = True)

    # descendant_counts is a dictionary that stored for each node the number of descendants
    descendant_counts = {}
    
    for node in tqdm.tqdm(tree_expansion):
        n_children = len(node.children)
        if n_children > 0:
            descendant_counts[node.id] = sum(descendant_counts[child.id] for child in node.children) + n_children
        else:
            descendant_counts[node.id] = 0
    
    df_descendant_counts = pd.DataFrame(list(descendant_counts.items()), columns=["strain", "n_descendants"])
    
    return df_descendant_counts


def get_df_children(tree):
    # Get all the nodes in the tree (starting from the tips)
    tree_expansion = tree.breadth_first_expansion(reverse = True)

    # descendant_counts is a dictionary that stored for each node the number of descendants
    children_counts = {}
    
    for node in tqdm.tqdm(tree_expansion):
        n_children = len(node.children)
        children_counts[node.id] = n_children
    
    df_children_counts = pd.DataFrame(list(children_counts.items()), columns=["strain", "n_children"])
    
    return df_children_counts

if __name__ == "__main__":

    # Load arguments
    parser = argparse.ArgumentParser(
        description="Get parent / child node relationships",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help='UShER tree in protobuf format')
    parser.add_argument('--output_df_n_descendants', type=str, required=True, help='TSV file for output df')
    parser.add_argument('--output_df_n_children', type=str, required=True, help='TSV file for output df children')
    args = parser.parse_args()


    tree = read_tree(args.tree)

    df_descendants = get_df_descendant_counts(tree)
    df_children = get_df_children(tree)

    # Write results to output file
    df_descendants.to_csv(args.output_df_n_descendants, sep = "\t", index = False)
    df_children.to_csv(args.output_df_n_children, sep = "\t", index = False)