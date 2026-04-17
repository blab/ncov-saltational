import bte
import argparse
import pandas as pd
import tqdm

def read_tree(tree_file):
    tree = bte.MATree(tree_file)
    return tree

def get_clade_per_node(tree):
    # Get all the nodes in the tree (starting from the root)
    tree_expansion = tree.breadth_first_expansion(reverse = False)

    # clade_dict is a dictionary that stored for each node the corresponding clade
    clade_dict = {}
    
    for node in tqdm.tqdm(tree_expansion):
        curr_clade = node.annotations[0]
        if curr_clade == "":
            clade_dict[node.id] = clade_dict[node.parent.id]
        else:
            clade_dict[node.id] = curr_clade
    
    df_clades = pd.DataFrame(list(clade_dict.items()), columns=["strain", "clade"])
    
    return df_clades


if __name__ == "__main__":

    # Load arguments
    parser = argparse.ArgumentParser(
        description="Get clade from viridian tree",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help='UShER tree in protobuf format')
    parser.add_argument('--output_df_clades', type=str, required=True, help='TSV file for output df clades')
    args = parser.parse_args()


    tree = read_tree(args.tree)

    df_clades = get_clade_per_node(tree)

    
    # Write results to output file
    df_clades.to_csv(args.output_df_clades, sep = "\t", index = False)