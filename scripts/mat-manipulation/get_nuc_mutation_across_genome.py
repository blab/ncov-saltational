import bte
import argparse
import sys
from collections import defaultdict, Counter
import pandas as pd
import tqdm

def read_tree(tree_file):
    tree = bte.MATree(tree_file)
    return tree

def get_nuc_mutation_info(tree_expansion, node_id):
    if tree_expansion[node_id].branch_length > 0:
        curr_mutations = tree_expansion[node_id].mutations
        df_mutations = pd.DataFrame({
            "strain": [tree_expansion[node_id].id] * len(curr_mutations),
            "nuc_mutation": curr_mutations
        })
    else:
        df_mutations = pd.DataFrame(columns=["strain", "nuc_mutation"])
    return df_mutations

def get_nuc_mutation_across_tree(tree_expansion):
    rows = []
    for node_id in range(len(tree_expansion)):
        df_mutations = get_nuc_mutation_info(tree_expansion, node_id)
        if not df_mutations.empty:
            rows.extend(zip(df_mutations["strain"].array, df_mutations["nuc_mutation"].array))

    df_mutations_all = pd.DataFrame.from_records(rows, columns = ["strain", "nuc_mutation"])
    return df_mutations_all

if __name__ == "__main__":

    # Load arguments
    parser = argparse.ArgumentParser(
        description="Get parent / child node relationships",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help='UShER tree in protobuf format')
    parser.add_argument('--output_df_nuc_mutations', type=str, required=True, help='TSV output file with nuc mutation by branch')
    args = parser.parse_args()


    tree = read_tree(args.tree)

    tree_expansion = tree.breadth_first_expansion(reverse = True)

    df_mutations = get_nuc_mutation_across_tree(tree_expansion)

    # Write results to output file
    df_mutations.to_csv(args.output_df_nuc_mutations, sep = "\t", index = False)
