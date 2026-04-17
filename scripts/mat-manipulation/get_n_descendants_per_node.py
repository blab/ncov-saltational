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

def get_df_terminal_descendant_counts(tree, dict_date_tips):
    # Get all the nodes in the tree (starting from the tips)
    tree_expansion = tree.breadth_first_expansion(reverse = True)

    # descendant_counts is a dictionary that stored for each node the number of descendants
    # earliest_date_descendant is a dictionary that stored for each node the earliest date of descendant tips
    # latest_date_descendant is a dictionary that stored for each node the latest date of descendant tips

    descendant_counts = {}
    earliest_date_descendant = {}
    latest_date_descendant = {}

    
    for node in tqdm.tqdm(tree_expansion):
        n_terminal_children = sum(1 for child in node.children if len(child.children) == 0)
        n_children = len(node.children)
        if n_children > 0:
            descendant_counts[node.id] = sum(descendant_counts[child.id] for child in node.children) + n_terminal_children

            child_earliest_dates = [earliest_date_descendant[child.id] for child in node.children if pd.notna(earliest_date_descendant[child.id])]
            child_latest_dates = [latest_date_descendant[child.id] for child in node.children if pd.notna(latest_date_descendant[child.id])]
            earliest_date_descendant[node.id] = min(child_earliest_dates) if child_earliest_dates else pd.NaT
            latest_date_descendant[node.id] = max(child_latest_dates) if child_latest_dates else pd.NaT

        else:
            descendant_counts[node.id] = 0
            earliest_date_descendant[node.id] = dict_date_tips.get(node.id, pd.NaT)
            latest_date_descendant[node.id] = dict_date_tips.get(node.id, pd.NaT)
    
    df_terminal_descendants = pd.DataFrame({
        "strain": descendant_counts.keys(),
        "n_terminal_descendants": descendant_counts.values(),
         "earliest_date_descendants": earliest_date_descendant.values(),
        "latest_date_descendants": latest_date_descendant.values()})
    
    return df_terminal_descendants


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
    parser.add_argument('--metadata_dates', type=str, required=True, help='TSV file containing date information for each tip')
    parser.add_argument('--output_df_n_descendants', type=str, required=True, help='TSV file for output df descendants')
    args = parser.parse_args()


    tree = read_tree(args.tree)

    metadata_with_dates = pd.read_csv(args.metadata_dates, sep='\t', parse_dates=['date'])
    dict_date_tips = metadata_with_dates.set_index("strain")["date"].to_dict()

    df_descendants = get_df_descendant_counts(tree)
    df_terminal_descendants = get_df_terminal_descendant_counts(tree, dict_date_tips)
    df_children = get_df_children(tree)

    # Merge the dataframes
    df_descendants_merge = df_descendants.merge(
        df_terminal_descendants, on="strain", how="left").merge(
            df_children, on="strain", how="left"
        )
    
    # Write results to output file
    df_descendants_merge.to_csv(args.output_df_n_descendants, sep = "\t", index = False)