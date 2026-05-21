import bte
import argparse
import sys
from collections import defaultdict, Counter
import pandas as pd
import tqdm

def read_metadata(metadata_file, country_mapping_file, location_col):
    metadata = pd.read_csv(metadata_file, sep='\t', usecols = ['strain', 'country_usher'])

    country_mapping = pd.read_csv(country_mapping_file, usecols = ['country_usher', location_col])
    metadata = metadata.merge(country_mapping, how='left', on=['country_usher'])

    metadata = metadata.rename(columns={location_col: "location"})
    return metadata

def read_tree(tree_file):
    tree = bte.MATree(tree_file)
    return tree

def get_df_nodes_from_tree(tree, metadata):
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
        curr_parent = tree_expansion[node_id].parent
        if curr_parent is None:
            list_parents.append("None")
        else:
            list_parents.append(tree_expansion[node_id].parent.id)
        list_branch_lengths.append(tree_expansion[node_id].branch_length)

    # Create a dataframe with node information
    df_nodes = pd.DataFrame({'node_id': list_tree_expansion_ids, 'strain': list_nodes_ids, 'parent_name': list_parents,
                            'branch_length': list_branch_lengths})
    # Add metadata fields
    df_nodes = df_nodes.merge(metadata, how='left',on=['strain'])

    ## Add whether location is missing
    df_nodes['location_missing'] = pd.isna(df_nodes['location'])
    
    return df_nodes

def infer_location(df_nodes, tree):
    # Get all the nodes in the tree (starting from the tips)
    tree_expansion = tree.breadth_first_expansion(reverse = True)

    # Create a strain to location dictionary 
    strain_to_location = dict(zip(df_nodes['strain'], df_nodes['location']))
    strain_to_location_missing = dict(zip(df_nodes['strain'], df_nodes['location_missing']))

    # Infer location of ancestral nodes (if all the children are in the same location) 
    # This replicates the logic used in Sanderson et al.

    ## location_descendant_counts is a dictionary that will store for each node the number of descendant in each location
    location_descendant_counts = defaultdict(Counter) 
    for node in tqdm.tqdm(tree_expansion):
        if not strain_to_location_missing[node.id]:
            location = strain_to_location[node.id]
            location_descendant_counts[node.id][location] = 1
        else:
            for child in node.children:
                for location, count in location_descendant_counts[child.id].items():
                    location_descendant_counts[node.id][location] += count

    # Add inferred locations to df_nodes DF
    list_inferred_locations = []
    for node in tqdm.tqdm(tree_expansion):
        n_locations = len(location_descendant_counts[node.id])
        if n_locations == 1:
            list_inferred_locations.append(list(location_descendant_counts[node.id].keys())[0])
        elif n_locations == 0:
            list_inferred_locations.append(None)
        else:
            list_inferred_locations.append('Multiple descendant locations')

    return list_inferred_locations

if __name__ == "__main__":

    # Load arguments
    parser = argparse.ArgumentParser(
        description="Infer location for internal nodes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help='UShER tree in protobuf format')
    parser.add_argument('--metadata', type=str, required=True, help='TSV metadata containing location for tips')
    parser.add_argument('--country_mapping', type=str, required=True, help='TSV metadata containing country mapping')
    parser.add_argument('--location_col', type=str, default = 'country', required=True, help='Column name for location in metadata file')
    parser.add_argument('--output_locations', type=str, required=True, help='TSV output file with inferred locations')
    args = parser.parse_args()


    tree = read_tree(args.tree)
    metadata = read_metadata(args.metadata, args.country_mapping, args.location_col)

    df_nodes = get_df_nodes_from_tree(tree, metadata)

    df_nodes["inferred_location"] = infer_location(df_nodes, tree)
    results_df = df_nodes[["strain", "location", "inferred_location"]]

    # Write results to output file
    results_df.to_csv(args.output_locations, sep = "\t", index = False)
