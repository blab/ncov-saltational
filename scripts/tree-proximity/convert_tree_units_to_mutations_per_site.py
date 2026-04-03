import bte
import argparse
import sys
from collections import defaultdict, Counter

def read_tree(tree_file):
    tree = bte.MATree(tree_file)
    return tree

def convert_tree_units_to_mutations_per_site(tree, genome_length):
    # Get all the nodes in the tree (starting from the root)
    tree_expansion = tree.breadth_first_expansion(reverse = False)

    # Update tree branch length
    for node in tree_expansion:
        node.set_branch_length((node.branch_length )/genome_length)

    print('Finished updating tree branch lengths.')


if __name__ == "__main__":

    # Load arguments
    parser = argparse.ArgumentParser(
        description="Infer location for internal nodes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help='UShER tree in protobuf format')
    parser.add_argument('--genome_length', type=int, required=True, help='Genome length for converting branch lengths')
    parser.add_argument('--output_newick', type=str, required=True, help='Output UShER tree in newick format with updated branch lengths')
    args = parser.parse_args()


    tree = read_tree(args.tree)
    convert_tree_units_to_mutations_per_site(tree, args.genome_length)

    # Write results to output file
    tree.write_newick(args.output_newick)
