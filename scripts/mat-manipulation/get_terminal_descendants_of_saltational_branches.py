import bte
import argparse
import pandas as pd
import tqdm


def read_tree(tree_file):
    tree = bte.MATree(tree_file)
    return tree


def get_terminal_descendants_of_saltational_branches(tree, saltational_branch_ids):
    # Build a lookup dict from node ID to node object
    node_index = {node.id: node for node in tqdm.tqdm(tree.breadth_first_expansion(), desc="Indexing nodes")}

    rows = []
    for salt_id in tqdm.tqdm(saltational_branch_ids, desc="Processing saltational branches"):
        if salt_id not in node_index:
            print(f"Warning: saltational branch '{salt_id}' not found in tree, skipping.")
            continue
        start_node = node_index[salt_id]
        stack = list(start_node.children)
        while stack:
            node = stack.pop()
            if len(node.children) == 0:
                rows.append((salt_id, node.id))
            else:
                stack.extend(node.children)

    df = pd.DataFrame(rows, columns=["strain_saltational", "strain_descendant"])
    return df


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Get terminal descendants of saltational branches",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help='UShER tree in protobuf format')
    parser.add_argument('--name_saltational_branches', type=str, required=True,
                        help='File with one saltational branch node ID per line')
    parser.add_argument('--output', type=str, required=True, help='TSV output file')
    args = parser.parse_args()

    tree = read_tree(args.tree)

    with open(args.name_saltational_branches) as f:
        saltational_branch_ids = [line.strip() for line in f if line.strip()]

    df = get_terminal_descendants_of_saltational_branches(tree, saltational_branch_ids)

    df.to_csv(args.output, sep="\t", index=False)
