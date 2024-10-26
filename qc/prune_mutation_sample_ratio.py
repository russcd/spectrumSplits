import bte
import sys
import argparse
import numpy as np

# Command-line argument parsing
def parse_args():
    parser = argparse.ArgumentParser(description="Process a phylogenetic tree for changepoint detection in mutation/descendant ratios.")
    parser.add_argument("--input_tree", type=str, default="public-2024-08-06.masked.pb.gz", help="Input tree file (protobuf format)")
    parser.add_argument("--output_tree", type=str, default="pruned_tree.pb.gz", help="Output tree file (protobuf format)")
    parser.add_argument("--threshold", type=float, default=0, help="Mutation:Leaf Ratio to prune")
    return parser.parse_args()

# Compute the mutation-to-descendant ratio for each node and store it
def compute_descendants_mutations_ratio(node, mutation_ratio):
    # Base case: If node is a tip (no children), count it as one descendant
    if not node.children:
        mutation_ratio[node.id] = len(node.mutations) / 1  # 1 because it's a tip itself
        return 1, len(node.mutations)

    total_tips = 0
    total_mutations = len(node.mutations)

    # Recursive case: Traverse children to count tips and mutations
    for child in node.children:
        child_tips, child_mutations = compute_descendants_mutations_ratio(child, mutation_ratio)
        total_tips += child_tips  # Only count tips
        total_mutations += child_mutations

    # Calculate the mutation-to-tip ratio
    if total_tips > 0:
        ratio = total_mutations / total_tips
    else:
        ratio = float('inf')  # Handle case if somehow there are no tips

    # Store the ratio for this node
    mutation_ratio[node.id] = ratio

    return total_tips, total_mutations  # Only return the number of tips

# Traverse the tree and detect changepoints based on mutation/descendant ratio changes
def detect_changepoints(node, mutation_ratio, threshold, changepoints, to_prune):
    for child in node.children:
        if mutation_ratio[child.id] >= threshold:
            changepoints.append((node.id, child.id, mutation_ratio[child.id]))
            # Mark this child for removal
            to_prune.add(child)
            continue
        detect_changepoints(child, mutation_ratio, threshold, changepoints, to_prune)

# Determine the threshold based on the overall distribution of ratios
def compute_threshold(mutation_ratios):
    ratios = np.array(list(mutation_ratios.values()))
    return np.mean(ratios) + np.std(ratios) * 2  ### seems reasonable, all in all, but probably will be a low number in general. 

# Function to prune marked nodes from the tree
def prune_tree(tree, to_prune):
    for node in to_prune :
        tree.remove_node(node.id)

def main():
    args = parse_args()
    tree = bte.MATree(args.input_tree)

    mutation_ratio = {}
    compute_descendants_mutations_ratio(tree.root, mutation_ratio)

    ## compute it from this
    if args.threshold == 0 :
        args.threshold = compute_threshold(mutation_ratio)

    changepoints = []
    to_prune = set()  # Set for nodes to be removed
    detect_changepoints(tree.root, mutation_ratio, args.threshold, changepoints, to_prune)

    print("Changepoints Detected:")
    for parent_id, child_id, ratio in changepoints:
        print(f"Parent {parent_id}, Child {child_id}, Ratio of Child: {ratio}")

    # Prune the marked nodes from the tree
    prune_tree(tree, to_prune)

    # Save the modified tree to the output file
    tree.save_pb(args.output_tree)  # Use the BTE library's method to save the modified tree

if __name__ == "__main__":
    main()
