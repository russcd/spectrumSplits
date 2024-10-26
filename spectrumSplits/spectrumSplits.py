import bte
import sys
import argparse
from scipy.stats import chi2_contingency
from collections import defaultdict
from multiprocessing import Process, Manager
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Process a phylogenetic tree to find splits, compute spectra, and get representative tips.")
    parser.add_argument("--input_tree", type=str, default="public-latest.all.masked.pb.gz", help="Input tree file (protobuf format)")
    parser.add_argument("--output_tree", type=str, default="masked_sites.pb.gz", help="Output tree file (protobuf format)")
    parser.add_argument("--min_total", type=int, default=500, help="Minimum mutation count to accept a split")
    parser.add_argument("--min_count", type=int, default=50, help="Minimum mutation count to accept a split")
    parser.add_argument("--nthreads", type=int, default=100, help="Number of concurrent threads for processing")
    parser.add_argument("--mask_chi", type=float, default=5000, help="Minimum chi2 value for masking mutation below a node (defaults to off)")
    return parser.parse_args()

def get_position_from_mutation(mutation):
    """
    Extracts the position from a mutation string.
    Mutation strings are of the form \D\d+\D, where \d+ is the numeric position.
    """
    match = re.search(r'\d+', mutation)
    if match:
        return int(match.group(0))  # Return the numeric part as an integer
    return None

def get_mutation_counts(node):
    mutation_counts = defaultdict(int)

    def traverse_tree(node):
        for mutation in node.mutations:
            position = get_position_from_mutation(mutation)
            if position is not None:
                mutation_counts[position] += 1
        for child in node.children:
            traverse_tree(child)

    traverse_tree(node)
    return dict(mutation_counts)

def process_mutation(tree, position, count, total_mutations, args, mask_dict, chi_list):
    """Process one mutation split."""
    find_site_splits(position, count, total_mutations, tree.root, args, mask_dict, chi_list)

def run_in_process(tree, position, count, total_mutations, args, mask_dict, chi_list):
    """Helper function to run in a separate process."""
    p = Process(target=process_mutation, args=(tree, position, count, total_mutations, args, mask_dict, chi_list))
    p.start()
    return p   

def find_site_splits(position, mutation_count, total_mutations, root, args, mask_dict, chi_list):

    mutation_memo = {}
    total_memo = {}
    max_chi = 0
    max_node = root

    def traverse_and_count(node):
        nonlocal max_chi, max_node  # Declare that we want to modify the outer scope variables so we can do the recursion

        if node.id in mutation_memo:
            return mutation_memo[node.id], total_memo[node.id]

        mutation_occurrences = 0
        total_descendant_mutations = len(node.mutations)

        for mutation in node.mutations:
            if get_position_from_mutation(mutation) == position:
                mutation_occurrences += 1

        for child in node.children:
            child_mutation_count, child_total_count = traverse_and_count(child)
            mutation_occurrences += child_mutation_count
            total_descendant_mutations += child_total_count

        mutation_memo[node.id] = mutation_occurrences
        total_memo[node.id] = total_descendant_mutations

        snps_above = mutation_count - mutation_memo[node.id] 
        total_above = total_mutations - total_memo[node.id] - snps_above

        if total_above > args.min_total and total_memo[node.id] > args.min_total:
            observed = [[total_above, snps_above], [total_memo[node.id]-mutation_memo[node.id], mutation_memo[node.id]]]
            chi2, p, dof, expected = chi2_contingency(observed)

            if chi2 > max_chi:
                max_chi = chi2
                max_node = node

        return mutation_occurrences, total_descendant_mutations

    traverse_and_count(root)
    
    # Add to mask_dict if max_chi exceeds the threshold
    if max_chi > args.mask_chi:
        current_nodes = mask_dict.get(max_node.id, [])  # Get the list if it exists, otherwise get an empty list
        current_nodes.append(position)
        mask_dict[max_node.id] = current_nodes  # Reassign the updated list back to the dictionary

    # Append chi result to chi_list
    chi_list.append((position, max_chi, max_node.id))

def mask_mutations(root, mask_dict):
    """
    Perform a DFS traversal of the tree and for each node in mask_dict,
    remove the listed mutations from all descendant nodes using update_mutations.
    """
    def dfs(node):
        if node.id in mask_dict:
            positions_to_mask = mask_dict[node.id]
            print(f"Masking mutations at positions {positions_to_mask} in descendants of node {node.id}", file = sys.stderr)

            # Recursively mask mutations in descendants
            def remove_mutations(descendant):
                # Create a set of positions of current mutations for the descendant
                current_positions = set(get_position_from_mutation(m) for m in descendant.mutations)
                # Calculate remaining positions after masking
                remaining_positions = current_positions - set(positions_to_mask)
                # Update the descendant's mutations with the remaining ones
                remaining_mutations = [m for m in descendant.mutations if get_position_from_mutation(m) in remaining_positions]
                descendant.update_mutations(remaining_mutations, update_branch_length=True)

                # Continue removing mutations for each child of this descendant
                for child in descendant.children:
                    remove_mutations(child)

            # Start masking from the children of the current node
            for child in node.children:
                remove_mutations(child)

        # Continue DFS to check children nodes
        for child in node.children:
            dfs(child)

    # Start the DFS traversal from the root node
    dfs(root)

def main():
    args = parse_args()
    tree = bte.MATree(args.input_tree)

    # Get mutation counts
    mutation_counts = get_mutation_counts(tree.root)

    print("Counting mutations", file=sys.stderr)
    total_mutations = sum(mutation_counts.values())
    
    # Prune mutations with fewer occurrences than the minimum count
    mutation_counts = {k: v for k, v in mutation_counts.items() if v >= args.min_count}

    ## iteratively check for unbalanced splits
    iteration = 1
    while len(mutation_counts.keys()) > 0 :
        print(f"Finding splits. Iteration no: {iteration}", file=sys.stderr)
    
        # Use Manager to create shared data structures
        manager = Manager()
        mask_dict = manager.dict()  # Shared dictionary for storing mutations and node_ids
        chi_list = manager.list()  # Shared list for storing chi results

        # Limit the number of concurrent processes
        processes = []
        for position, count in mutation_counts.items():

            print(f"\tPosition: {position}\tOccurrences: {count}", file=sys.stderr)
            
            # Start a new process for each mutation
            p = run_in_process(tree, position, count, total_mutations, args, mask_dict, chi_list)
            processes.append(p)
            
            # Ensure we don't exceed the specified number of threads
            if len(processes) >= args.nthreads:
                for proc in processes:
                    proc.join()
                processes = []  # Reset the list after joining

        # Wait for remaining processes to finish
        for proc in processes:
            proc.join()

        # Print mutations with chi-square values exceeding mask_chi
        print(f"Mutations checked: ", file=sys.stderr)
        for position, chi, node_id in sorted(chi_list, key=lambda x: -x[1]):
            print(f"{position}\t{chi}\t{node_id}", file=sys.stderr)
        
        # Mask mutations if mask_chi is greater than 0
        if args.mask_chi > 0 and len( mask_dict.keys() ) > 0 :
            print("Masking mutations: ", len(mask_dict.keys()), file=sys.stderr)
            mask_mutations(tree.root, mask_dict)
        
            print("Recounting mutations", file=sys.stderr)
            mutation_counts = get_mutation_counts(tree.root)

            # Flatten the mutation positions from mask_dict into a single set for faster lookup
            masked_positions = set(pos for mutations in mask_dict.values() for pos in mutations)

            # Filter mutation_counts based on whether they are present in the flattened mask positions
            filtered_counts = {k: v for k, v in mutation_counts.items() if k in masked_positions}

            print("Filtered mutation counts:", filtered_counts, file=sys.stderr)
            print("Sites to recheck: ", len(filtered_counts.keys()), file=sys.stderr)
            mutation_counts = filtered_counts

        else :
            mutation_counts = {}

        iteration += 1

    print("Saving tree to: ", args.output_tree, file=sys.stderr)
    tree.save_pb(args.output_tree)  # Use the BTE library's method to save the modified tree

if __name__ == "__main__":
    main()
