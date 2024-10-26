import os
import re
import argparse
import bte
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Process bootstrap spectra output files.")
    parser.add_argument('--spectrum_file', type=str, default="spectra_output.tsv", help="Spectrum file computed from the tree.")
    parser.add_argument("--input_tree", type=str, default="public-latest.all.masked.pb.gz", help="Input tree file (protobuf format)")
    parser.add_argument("--annotate_nodes_output_file", type=str, default="annotated_nodes_output.tsv", help="File to save annotated node data.")
    parser.add_argument("--metadata_output", type=str, default="metadata_output.tsv", help="File to save post-processed metadata.")
    return parser.parse_args()

def import_tsv_to_dict(tsv_file):
    data = {}
    with open(tsv_file, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            values = line.strip().split('\t')
            # because bootstraps won't have representative nodes, skip that.
            row_dict = {header[i]: values[i] for i in range(len(header)-1)}
            data[row_dict[header[0]]] = row_dict
    return data

def get_spectrum_roots(tree, subset_nodes_set):
    annotations = {}

    def traverse_and_annotate(node, current_ancestor):
        # If the node is in the subset, update the current ancestor
        if node.id in subset_nodes_set:
            current_ancestor = node.id
        # Annotate the node and its tips with the current ancestor
        annotations[node.id] = current_ancestor

        # Recursively traverse children
        for child in node.children:
            traverse_and_annotate(child, current_ancestor)

    # Start the traversal from the root
    traverse_and_annotate(tree.root, None)

    return annotations

if __name__ == "__main__":
    # Take input data
    args = parse_args()
    spectra_data = import_tsv_to_dict(args.spectrum_file)
    splits = set(spectra_data.keys())
    tree = bte.MATree(args.input_tree)

    # Get the parent in the splits
    spectrum_roots = get_spectrum_roots(tree, splits)

    # Create a DataFrame to store the results
    result_data = {
        "NodeID": [],
        "SepctrumRoot": [],
        "AC": [], "AG": [], "AT": [],
        "CA": [], "CG": [], "CT": [],
        "GA": [], "GC": [], "GT": [],
        "TA": [], "TC": [], "TG": []
    }

    # Populate the data
    for node, s in spectrum_roots.items():
        result_data["NodeID"].append(node)
        result_data["SepctrumRoot"].append(s)
        result_data["AC"].append(spectra_data[s]["AC"])
        result_data["AG"].append(spectra_data[s]["AG"])
        result_data["AT"].append(spectra_data[s]["AT"])
        result_data["CA"].append(spectra_data[s]["CA"])
        result_data["CG"].append(spectra_data[s]["CG"])
        result_data["CT"].append(spectra_data[s]["CT"])
        result_data["GA"].append(spectra_data[s]["GA"])
        result_data["GC"].append(spectra_data[s]["GC"])
        result_data["GT"].append(spectra_data[s]["GT"])
        result_data["TA"].append(spectra_data[s]["TA"])
        result_data["TC"].append(spectra_data[s]["TC"])
        result_data["TG"].append(spectra_data[s]["TG"])

    # Convert result data into a pandas DataFrame
    df = pd.DataFrame(result_data)

    # Save the annotated node data
    df.to_csv(args.annotate_nodes_output_file, sep='\t', index=False)

    # Post-process to rename 'NodeID' to 'strain' and save the output
    df.rename(columns={"NodeID": 'strain'}, inplace=True)  # rename the 'NodeID' column to 'strain'
    df.to_csv(args.metadata_output, sep='\t', index=False)  # save it as a .tsv extension
