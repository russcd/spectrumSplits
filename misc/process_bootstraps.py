import os
import re
import argparse
import bte
import pandas as pd
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Process bootstrap spectra output files.")
    parser.add_argument('--bootstrap_directory', type=str, default="./", help="Directory containing the bootstrap spectra output files.")
    parser.add_argument('--spectrum_file', type=str, default="spectra_output.tsv", help="Spectrum file computed from the tree.")
    parser.add_argument("--input_tree", type=str, default="public-2024-08-06.masked.pb.gz", help="Input tree file (protobuf format)")
    return parser.parse_args()

def import_tsv_to_dict(tsv_file):
    """Read TSV file and store data in a dictionary."""
    data = {}
    with open(tsv_file, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            values = line.strip().split('\t')
            row_dict = {header[i]: values[i] for i in range(len(header)-1)}
            data[row_dict[header[0]]] = row_dict
    return data

def process_bootstrap_files(bootstrap_dir):
    """Process bootstrap files from the specified directory."""
    regex = re.compile(r'bootstrap_(\d+)_splits_output.tsv')
    spectra_dict = {}
    for filename in os.listdir(bootstrap_dir):
        if regex.match(filename):
            file_path = os.path.join(bootstrap_dir, filename)
            spectra_dict[filename] = import_tsv_to_dict(file_path)
    return spectra_dict

def bootstrapSplits(bootstrap_data, spectra_data):
    """Calculate support rates for each node in the bootstrap data."""
    rates = {}
    for replicate in bootstrap_data.keys():
        for node in bootstrap_data[replicate].keys():
            if node not in rates:
                rates[node] = 0
            rates[node] += 1
    recovery = {node: rates.get(node, 0) for node in spectra_data}
    return recovery

def nodeDistance(tree, node1, node2):
    """Calculate the distance between two nodes in the tree."""
    possible_lcas_order = [node1]
    if node1 != "node_1":
        possible_lcas_order.extend([anc.id for anc in tree.rsearch(node1)])
    possible_lcas = set(possible_lcas_order)

    new_ancestors_order = [node2]
    if node2 != "node_1":
        new_ancestors_order.extend([anc.id for anc in tree.rsearch(node2)])
    new_ancestors = set(new_ancestors_order)

    possible_lcas = possible_lcas.intersection(new_ancestors)
    LCA = next(pl for pl in possible_lcas_order if pl in possible_lcas)
    
    return possible_lcas_order.index(LCA) + new_ancestors_order.index(LCA)

def getDistances(tree, splits, bootstrap_spectra_data):
    """Get distance to the nearest split for each node."""
    distanceDict = {}
    for node in splits:
        distances = []
        for replicate in bootstrap_spectra_data.keys():
            minDist = 1000000
            if node in bootstrap_spectra_data[replicate].keys():
                minDist = 0
            else:
                for replicateNode in bootstrap_spectra_data[replicate].keys():
                    dist = nodeDistance(tree, node, replicateNode)
                    minDist = min(minDist, dist)
            distances.append(minDist)
        distanceDict[node] = distances
    return distanceDict

def get_spectrum_roots(tree, subset_nodes_set):
    """Get the spectrum roots by annotating the tree nodes."""
    annotations = {}
    def traverse_and_annotate(node, current_ancestor):
        if node.id in subset_nodes_set:
            current_ancestor = node.id
        if current_ancestor not in annotations:
            annotations[current_ancestor] = []
        annotations[current_ancestor].append(node.id)
        for child in node.children:
            traverse_and_annotate(child, current_ancestor)
    traverse_and_annotate(tree.root, None)
    return annotations

def jaccard_similarity(set1, set2):
    """Calculate Jaccard similarity between two sets."""
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    if not union:
        return 0.0
    return len(intersection) / len(union)

### consider minhashing this
def max_jaccard_similarity(dict1, dict2, result):
    """Find the maximum Jaccard similarity for each key in dict1 against all keys in dict2."""
    for key1, set1 in dict1.items():
        max_similarity = 0.0
        best_match = None
        for key2, set2 in dict2.items():
            similarity = jaccard_similarity(set1, set2)
            if similarity > max_similarity:
                max_similarity = similarity
                best_match = key2
        if key1 not in result :
            result[key1] = [max_similarity]
        else :
            result[key1].append( max_similarity )
    return result

if __name__ == "__main__":
    args = parse_args()
    bootstrap_spectra_data = process_bootstrap_files(args.bootstrap_directory)
    spectra_data = import_tsv_to_dict(args.spectrum_file)
    tree = bte.MATree(args.input_tree)

    # Calculate support probabilities and distances
    supportProps = bootstrapSplits(bootstrap_spectra_data, spectra_data)
    nearestSplitDistance = getDistances(tree, supportProps.keys(), bootstrap_spectra_data)

    # Get spectrum roots as sets for faster Jaccard similarity calculations
    bootstrap_jaccard = {}
    spectrum_tips = {key: set(values) for key, values in get_spectrum_roots(tree, spectra_data.keys()).items()}
    for bootstrap in bootstrap_spectra_data.keys():
        print( "computing jaccard for bootstrap: ", bootstrap, file=sys.stderr )
        bootstrap_spectrum_tips = {key: set(values) for key, values in get_spectrum_roots(tree, bootstrap_spectra_data[bootstrap].keys()).items()}
        bootstra_jaccard = max_jaccard_similarity(spectrum_tips, bootstrap_spectrum_tips, bootstrap_jaccard) 

    for node, prob in supportProps.items():
        print(node, prob / len(bootstrap_spectra_data.keys()), sum(nearestSplitDistance[node]) / len(nearestSplitDistance[node]), sum(bootstrap_jaccard[node])/len(bootstrap_jaccard[node]), ','.join(map(str, nearestSplitDistance[node])), ','.join(map(str,bootstrap_jaccard[node])))
