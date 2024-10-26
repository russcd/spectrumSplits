# Miscellaneous scripts
## PCA
```
usage: PCA.py [-h] -i INPUT -o OUTPUT [-r RESULTS]
PCA.py: error: the following arguments are required: -i/--input, -o/--output
```
## Post-process bootstraps
```
usage: process_bootstraps.py [-h] [--bootstrap_directory BOOTSTRAP_DIRECTORY] [--spectrum_file SPECTRUM_FILE] [--input_tree INPUT_TREE]

Process bootstrap spectra output files.

options:
  -h, --help            show this help message and exit
  --bootstrap_directory BOOTSTRAP_DIRECTORY
                        Directory containing the bootstrap spectra output files.
  --spectrum_file SPECTRUM_FILE
                        Spectrum file computed from the tree.
  --input_tree INPUT_TREE
                        Input tree file (protobuf format)
```
## Annotate nodes
```
usage: annotate_nodes.py [-h] [--spectrum_file SPECTRUM_FILE] [--input_tree INPUT_TREE] [--annotate_nodes_output_file ANNOTATE_NODES_OUTPUT_FILE] [--metadata_output METADATA_OUTPUT]

Annotate nodes with spectrum splits.

options:
  -h, --help            show this help message and exit
  --spectrum_file SPECTRUM_FILE
                        Spectrum file computed from the tree.
  --input_tree INPUT_TREE
                        Input tree file (protobuf format)
  --annotate_nodes_output_file ANNOTATE_NODES_OUTPUT_FILE
                        File to save annotated node data.
  --metadata_output METADATA_OUTPUT
                        File to save post-processed metadata.
```
