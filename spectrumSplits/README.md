# Spectrum splits

```
usage: spectrumSplits.py [-h] [--input_tree INPUT_TREE] [--output_spectrum OUTPUT_SPECTRUM] [--min_chi MIN_CHI] [--min_mutations MIN_MUTATIONS] [--ntips NTIPS]
                         [--bootstrap_splits BOOTSTRAP_SPLITS] [--bootstrap_spectra BOOTSTRAP_SPECTRA] [--nthreads NTHREADS] [--max_branch_length MAX_BRANCH_LENGTH]

Process a phylogenetic tree to find splits, compute spectra, and get representative tips.

options:
  -h, --help            show this help message and exit
  --input_tree INPUT_TREE
                        Input tree file (protobuf format)
  --output_spectrum OUTPUT_SPECTRUM
                        Output TSV file for spectra
  --min_chi MIN_CHI     Minimum Chi-square value to accept a split
  --min_mutations MIN_MUTATIONS
                        Minimum number of mutations required for a split
  --ntips NTIPS         Number of tips to retrieve for each split
  --bootstrap_splits BOOTSTRAP_SPLITS
                        Number of bootstrap replicates to attempt in defining splits
  --bootstrap_spectra BOOTSTRAP_SPECTRA
                        Number of bootstrap replicates to attempt in defining spectra
  --nthreads NTHREADS   Number of threads for concurrent bootstrapping
  --max_branch_length MAX_BRANCH_LENGTH
                        Maximum branch length to include in spectrum calculations
```
