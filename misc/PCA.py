import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
import argparse

# Set up argparse to handle command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Perform PCA on input spectrum data and save the results.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input spectrum file (TSV format)')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output PDF file for the PCA plot')
    parser.add_argument('-r', '--results', type=str, default='pca_results.tsv', help='Output TSV file for the PCA results')
    return parser.parse_args()

def main():
    # Parse the command-line arguments
    args = parse_args()

    # Load the TSV file into a DataFrame
    input_file = args.input
    data = pd.read_csv(input_file, sep='\t', index_col=0)

    # Select columns 5 to 15 (4 to 14 in Python's zero-based index)
    data_subset = data.iloc[:, 4:15]

    # Perform PCA on the selected columns
    pca = PCA(n_components=10)
    pca_result = pca.fit_transform(data_subset)

    # Convert the result into a DataFrame for easier handling
    pca_df = pd.DataFrame(pca_result, columns=[f'PC{i+1}' for i in range(pca.n_components_)], index=data.index)

    # Save the PCA results to a new TSV file
    pca_df.to_csv(args.results, sep='\t')

    # Get the loadings (contributions of each feature to the principal components)
    loadings = pd.DataFrame(pca.components_.T, columns=[f'PC{i+1}' for i in range(pca.n_components_)], index=data_subset.columns)

    # Explained variance
    explained_variance = pca.explained_variance_ratio_

    # Create the figure and three subplots, smaller size (half)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

    # Scatterplot of the first two principal components
    ax1.scatter(pca_df['PC1'], pca_df['PC2'], c='blue', alpha=0.5)
    ax1.set_xlabel('Principal Component 1 (78%)')
    ax1.set_ylabel('Principal Component 2 (21%)')
    ax1.set_title('')
    ax1.grid(False)  # Remove gridlines

    # Plot the loadings of the first three principal components
    bar_width = 0.25  # Adjusted width for three sets of bars
    index = np.arange(len(loadings))  # Position of the bars on the x-axis

    # Plot PC1, PC2, and PC3 loadings side by side
    ax2.bar(index, loadings['PC1'], bar_width, label='PC1', alpha=0.6)
    ax2.bar(index + bar_width, loadings['PC2'], bar_width, label='PC2', alpha=0.6)

    # Add labels, title, and legend for the loadings plot
    ax2.set_xlabel('Features')
    ax2.set_ylabel('Loadings')
    ax2.set_title('')
    ax2.set_xticks(index + bar_width)  # Adjusted to center labels between the bars
    ax2.set_xticklabels(loadings.index, rotation=45, ha='right')
    ax2.legend()
    ax2.grid(False)  # Remove gridlines

    # Adjust layout and save the plot to a PDF file
    plt.tight_layout()
    plt.savefig(args.output, format='pdf')

    # Optionally, show the plot
    plt.show()

if __name__ == '__main__':
    main()
