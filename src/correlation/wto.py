"""
wto.py
======

Module for computing the Weighted Topological Overlap (wTO) matrix from a given
gene expression array. The wTO quantifies the similarity of interaction profiles
between genes in an adjacency matrix derived from correlation measures.

This implementation supports:
- Distance correlation ("dc")
- Pearson correlation ("pearson")

The wTO metric is widely used in gene co-expression network analysis and serves
as a biologically meaningful adjacency refinement before downstream topology-
based analyses such as persistent homology or simplicial complex construction.

References
----------
Zhang, B. and Horvath, S. (2005).
A General Framework for Weighted Gene Co-Expression Network Analysis (WGCNA).
Statistical Applications in Genetics and Molecular Biology, 4(1).
"""

import numpy as np
from src.correlation.distance_corr import compute_distance_correlation_matrix
from src.correlation.pearson import compute_pearson_correlation_matrix

def compute_wto_matrix(df: np.ndarray, adj_matrix: str ="dc")-> np.ndarray:
    """
    Compute the Weighted Topological Overlap (wTO) matrix.

    Parameters:
    - DF: np.ndarray of shape (n_samples, n_genes)
        Raw gene expression matrix with samples as rows and genes as columns.

    - adj_matrix : {"dc", "pearsons"}, default="dc"
        Specifies which adjacency (correlation) matrix to compute:
        
        - "dc"        : Distance correlation
        - "pearsons"  : Pearson correlation

     Returns
    -------
    wto_matrix : np.ndarray of shape (n_genes, n_genes)
        Symmetric weighted topological overlap matrix where each entry wTO(i, j)
        indicates the neighborhood similarity between gene i and gene j.
    """

    if not isinstance(df, np.ndarray):
        raise TypeError("Input 'df' must be a NumPy ndarray.")
    
    if adj_matrix == "dc":
        adjacency_matrix = compute_distance_correlation_matrix(df=df)
    elif adj_matrix == "pearson":
        adjacency_matrix = compute_pearson_correlation_matrix(df=df)
    else:
        raise ValueError("Invalid adjacency matrix type. Choose 'dc' or 'pearson.")

    num_genes = adjacency_matrix.shape[0]
    wto_matrix = np.zeros((num_genes, num_genes))

    for i in range(num_genes):
        for j in range(num_genes):
            if i != j:
                min_ki_kj = (
                    min(
                        np.sum(np.abs(adjacency_matrix[i, :])),
                        np.sum(np.abs(adjacency_matrix[j, :])),
                    )
                    + 1
                    - np.abs(adjacency_matrix[i, j])
                )
                wto_matrix[i, j] = (
                    np.dot(adjacency_matrix[i, :], adjacency_matrix[:, j])
                    + adjacency_matrix[i, j]
                ) / min_ki_kj

    return wto_matrix