"""
pearson.py
==========

Module for computing Pearson correlation matrix. This module uses
 `scipy.stats.pearsonr` to compute correlation coefficients 
 between each pair of genes

 This module is also used to compute the adjacency matrix of the 
 Weighted Topological Overlap (wTO)
"""
import numpy as np
from scipy.stats import pearsonr

def compute_pearson_correlation_matrix(df: np.ndarray) -> np.ndarray:
    """
    Compute the Pearson correlation matrix.

    Parameters:
    - df: np.ndarray of shape (n_samples, n_genes)
        Raw gene expression matrix with samples as rows and genes as columns.

    Returns:
    - corr_matrix: np.ndarray of shape (n_genes, n_genes)
    Symmetric Pearson correlation matrix where entry (i, j) is the Pearson
        correlation coefficient between gene i and gene j.
    """

    if not isinstance(df, np.ndarray):
        raise TypeError("Input 'df' must be a NumPy ndarray.")

    num_genes = df.shape[1]
    corr_matrix = np.zeros((num_genes, num_genes))

    for i in range(num_genes):
        for j in range(i, num_genes):
            if i == j:
                corr_matrix[i, j] = 1.0  # Perfect correlation with itself
            else:
                corr, _ = pearsonr(df[:, i], df[:, j])
                corr_matrix[i, j] = corr
                corr_matrix[j, i] = corr  # Symmetric matrix


    return corr_matrix