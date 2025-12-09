"""
distance_corr.py
================

Module for computing a distance correlation matrix. This implementation uses the `dcor` package.

Distance correlation captures both linear and non-linear dependencies between
two random variables, making it a more flexible alternative to Pearson
correlation.
"""
import numpy as np
import dcor

def compute_distance_correlation_matrix(df: np.ndarray) -> np.ndarray:
    """
    Compute the Pearson correlation matrix. 

    Parameters:
    - df: np.ndarray of shape (n_samples, n_genes)
        Raw gene expression matrix with samples as rows and genes as columns.

    Returns:
    - corr_matrix: np.ndarray of shape (n_genes, n_genes)
    Symmetric distance correlation matrix where entry (i, j) is the Pearson
        correlation coefficient between gene i and gene j.
    """
    
    num_genes = df.shape[1]
    dist = np.zeros((num_genes, num_genes))
    
    for i in range(num_genes):
        
        for j in range(i+1, num_genes):
            
            dist[i,j] = dcor.distance_correlation(df[:,i], df[:,j]) #Distance Correlations 
    
    dist_matrix = dist + dist.T + np.eye(num_genes)
    
    return  dist_matrix
