import dcor
import numpy as np


def compute_distance_correlation_matrix(gene_exp_arr: np.ndarray) -> np.ndarray:
    """
    Compute the distance correlation matrix for a given DataFrame.

    Parameters:
    - df: np.ndarray input data.

    Returns:
    - dist_corr_matrix: np.ndarray, the distance correlation matrix.
    """
    num_genes = gene_exp_arr.shape[1]
    dist_matrix = np.zeros((num_genes, num_genes))

    for i in range(num_genes):
        for j in range(i + 1, num_genes):
            dist_matrix[i, j] = dcor.distance_correlation(
                gene_exp_arr[:, i], gene_exp_arr[:, j]
            )

    # Symmetrize the matrix and set diagonal elements to 1
    dist_matrix = dist_matrix + dist_matrix.T + np.eye(num_genes)

    # Convert to distance measure
    dist_corr_matrix = 1 - dist_matrix

    return dist_corr_matrix


def compute_wto_matrix(gene_exp_arr: np.ndarray)-> np.ndarray:
    """
    Compute the Signed Weighted Topological Overlap (wTO) matrix.

    Parameters:
    - df: np.ndarray, input data.

    Returns:
    - wto_matrix: np.ndarray, Signed Topological Overlap matrix
    """
    adjacency_matrix = compute_distance_correlation_matrix(gene_exp_arr=gene_exp_arr)

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
