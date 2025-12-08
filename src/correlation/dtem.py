import numpy as np

def compute_dtem_matrix(df: np.ndarray, k: int=4, r: float=2.0) -> np.ndarray:
    """
    Compute the Distance to Empirical Measure (DTEM) matrix.

    Parameters:
    - df: np.ndarray of shape (n_samples, n_genes)
        Raw gene expression matrix with samples as rows and genes as columns.

    - k: int, number of nearest neighbors.
    
    - r: float, distance metric (Euclidean distance by default).

    Returns:
    - dtem_matrix: np.ndarray, the DTEM matrix.
    """

    if not isinstance(df, np.ndarray):
        raise TypeError("Input 'df' must be a NumPy ndarray.")

    num_genes = df.shape[1]
    dtem_matrix = np.zeros((num_genes, num_genes))

    for i in range(num_genes):
        for j in range(num_genes):
            if i == j:
                dtem_matrix[i, j] = 0.0  # Perfect distance with itself
            else:
                distances = np.linalg.norm(df[i] - df[j], ord=r)
                dtem_matrix[i, j] = distances.sum() / k
                dtem_matrix[j, i] = dtem_matrix[i, j]  # Symmetric matrix

    return dtem_matrix