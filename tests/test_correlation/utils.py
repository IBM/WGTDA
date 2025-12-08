import numpy as np

def sample_data(n_samples: int = 20, n_genes: int = 5, seed: int = 0) -> np.ndarray:
    """
    Generate a random gene expression matrix for testing.

    Parameters
    ----------
    n_samples : int, default=20
    n_genes : int, default=5
    seed : int, default=0

    Returns
    -------
    df : np.ndarray of shape (n_samples, n_genes)
    """
    rng = np.random.default_rng(seed)
    return rng.standard_normal(size=(n_samples, n_genes))