"""
factory.py
=======================

Factory module for constructing gene-to-gene matrices using different
correlation, similarity, or dissimilarity measures.

This provides a *single entry point* to all matrix types used in WGTDA:
    - Pearson correlation
    - Distance correlation
    - Weighted Topological Overlap (WTO)
    - DTEM-based dissimilarity

Using this factory keeps the WGTDA pipeline clean and enables easy
experimentation or swapping of methods via configuration files.

"""

import numpy as np

from src.correlation.pearson import compute_pearson_correlation_matrix
from src.correlation.distance_corr import compute_distance_correlation_matrix
from src.correlation.wto import compute_wto_matrix
from src.correlation.dtem import compute_dtem_matrix

SUPPORTED_METHODS = {
    "pearson",
    "distance",
    "wto_dc",
    "wto_pearsons",
    "dtem",
}

def compute_gene_to_gene_matrix(
    df: np.ndarray,
    method: str,
    **kwargs,
) -> np.ndarray:
    """
    Compute a gene-to-gene matrix using a specified correlation or similarity method.

    Parameters
    ----------
    df : np.ndarray of shape (n_samples, n_genes)
        Gene expression matrix with samples as rows and genes as columns.

    method : {"pearson", "distance", "wto_dc", "wto_pearsons", "dtem"}
        Which matrix to compute:

        - "pearson"       : Pearson correlation
        - "distance"      : Distance correlation
        - "wto_dc"        : WTO using distance correlation adjacency
        - "wto_pearsons"  : WTO using Pearson correlation adjacency
        - "dtem"          : DTEM-like dissimilarity matrix

    **kwargs :
        Extra arguments passed to individual methods.
        For DTEM:
            - k : int
            - r : float

    Returns
    -------
    matrix : np.ndarray of shape (n_genes, n_genes)
        The gene-to-gene matrix computed by the chosen method.

    Raises
    ------
    ValueError
        If an unsupported method is chosen.

    Examples
    --------
    >>> A = compute_gene_to_gene_matrix(df, method="pearson")

    >>> A = compute_gene_to_gene_matrix(df, method="distance")

    >>> A = compute_gene_to_gene_matrix(df, method="dtem", k=4, r=2.0)
    """

    if method not in SUPPORTED_METHODS:
        raise ValueError(
            f"Unsupported method '{method}'. "
            f"Supported methods: {sorted(SUPPORTED_METHODS)}"
        )

    if method == "pearson":
        return compute_pearson_correlation_matrix(df)

    if method == "distance":
        return compute_distance_correlation_matrix(df)

    if method == "wto_dc":
        return compute_wto_matrix(df, adj_matrix="dc")

    if method == "wto_pearsons":
        return compute_wto_matrix(df, adj_matrix="pearsons")

    if method == "dtem":
        return compute_dtem_matrix(df, **kwargs)

    #####EXAMPLE TO ADD NEW METHOD 
    #if method == "mi":
    #return compute_mutual_information_matrix(df)


    # Should never happen because of SUPPORTED_METHODS check
    raise RuntimeError(f"Unhandled method '{method}'")