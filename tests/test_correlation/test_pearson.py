import numpy as np
from src.correlation.pearson import compute_pearson_correlation_matrix
from test_correlation.utils import sample_data

def test_pearson_shape():
    df = sample_data()
    dc = compute_pearson_correlation_matrix(df)
    assert dc.shape == (df.shape[1], df.shape[1])


def test_distance_corr_symmetry():
    df = sample_data()
    dc = compute_pearson_correlation_matrix(df)
    assert np.allclose(dc, dc.T)
