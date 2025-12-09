import numpy as np
from src.correlation.dtem import compute_dtem_matrix
from test_correlation.utils import sample_data


def test_dtem_shape():
    df = sample_data()
    dtem = compute_dtem_matrix(df, k=4, r=2.0)
    assert dtem.shape == (df.shape[1], df.shape[1])


def test_dtem_symmetry():
    df = sample_data()
    dtem = compute_dtem_matrix(df, k=4, r=2.0)
    assert np.allclose(dtem, dtem.T)


def test_dtem_diagonal_zero():
    df = sample_data()
    dtem = compute_dtem_matrix(df, k=4, r=2.0)
    assert np.allclose(np.diag(dtem), 0.0)
