import numpy as np
from src.correlation.wto import compute_wto_matrix
from test_correlation.utils import sample_data


def test_wto_shape_pearson():
    df = sample_data()
    wto = compute_wto_matrix(df, adj_matrix="pearson")
    assert wto.shape == (df.shape[1], df.shape[1])


def test_wto_off_diagonal_not_all_zero_pearson():
    df = sample_data()
    wto = compute_wto_matrix(df, adj_matrix="pearson")
    off_diag = wto[~np.eye(wto.shape[0], dtype=bool)]
    assert np.any(off_diag != 0.0)

def test_wto_shape_dc():
    df = sample_data()
    wto = compute_wto_matrix(df, adj_matrix="dc")
    assert wto.shape == (df.shape[1], df.shape[1])


def test_wto_off_diagonal_not_all_zero_dc():
    df = sample_data()
    wto = compute_wto_matrix(df, adj_matrix="dc")
    off_diag = wto[~np.eye(wto.shape[0], dtype=bool)]
    assert np.any(off_diag != 0.0)
