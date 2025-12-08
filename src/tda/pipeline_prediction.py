import numpy as np
from tqdm import tqdm
from .tda_complex import (
    build_rips_complex,
    compute_persistence,
  
)
import src.correlation.patient_correlation_measure as pcm
from gudhi.representations import Landscape

def run_prediction_pipeline(
        patient_data: np.ndarray,
        population_matrix: np.ndarray,
        num_layers=2,
        resolution=50,
        max_dim=2,
):
    """
    Generate persistence landscapes for ML.
    """

    landscape = Landscape(num_landscapes=num_layers, resolution=resolution)

    diag0, diag1, diag2 = [], [], []

    for exp in tqdm(patient_data, desc="Computing diagrams"):
        mat = pcm(exp, population_matrix)
        st = build_rips_complex(mat, max_dim)
        _ = compute_persistence(st, max_dim)

        diag0.append(st.persistence_intervals_in_dimension(0))
        diag1.append(st.persistence_intervals_in_dimension(1))
        if max_dim >= 2:
            diag2.append(st.persistence_intervals_in_dimension(2))


    l0 = landscape.fit_transform(diag0)
    l1 = landscape.fit_transform(diag1)
    if max_dim >= 2:
        l2 = landscape.fit_transform(diag2)

        landscapes = np.column_stack([l0, l1, l2])
        return landscapes
    
    else:
        landscapes = np.column_stack([l0, l1])

    return landscapes