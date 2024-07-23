
import matilda
import numpy as np
import pandas as pd
from matilda import FilteredSimplicialComplex, PersistentHomologyComputer
from numpy import ndarray
from typing import Tuple


def construct_vr_complex_rna_matrix(
    preprocessing: ndarray, dimensions=3
) -> Tuple[PersistentHomologyComputer, FilteredSimplicialComplex]:
    """
    Construct a Vietoris-Rips complex from the specified preprocessing matrix and compute its persistent homology.

    Parameters:
    - preprocessing : ndarray
        A square matrix where element [i, j] represents the distance between the i-th and j-th elements.
    - dimension : int, default = 3
        The maximum dimension of simplices to be considered in the Vietoris-Rips complex. Default is 3.

    Returns:
    - tuple[PersistentHomologyComputer, FilteredSimplicialComplex]
        A tuple containing the persistent homology computer and the Vietoris-Rips complex.
        The PersistentHomologyComputer object has methods to access and manipulate the computed homology.
        The FilteredSimplicialComplex object represents the simplicial complex constructed from the given distance matrix.

    Notes:
    - The function uses the 'matilda' library, which must be installed and properly configured in the environment.
    - This function is intended for use in computational topology, particularly in the analysis of high-dimensional data.

    Examples:
    - Constructing the complex and computing homology for a given distance matrix:
         matrix = np.array([[0, 1, 1.5], [1, 0, 2], [1.5, 2, 0]])
     persistence, rips_complex = VRcomplex(matrix, dimension=2)
    """
    # Create a FilteredSimplicialComplex object
    rips_complex = matilda.FilteredSimplicialComplex()

    # Construct the Vietoris-Rips complex from the preprocessed distance matrix
    rips_complex.construct_vietoris_from_metric(
        matrix=preprocessing, dimension=dimensions, upper_bound=np.inf
    )

    persistence = matilda.PersistentHomologyComputer()

    # Compute the persistent homology of the VR complex
    persistence.compute_persistent_homology(
        rips_complex, with_representatives=True, verbose=True
    )

    return persistence, rips_complex


def interactions_dataframe(
    persistence: PersistentHomologyComputer,
    rips_complex: FilteredSimplicialComplex,
    gene_dict: dict,
) -> pd.DataFrame:
    """
    Extract and put in the topological features into a dataframe and csv from persistent homology computations in a dataframe.

    Parameters:
    - persistence : PersistentHomologyComputer
        An object containing computed homological features and persistence bars.
    - rips_complex : FilteredSimplicialComplex
        The simplicial complex from which these homological features are computed.
    - gene_dict : dict
        A dictionary mapping vertex indices to gene names.

    Returns:
    - pd.DataFrame
        A DataFrame containing information about each topological interaction, including:
        interaction_id, Betti number, birth and death times of features, lifespan,
        vertices, and the gene sets.

    This function processes the results of persistent homology calculations to extract
    meaningful biological interactions, which can help in understanding the connectivity
    and interaction between different genes within the analyzed dataset.
    """
    interactions_list = []
    interactions = pd.DataFrame(
        columns=[
            "interaction_id",
            "betti_number",
            "birth",
            "death",
            "lifespan",
            "vertices",
            "vertices_set",
        ]
    )

    # Organize persistent cycles and betti pairs for easy access
    betti_pairs = [
        (outer_key, (inner_key, tuple(inner_value)))
        for outer_key, inner_dict in persistence.bars.items()
        for inner_key, inner_value in inner_dict.items()
    ]
    persistence_pairs = {
        outer_key: {
            middle_key: list(inner_dict.keys())
            for middle_key, inner_dict in middle_dict.items()
        }
        for outer_key, middle_dict in persistence.persistent_cycles.items()
    }
    persistence_pairs = [
        (outer_key, (inner_key, value))
        for outer_key, middle_dict in persistence_pairs.items()
        for inner_key, value in middle_dict.items()
    ]

    # Process each pair and populate the dataframe
    for index, (
        (betti_number, (interaction_id, (birth, death))),
        (betti_number2, (interaction_id2, vertices)),
    ) in enumerate(zip(betti_pairs, persistence_pairs)):
        lifespan = death - birth
        vertices_id = [rips_complex.simplices[i] for i in vertices]
        vertices_id = [list(item) for item in vertices_id]
        vertices_set = [
            [gene_dict[item] for item in sublist] for sublist in vertices_id
        ]
        interaction = [
            index,
            betti_number,
            birth,
            death,
            lifespan,
            vertices,
            vertices_set,
        ]
        interactions_list.append(interaction)
    new_interactions = pd.DataFrame(interactions_list, columns=interactions.columns)

    interactions = pd.concat([interactions, new_interactions], ignore_index=True)

    return interactions
