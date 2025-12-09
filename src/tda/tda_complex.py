import numpy as np
import pandas as pd
import gudhi as gd
import matplotlib.pyplot as plt


def build_rips_complex(matrix: np.ndarray, max_dim=2):
    rips = gd.RipsComplex(distance_matrix=matrix, max_edge_length=float("inf"))
    st = rips.create_simplex_tree(max_dimension=max_dim)
    st.collapse_edges()
    st.expansion(max_dim + 1)
    return st


def compute_persistence(st, max_dim=2):
    return st.persistence(persistence_dim_max=max_dim)



def save_persistence_diagram(diagrams, filename="pd.png", title="Persistence Diagram"):
    plt.figure(figsize=(8, 6))
    gd.plot_persistence_diagram(diagrams)
    plt.title(title)
    plt.xlabel("Birth")
    plt.ylabel("Death")
    plt.savefig(filename, dpi=400, bbox_inches="tight")
    plt.close()


def extract_gene_cycles(st, gene_dict, max_dim: int | None = None) -> pd.DataFrame:
    """
    Extract topological interactions from a Gudhi SimplexTree and
    return an `interactions` DataFrame matching the original WGTDA format.

    Columns:
        interaction_id, betti_number, birth, death, lifespan,
        birth_nodes, death_nodes, birth_geneset, death_geneset, geneset
    """

    interactions = pd.DataFrame(
        columns=[
            "interaction_id",
            "betti_number",
            "birth",
            "death",
            "lifespan",
            "birth_nodes",
            "death_nodes",
            "birth_geneset",
            "death_geneset",
            "geneset",
        ]
    )
    # Get (betti_number, (birth, death))
    if max_dim is not None:
        betti_pairs = st.persistence(persistence_dim_max=max_dim)
    else:
        betti_pairs = st.persistence()

    # Get the simplex pairs (birth simplex, death simplex)
    persistence_pairs = st.persistence_pairs()

    interactions_list = []

    for index, (
        (betti_number, (birth, death)),
        (birth_nodes, death_nodes),
    ) in enumerate(zip(betti_pairs, persistence_pairs)):
        birth_genes = [gene_dict[item] for item in birth_nodes]
        death_genes = [gene_dict[item] for item in death_nodes]
        gene_set = birth_genes + death_genes
        lifespan = death - birth
        interaction = [
            index,
            betti_number,
            birth,
            death,
            lifespan,
            birth_nodes,
            death_nodes,
            birth_genes,
            death_genes,
            gene_set,
        ]
        interactions_list.append(interaction)

    interactions = pd.DataFrame(interactions_list, columns=interactions.columns)
    return interactions


