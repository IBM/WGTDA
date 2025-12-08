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


def extract_gene_cycles(st, gene_dict):
    pairs = st.persistence_pairs()
    interactions = []

    for idx, ((dim, (birth, death)), (birth_nodes, death_nodes)) in enumerate(pairs):
        birth_genes = [gene_dict[i] for i in birth_nodes]
        death_genes = [gene_dict[i] for i in death_nodes]

        interactions.append({
            "interaction_id": idx,
            "betti_number": dim,
            "birth": float(birth),
            "death": float(death),
            "lifespan": float(death - birth),
            "birth_nodes": list(birth_nodes),
            "death_nodes": list(death_nodes),
            "birth_genes": birth_genes,
            "death_genes": death_genes,
            "geneset": birth_genes + death_genes
        })

    return pd.DataFrame(interactions)


