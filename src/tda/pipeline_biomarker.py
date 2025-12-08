from .tda_complex import (
    build_rips_complex,
    compute_persistence,
    save_persistence_diagram,
    extract_gene_cycles,
)

def run_biomarker_discovery(matrix, gene_dict, max_dim=2,
                            save_plot=True, filename="pd.png",
                            title="WGTDA Biomarker Discovery"):

    st = build_rips_complex(matrix, max_dim)
    diagrams = compute_persistence(st, max_dim)

    if save_plot:
        save_persistence_diagram(diagrams, filename, title)

    interactions_df = extract_gene_cycles(st, gene_dict)

    return {
        "diagrams": diagrams,
        "interactions": interactions_df
    }