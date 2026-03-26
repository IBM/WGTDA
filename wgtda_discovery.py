import src.tda.pipeline_biomarker as tda_bio
import argparse
import pandas as pd
import src.correlation.dtem as dtem
import src.correlation.wto as wto
import src.correlation.pearson as pearson
import src.correlation.distance_corr as dc
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="Run the prediction pipeline.")
    parser.add_argument("--gene_expression_matrix", "-p", type=str, required=True,
                        help="Path to the gene expression matrix file.")
    parser.add_argument("--preselection_genes", "-pp", type=str, default=None, required=False,
                        help="Path to the preselection genes file.")
    parser.add_argument("--coexpression", "-c", type=str, default="wto_pearson",
                        choices=["wto_dc", "wto_pearson", "dc", "pearson", "dtem"])
    
    parser.add_argument("--padj_threshold", "-padj", type=float, default=None,
                        help="Adjusted p-value threshold for gene selection."
                )
    parser.add_argument("--log2fc", "-l", type=int, default=None,
                        help="Log2 fold change threshold for gene selection."
                )
    parser.add_argument("--max_dim", type=int, default=2,
                        help="Maximum homology dimension to compute.")
    parser.add_argument("--interactions_output", "-o", default="interactions.csv", type=str,
                        help="Path to save the interactions DataFrame.")
    
    return parser.parse_args() 



def main():
    args = parse_args()
    # Here you would typically load your data using the provided file paths
    # For example:
    expression_matrix_path = args.gene_expression_matrix
    preselection_genes_path = args.preselection_genes
    coexpression_method = args.coexpression
    max_dim = args.max_dim

    
    # ---------------------------------------------------------
    # OPTIONAL: Apply DEG filtering only if BOTH thresholds are provided
    # ---------------------------------------------------------
    if args.padj_threshold is not None and args.log2fc is not None:

        print(f"Applying DEG filtering: padj < {args.padj_threshold}, |log2FC| > {args.log2fc}")

    expr_df = pd.read_csv(args.gene_expression_matrix)
    if preselection_genes_path is not None:
        preselection_df = pd.read_csv(preselection_genes_path, index_col=0)


    # --- Optional DEG filtering ---
    if args.padj_threshold is not None and args.log2fc is not None:
        if not {"padj", "log2FoldChange"}.issubset(preselection_df.columns):
            raise ValueError(
                "Preselection file must contain 'padj' and 'log2FoldChange' "
                "columns to use padj/log2fc filtering."
            )

        print(
            f"Applying DEG filtering: padj < {args.padj_threshold}, "
            f"|log2FC| > {args.log2fc}"
        )

        deg_filtered = preselection_df[
            (preselection_df["padj"] < args.padj_threshold)
            & (preselection_df["log2FoldChange"].abs() > args.log2fc)
        ]
        gene_list = deg_filtered.index.tolist()
        print(f"✅ Genes after DEG filtering: {len(gene_list)}")
        expr_sub = expr_df[gene_list].copy()
    else:
            print("ℹ️ No DEG filtering applied; using full preselection gene list.")
            expr_sub = expr_df
            gene_list = expr_df.columns.tolist()  # Use all genes from the expression matrix if no preselection provided
            print(f"✅ Using {len(gene_list)} genes from the expression matrix.")

    

    X = expr_sub.values  # shape (n_samples, n_genes)
    n_genes = X.shape[1]
    gene_dict = {i: gene_list[i] for i in range(n_genes)}


    print(f"✅ Loaded expression matrix with shape {expr_df.shape}")
    print(f"✅ Preselection genes: {len(gene_list)}")

    if coexpression_method == "dtem":
        coexpression_matrix = dtem.compute_dtem_matrix(X)
    elif coexpression_method == "wto_dc":
        coexpression_matrix = wto.compute_wto_matrix(X, adj_matrix="dc")

    elif coexpression_method == "wto_pearson":
        coexpression_matrix = wto.compute_wto_matrix(X, adj_matrix="pearson")

    elif coexpression_method == "dc":
        coexpression_matrix = dc.compute_distance_correlation_matrix(X)
    elif coexpression_method == "pearson":
        coexpression_matrix = pearson.compute_pearson_correlation_matrix(X)

    else:
        raise ValueError(f"Unsupported coexpression_method: {coexpression_method}")


    print("\n🏗 Running WG-TDA biomarker discovery...")
    results = tda_bio.run_biomarker_discovery(
        matrix=coexpression_matrix,
        gene_dict=gene_dict,
        max_dim=args.max_dim,
        save_plot=True,
        filename="persistence_diagram_biomarker.png",
        title=f"WG-TDA Biomarker Discovery",
    )
    print(results)
    results['interactions'].to_csv(args.interactions_output, index=False)
    print(f"✅ Saved interactions DataFrame to {args.interactions_output}")

if __name__ == "__main__":
    main()