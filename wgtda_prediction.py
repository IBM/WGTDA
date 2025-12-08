import src.tda.pipeline_prediction as tda_pred
import argparse
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, f1_score
import src.correlation.dtem as dtem
import src.correlation.wto as wto
import src.correlation.pearson as pearson
import src.correlation.distance_corr as dc
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="Run the prediction pipeline.")
    parser.add_argument("--gene_expression_matrix", "-p", type=str, required=True,
                        help="Path to the gene expression matrix file.")
    parser.add_argument("--preselection_genes", "-pp", type=str, required=True,
                        help="Path to the preselection genes file.")
    parser.add_argument("--padj_threshold", type=float, default=0.05,
                        help="Adjusted p-value threshold for gene selection."
                )
    parser.add_argument("--log2fc", type=int, default=3.5, 
                        help="Log2 fold change threshold for gene selection."
                )
    parser.add_argument("--coexpression", "-c", type=str, default="dtem",
                        choices=["wto_dc", "wto_pearson", "dc", "pearson", "dtem"])
    
    parser.add_argument("--max_dim", type=int, default=2,
                        help="Maximum homology dimension to compute.")
    parser.add_argument("--save_landscapes", "-sl", default="landscapes", type=str,
                        help="Path to save the computed TDA landscapes.")
    parser.add_argument("--num_layers", type=int, default=2,
                        help="Number of landscape layers.")
    parser.add_argument("--resolution", type=int, default=50,
                        help="Resolution of the landscapes.")
    
    
    return parser.parse_args() 


def train_rf(X_train, y_train, X_test, y_test, label: str):
    """Train a Random Forest and print metrics."""
    clf = RandomForestClassifier()

    print(f"\n🧬 Training Random Forest on {label}...")
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)

    acc = accuracy_score(y_test, y_pred)
    f1w = f1_score(y_test, y_pred, average="weighted")

    print(f"➡ Results ({label}):")
    print(f"   Accuracy   : {acc:.3f}")
    print(f"   F1 weighted: {f1w:.3f}")

    return clf, acc, f1w


def main():
    args = parse_args()
    # Here you would typically load your data using the provided file paths
    # For example:
    expression_matrix_path = args.gene_expression_matrix
    preselection_genes_path = args.preselection_genes
    coexpression_method = args.coexpression
    padj_threshold = args.padj_threshold
    log2fc = args.log2fc
    save_landscapes = args.save_landscapes
    max_dim = args.max_dim
    num_layers = args.num_layers
    resolution = args.resolution


    # Load your data using the provided file paths
    expression_matrix = pd.read_csv(expression_matrix_path)
    preselection_genes = pd.read_csv(preselection_genes_path, index_col=0).index.tolist()
    print(preselection_genes)
    phenotype = expression_matrix["phenotype"]
    gene_list = preselection_genes
    print(gene_list)
    sig_exp_matrix = expression_matrix[gene_list].copy()
    sig_exp_matrix["phenotype"] = phenotype.values


    train_df, test_df = train_test_split(
        sig_exp_matrix,
        test_size=0.2,
        random_state=42,
        stratify=sig_exp_matrix["phenotype"]
    )

    print(train_df)

    X_train_expr = train_df.drop(columns=["phenotype"]).values
    X_test_expr = test_df.drop(columns=["phenotype"]).values
    y_train = train_df["phenotype"].values
    y_test = test_df["phenotype"].values

    le = LabelEncoder()
    y_train = le.fit_transform(y_train)
    y_test = le.transform(y_test)

     # 1) Random Forest on raw expression
    _expr_clf, expr_acc, expr_f1 = train_rf(
        X_train_expr, y_train, X_test_expr, y_test,
        label="Raw Gene Expression (DEG)"
    )


    # Compute coexpression matrix for the expression matrix
    print(f"\n Computing {coexpression_method} matrix on training data...")
    if coexpression_method == "dtem":
        coexpression_matrix = dtem.compute_dtem_matrix(X_train_expr)
    elif coexpression_method == "wto_dc":
        coexpression_matrix = wto.compute_wto_matrix(X_train_expr, adj_matrix="dc")

    elif coexpression_method == "wto_pearson":
        coexpression_matrix = wto.compute_wto_matrix(X_train_expr, adj_matrix="pearson")

    elif coexpression_method == "dc":
        coexpression_matrix = dc.compute_distance_correlation_matrix(X_train_expr)
    elif coexpression_method == "pearson":
        coexpression_matrix = pearson.compute_pearson_correlation_matrix(X_train_expr)


    print("\n🏗 Computing TDA landscapes (train)...")
    landscapes_train = tda_pred.run_prediction_pipeline(
        X_train_expr,
        coexpression_matrix,
        num_layers=num_layers,
        resolution=resolution,
        max_dim=max_dim,
    )
    
    print("\n🏗 Computing TDA landscapes (test)...")
    landscapes_test = tda_pred.run_prediction_pipeline(
        X_test_expr,
        coexpression_matrix,
        num_layers=num_layers,
        resolution=resolution,
        max_dim=max_dim,
    )
    landscapes_train = np.nan_to_num(landscapes_train)
    landscapes_test = np.nan_to_num(landscapes_test)

    # Optionally save landscapes
    np.save(f"{save_landscapes}_train.npy", landscapes_train)
    np.save(f"{save_landscapes}_test.npy", landscapes_test)

    # use landscapes to predict 
    print("\n🧬 Training on TDA Landscapes Features..."
        )
     # 3) Random Forest on TDA landscapes
    _tda_clf, tda_acc, tda_f1 = train_rf(
        landscapes_train, y_train, landscapes_test, y_test,
        label="TDA Landscapes (DEG)"
    )

    print("\n✅ Summary:")
    print(f"Raw Expression RF  -> Acc: {expr_acc:.3f}, F1w: {expr_f1:.3f}")
    print(f"TDA Landscapes RF  -> Acc: {tda_acc:.3f}, F1w: {tda_f1:.3f}")

if __name__ == "__main__":
    main()
                    
