import os
os.environ["LIGHTGBM_VERBOSE"] = "0"
import datetime
import sys
import copy
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import LabelEncoder

import traceback
from collections import defaultdict
from tqdm import tqdm
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, mean_squared_error, log_loss
from sklearn.model_selection import GridSearchCV
from helper import log_error, suppress_stdout_stderr,setup_dual_logging
setup_dual_logging("cancer_staging")

from src.functions import (compute_distance_correlation_matrix,
 compute_simplicial_complex_and_landscapes, compute_wto_matrix, compute_pearson_correlation_matrix,
  patient_correlation_measure)
from lightgbm import LGBMClassifier
from train_models import train_all_models

import warnings
warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in subtract",
    category=RuntimeWarning
)
import os
os.environ["LIGHTGBM_VERBOSE"] = "0"

cancer_types = ["BLCA", "BRCA", "COAD", "HNSC", "KIRC", "LUAD", "LUSC", "SKCM", "STAD", "THCA"]

THRESHOLD = {
    "BLCA": {"padj": 0.05, "log2fc": 2.5},
    "BRCA": {"padj": 0.05, "log2fc": 2.5},
    "COAD": {"padj": 0.05, "log2fc": 1.5},
    "HNSC": {"padj": 0.05, "log2fc": 2.0},
    "KIRC": {"padj": 0.05, "log2fc": 1.5},
    "LUAD": {"padj": 0.05, "log2fc": 2.0},
    "LUSC": {"padj": 0.05, "log2fc": 2.0},
    "SKCM": {"padj": 0.05, "log2fc": 1.5},
    "STAD": {"padj": 0.05, "log2fc": 2.0},
    "THCA": {"padj": 0.05, "log2fc": 1.5},
}

methods = {
    "pearson": compute_pearson_correlation_matrix,
    "distcorr": compute_distance_correlation_matrix,
    "wto_pearson": lambda X: compute_wto_matrix(X, "pearsons"),
    "wto_distcorr": lambda X: compute_wto_matrix(X, "dc")
}

stages = {
    "Stage I": 1,
    "Stage II": 2,
    "Stage III": 3,
    "Stage IV": 4
}

num_landscape = 2
resolution = 100

DEBUG = False

if DEBUG == True:
    print("DEBUG mode is enabled.")
    cancer_types = ["BRCA", "LUAD"]
    methods = {
    "pearson": compute_pearson_correlation_matrix}

all_tda_results_by_cancer = {}
all_expr_results = {}
results_list = []

if __name__ == "__main__":
    for cancer_type in cancer_types:
        all_tda_results = {}
        print(f"\n--- Processing: {cancer_type} ---")

        if cancer_type not in THRESHOLD:
            print(f"No threshold defined for {cancer_type}, skipping.")
            continue

        # Get thresholds for this cancer type
        padj = THRESHOLD[cancer_type]["padj"]
        log2fc = THRESHOLD[cancer_type]["log2fc"] + 0.5
        print(f"\nThresholds: padj={padj}, log2fc={log2fc}")

        # ---------- File Paths ----------
        base_path = os.path.join("TCGA_Processed", cancer_type)
        staging_file = os.path.join(base_path, f"{cancer_type}_staging_clean.csv")
        sig_genes_file = os.path.join(base_path, f"{cancer_type}_sig_genes.csv")
        fpkm_file = os.path.join(base_path, f"{cancer_type}_fpkm_counts.csv")

        # ---------- Load Data ----------
        staging_clean = pd.read_csv(staging_file)
        sig_genes = pd.read_csv(sig_genes_file, index_col=0)
        fpkm = pd.read_csv(fpkm_file)

        # ---------- Merge Phenotype ----------
        merger = staging_clean[["Unnamed: 0", "ajcc_pathologic_stage_grouped"]]
        expression_matrix = fpkm.merge(merger, on="Unnamed: 0", how="inner")

        # ---------- Extract and Clean ----------
        phenotype = expression_matrix["ajcc_pathologic_stage_grouped"]
        expression_matrix = expression_matrix.drop(columns=["Unnamed: 0", "ajcc_pathologic_stage_grouped"])

        # ---------- Filter Significant Genes ----------
        deg = sig_genes[(sig_genes["padj"] < padj) & (sig_genes["log2FoldChange"].abs() > log2fc)]
        gene_list = deg.index.intersection(expression_matrix.columns)

        if len(gene_list) == 0:
            raise ValueError(f"\nNo significant genes found for {cancer_type} with current thresholds.")

        # ---------- Subset and Append Phenotype ----------
        sig_exp_matrix = expression_matrix[gene_list].copy()
        sig_exp_matrix["phenotype"] = phenotype.values

        print(f"\n[{cancer_type}] Data loaded, filtered and phenotype assigned. {len(gene_list)} genes retained.")

        print(sig_exp_matrix.head(2))
        # Step 1: Split the entire dataset first
        # Step 1: Split the entire dataset
        train_df, test_df = train_test_split(
        sig_exp_matrix,
        test_size=0.2,
        random_state=42,
        stratify=sig_exp_matrix["phenotype"]
    )
        
        X_train_expr = train_df.drop(columns=["phenotype"]).values
        X_test_expr = test_df.drop(columns=["phenotype"]).values
        y_train_expr = train_df["phenotype"].map(stages).values
        y_test_expr = test_df["phenotype"].map(stages).values

        le = LabelEncoder()
        y_train_expr = le.fit_transform(y_train_expr)
        y_test_expr = le.transform(y_test_expr)

        # Train on gene expression
        print("\n🧬 Training on Raw Gene Expression Features...")
        expr_results = train_all_models(
        X_train_expr, y_train_expr, X_test_expr, y_test_expr,
        cancer_type, feature_type="raw_expression"
)

        # Step 2: Loop over each method
        for method_name, method_fn in methods.items():
            print(f"\n=== Now Running with WGTDA method: {method_name} ===\n")

            train_landscapes = []
            test_landscapes = []
            train_labels = []
            test_labels = []

            for stage, label in stages.items():
                # Stage-specific train/test sets
                stage_train = train_df[train_df["phenotype"] == stage].drop(columns=["phenotype"])
                stage_test = test_df[test_df["phenotype"] == stage].drop(columns=["phenotype"])


                if len(stage_train) < 2 or len(stage_test) < 2:
                    print(f" Skipping {stage} — too few samples.")
                    continue

                # Compute correlation matrix on training set
                corr_matrix = method_fn(stage_train.values)

                # Compute landscapes
                print(f"\nComputing simplicial complex, persistent homology and landscape for training and testing {stage}...")
                train_land = compute_simplicial_complex_and_landscapes(stage_train.values, corr_matrix, num_landscape, resolution)
                test_land = compute_simplicial_complex_and_landscapes(stage_test.values, corr_matrix, num_landscape, resolution)

                # after computing train_land, test_land
                train_land = np.nan_to_num(train_land, nan=0.0, posinf=0.0, neginf=0.0)
                test_land  = np.nan_to_num(test_land,  nan=0.0, posinf=0.0, neginf=0.0)

                # Append data
                train_landscapes.append(train_land)
                test_landscapes.append(test_land)
                train_labels.append(np.full(train_land.shape[0], label))
                test_labels.append(np.full(test_land.shape[0], label))

            # Combine all stages
            if train_landscapes and test_landscapes:
                X_train = np.vstack(train_landscapes)
                X_test = np.vstack(test_landscapes)
                y_train = np.concatenate(train_labels)
                y_test = np.concatenate(test_labels)

                le = LabelEncoder()
                y_train = le.fit_transform(y_train)
                y_test = le.transform(y_test)


                print(f"\nMethod: {method_name}")
                print("Training Landscapes shape:", X_train.shape)
                print("Training Labels shape    :", y_train.shape)
                print("Testing Landscapes shape :", X_test.shape)
                print("Testing Labels shape     :", y_test.shape)

                # Train on TDA features
                print("\n🏗️ Training on TDA Landscape Features...")
                tda_results = train_all_models(
                 X_train, y_train, X_test, y_test,
                cancer_type, feature_type=method_name
)
                all_tda_results[method_name] = tda_results


            

        print(f"\n Raw Gene Expression Results for {cancer_type}:")
        for model_name, metrics in expr_results.items():
            print(f"\n--- {model_name.upper()} ---")
            print(f"Cancer Type: {cancer_type}")
            print("Accuracy:", metrics["accuracy"])

        print(f"\n All TDA Feature Results for {cancer_type}:")
        for method_name, tda_results in all_tda_results.items():
            print(f"\n Method: {method_name}")
            for model_name, metrics in tda_results.items():
                print(f"--- {model_name.upper()} ---")
                print(f"Cancer Type: {cancer_type}")
                print("Accuracy:", metrics["accuracy"])

        all_expr_results[cancer_type] = copy.deepcopy(expr_results)
        all_tda_results_by_cancer[cancer_type] = copy.deepcopy(all_tda_results)

    summary_rows = []
    for cancer, expr in all_expr_results.items():
        for model_name, metrics in expr.items():
            summary_rows.append({
                "Cancer": cancer,
                "Feature": "raw_expression",
                "Model": model_name,
                "Accuracy": metrics["accuracy"],
                "F1_weighted": metrics["f1_weighted"]
            })
    for cancer, method_dict in all_tda_results_by_cancer.items():
        for feature, models in method_dict.items():
            for model_name, metrics in models.items():
                summary_rows.append({
                    "Cancer": cancer,
                    "Feature": feature,
                    "Model": model_name,
                    "Accuracy": metrics["accuracy"],
                    "F1_weighted": metrics["f1_weighted"]
                })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(f"results/final_summary_with_f1{datetime.datetime.now().strftime('%Y%m%d_%H%M')}.csv", index=False)

    print("\n FINAL SUMMARY")
    print(summary_df.to_string(index=False))



        




                



                    


