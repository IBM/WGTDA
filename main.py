import argparse
import os
import warnings

from wgtda import (construct_vr_complex_rna_matrix,
                   convert_gene_exp_to_array_and_dict, filter_genes,
                   flatten_gene_list, interactions_dataframe,
                   load_gene_expression_data)
from wgtda.correlation import (compute_distance_correlation_matrix,
                               compute_wto_matrix)
from wgtda.filters import extract_top_n_persistent_holes, remove_infinite_holes

warnings.simplefilter(action="ignore", category=FutureWarning)


def parse_args():
    """
    Function to parse the arguments for wgtda
    """
    parser = argparse.ArgumentParser(description="Main program for WGTDA")

    # data rel
    parser.add_argument(
        "--file_path",
        "-p",
        type=str,
        default="data/TCGA/BRCA.pkl",
        help="The path to the data file. Supported file types are CSV (.csv), Excel (.xls, .xlsx, .xlsm, "
        ".xlsb), Pickle (.pkl, .pickle), and TSV/Text (.tsv, .txt).",
    )

    parser.add_argument(
        "--filter_genes_path",
        "-fg",
        type=str,
        default="data/preselection/cancer_genes.txt",
        help="The path to a txt file containing selected genes to use in WGTDA.",
    )

    parser.add_argument(
        "--preprocessing",
        "-pp",
        type=str,
        default="dc",
        help="The path to select which method to use for "
        "preprocessing of gene expression data"
        "('dc' for distance correlation or 'stom' for signed TOMs",
    )

    parser.add_argument(
        "--outputdir",
        "-o",
        type=str,
        default="./output/",
        help="The path for the output interactions.csv",
    )

    parser.add_argument(
        "--dimensions",
        "-d",
        type=int,
        default=3,
        help="specify how many dimensions that user wants to input",
    )

    # Topological Filters
    parser.add_argument(
        "--remove_inf_values",
        "-inf",
        type=bool,
        default=True,
        help="Remove holes that do not close",
    )

    parser.add_argument(
        "--filter_persistence",
        "-fp",
        type=int,
        default=10,
        help="Filter top n% of persistent features",
    )

    return parser.parse_args()


def main():
    """
    Function to execute WGTDA
    """
    args = parse_args()
    output = args.outputdir
    remove_inf_values = args.remove_inf_values
    dimensions = args.dimensions

    # Check args
    preprocessing_funcs = ["dc", "stom"]
    if args.preprocessing not in preprocessing_funcs:
        return

    print("Loading Gene Expression Data")

    df = load_gene_expression_data(args.file_path)
    print("Preselecting Genes from " + args.filter_genes_path)
    gene_exp_df = filter_genes(df, args.filter_genes_path)
    gene_exp_arr, gene_dict = convert_gene_exp_to_array_and_dict(gene_exp_df)

    if args.preprocessing == "dc":
        print("Computing the distance correlation matrix")
        dist_matrix = compute_distance_correlation_matrix(gene_exp_arr=gene_exp_arr)
    elif args.preprocessing == "stom":
        print("Computing the weighted signed topological overlapping matrix")
        dist_matrix = compute_wto_matrix(gene_exp_arr=gene_exp_arr)
    else:
        raise ValueError("Unsupported or Unknown preprocessing method.")

    print("Constructing the Vietoris Rips Complex")
    persistence, rips_complex = construct_vr_complex_rna_matrix(dist_matrix, dimensions)
    interactions = interactions_dataframe(persistence, rips_complex, gene_dict)

    if remove_inf_values:
        print("Filtering  Holes")
        interactions = remove_infinite_holes(interactions)

    interactions = extract_top_n_persistent_holes(interactions, args.filter_persistence)

    if not os.path.exists(output):
        # If the directory does not exist, create it
        os.makedirs(output)

    # Apply the function to the 'gene_interactions' column
    interactions["gene_set"] = interactions["vertices_set"].apply(flatten_gene_list)

    interactions.to_csv(output + "interactions.csv", index=True)

    print("Saved to " + output + "interactions.csv")


if __name__ == "__main__":
    main()
