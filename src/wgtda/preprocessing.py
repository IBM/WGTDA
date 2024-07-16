import ast
import os
from typing import TextIO

import numpy as np
import pandas as pd


def load_gene_expression_data(file_path: str) -> pd.DataFrame:
    """
    Load gene expression data from a file, automatically determining the file type.

    This function reads a file specified by `file_path` and loads it into a pandas DataFrame. The file type is inferred
    from the file extension. Supported file types are CSV, Excel, Pickle, and TSV/Text.

    Parameters:
    - file_path : str
        The path to the data file.

    Returns:
    - pandas.DataFrame
        The data loaded into a pandas DataFrame.

    Raises:
    - ValueError
        If the file extension is not recognized as a supported file type.
    """

    # Infer the file type from the file extension
    _, file_extension = os.path.splitext(file_path)
    file_extension = file_extension.lower()

    if file_extension == ".csv":
        df = pd.read_csv(file_path)
    elif file_extension in [".xls", ".xlsx", ".xlsm", ".xlsb"]:
        df = pd.read_excel(file_path)
    elif file_extension in [".pkl", ".pickle"]:
        df = pd.read_pickle(file_path)
    elif file_extension in ".tsv":
        df = pd.read_csv(file_path, sep="\t")
    else:
        raise ValueError(
            "Unsupported file extension. Supported extensions are: .csv, .xls(x), .pkl(pickle), .tsv"
        )
    return df


def filter_genes(
    gene_expression_df: pd.DataFrame, gene_list_file: TextIO
) -> tuple[np.ndarray, dict]:
    """
    Filter gene expression data based on a txt list of genes.

    Parameters:
    gene_expression_df (pd.DataFrame): Gene expression DataFrame where columns are gene names.
    gene_list_file (str): Path to the text file containing the list of genes to filter.


    Returns:
    - pandas.DataFrame
        The filtered DataFrame containing only the preselected genes.
    - dict
        A dictionary mapping column indices to gene names in the filtered DataFrame.
    """
    # Read the list of genes from the text file
    with open(gene_list_file, "r") as file:
        gene_list = file.read().splitlines()

    # Filter the DataFrame to include only the genes from the list
    filtered_df = gene_expression_df[
        gene_expression_df.columns.intersection(gene_list[0:50])
    ]

    # Count the number of genes that were successfully filtered
    num_filtered_genes = len(filtered_df.columns)

    print("Number of genes filtered: ", num_filtered_genes)

    return filtered_df


def convert_gene_exp_to_array_and_dict(
    filtered_df: pd.DataFrame,
) -> tuple[np.ndarray, dict[int, str]]:
    """
    Convert the filtered DataFrame to a NumPy array and create a dictionary mapping indices to gene names.

    Parameters:
    filtered_df (pd.DataFrame): Filtered gene expression DataFrame.

    Returns:
    Tuple containing:
    - np.ndarray: The DataFrame converted to a NumPy array.
    - Dict[int, str]: A dictionary mapping column indices to gene names.
    """
    # Create the dictionary mapping indices to gene names
    gene_dict = dict(enumerate(filtered_df.columns))

    # Convert the DataFrame to a NumPy array

    return np.array(filtered_df), gene_dict


def create_gene_list_file(
    gene_expression_df: pd.DataFrame, column_name: str, output_gene_list_file: TextIO
):
    """
    Create a text file with all the genes from a specified column in the DataFrame.

    Parameters:
    gene_expression_df (pd.DataFrame): Gene expression DataFrame.
    column_name (str): The name of the column containing gene names.
    output_gene_list_file (str): Path to the output text file.

    Returns:
    None
    """
    genes = gene_expression_df[column_name].dropna().unique()

    with open(output_gene_list_file, "w") as file:
        for gene in genes:
            file.write(f"{gene}\n")


# Function to convert sting['vertices_set'] to a flattened list of genes in ['gene_set']
def flatten_gene_list(gene_list_str):
    try:
        # Convert the string to an actual list of lists
        list_of_lists = ast.literal_eval(gene_list_str)
        # Flatten the list of lists into a single list
        flattened_list = list(
            set(gene for sublist in list_of_lists for gene in sublist)
        )
        return flattened_list
    except (ValueError, SyntaxError):
        return []
