import pandas as pd


def create_gene_list_file(gene_expression_df: pd.DataFrame, column_name: str, output_gene_list_file: str):
    """
    Create a text file with all the genes from a specified column in the DataFrame.

    Parameters:
    gene_expression_df (pd.DataFrame): Gene expression DataFrame.
    column_name (str): The name of the column containing gene names.
    output_gene_list_file (str): Path to the output text file.

    Returns:
    None
    """
    try:
        # Attempt to access the specified column
        genes = gene_expression_df[column_name].dropna().unique()
    except KeyError as e:
        # Handle the case where the specified column does not exist
        raise KeyError(
            f"The column '{column_name}' was not found in the DataFrame. Please check the column name and try again."
        ) from e

    with open(output_gene_list_file, "w") as file:
        for gene in genes:
            file.write(f"{gene}\n")


# Example usage
if __name__ == "__main__":
    # Load the gene expression data
    gene_expression_df = pd.read_pickle("data/TCGA/BRCA.pkl")
    preselection_df = pd.read_csv("data/preselection/drug_metabolism.csv", sep=",")

    print(preselection_df)

    # Create a gene list file from a column
    create_gene_list_file(
        preselection_df, "Symbol", "data/preselection/drug_metabolism.txt"
    )
