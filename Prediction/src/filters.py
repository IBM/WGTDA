import numpy as np
import pandas as pd

def remove_inf_values(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove rows with infinite values from the lifespan.

    Parameters:
    - df : pandas.DataFrame
        The DataFrame from which to remove rows.

    Returns:
    - pandas.DataFrame
        A new DataFrame with rows containing infinite values in the lifespan column.
    """
    # Exclude 'inf' values from the specified column
    filtered_df = df[df["lifespan"] != np.inf]

    return filtered_df

def extract_top_n_persistent_holes(df: pd.DataFrame, n=10)-> pd.DataFrame:
    """
    Extract the top n% rows for each Betti number based on the specified score column.

    Parameters:
    df (pd.DataFrame): The input DataFrame.
    n (float): The top percentage to extract (0 < n <= 100).
    score_column (str): The column name to base the extraction on.

    Returns:
    pd.DataFrame: A DataFrame containing the top n% rows for each Betti number.
    """
    if n <= 0 or n > 100:
        raise ValueError("n must be between 0 and 100.")

    # Create an empty DataFrame to store the results
    top_n_percent_df = pd.DataFrame()

    # Loop through each unique Betti number
    for betti in df["betti_number"].unique():
        # Filter the DataFrame for the current Betti number
        betti_df = df[df["betti_number"] == betti]

        # Calculate the number of top rows to extract
        top_n_count = int(len(betti_df) * (n / 100))

        # Extract the top n% rows based on the score column
        top_n_df = betti_df.nlargest(top_n_count, "lifespan")

        # Append the top n% rows to the results DataFrame
        top_n_percent_df = pd.concat([top_n_percent_df, top_n_df])

    return top_n_percent_df

def filter_rows_by_list_length(df, min_length):
    """
    Filters rows in a DataFrame based on the number of elements in a list within a specific column.
    
    Parameters:
    - df (pd.DataFrame): The input DataFrame.
    - column_name (str): The name of the column containing lists.
    - min_length (int): The minimum allowed number of elements in the list.
    
    Returns:
    - pd.DataFrame: A filtered DataFrame containing rows where the list length is <= max_length.
    """
    return df[df["geneset"].apply(len) >= min_length]


def filtering(df, top_n=10, min_length=3):
    filtered_df = remove_inf_values(df=df)
    filtered_df = filter_rows_by_list_length(df=filtered_df, min_length=min_length)
    top_df = extract_top_n_persistent_holes(df=filtered_df, n=top_n)

    return top_df
    
    