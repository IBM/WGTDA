import numpy as np
import pandas as pd


def remove_inf_rows(df: pd.DataFrame, column: str = "lifespan") -> pd.DataFrame:
    """
    Remove rows where a target column contains +inf or -inf.

    Parameters
    ----------
    df : pd.DataFrame
        Input interactions DataFrame.
    column : str
        Name of the column to check for infinite values.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame without infinite values.
    """
    return  df[df["lifespan"] != np.inf]


def filter_by_list_length(
    df: pd.DataFrame,
    column: str = "geneset",
    min_length: int = 3,
) -> pd.DataFrame:
    """
    Keep rows where the list inside a column has >= min_length elements.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    column : str
        Name of the column containing lists.
    min_length : int
        Minimum length required.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame.
    """
    return df[df[column].apply(lambda x: len(x) >= min_length)]


def extract_top_n_by_column(
    df: pd.DataFrame,
    group_col: str = "betti_number",
    score_col: str = "lifespan",
    top_n_percent: float = 10.0,
) -> pd.DataFrame:
    """
    Extract top N% rows *within each group* based on a score column.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame.
    group_col : str
        Column that defines groups (e.g., Betti numbers).
    score_col : str
        Column to score rows by (sorted descending).
    top_n_percent : float
        Percentage (0 < n ≤ 100).

    Returns
    -------
    pd.DataFrame
        DataFrame containing top rows from each group.
    """
    if not (0 < top_n_percent <= 100):
        raise ValueError("top_n_percent must be between 0 and 100.")

    result = []
    for val, group_df in df.groupby(group_col):
        k = max(1, int(len(group_df) * (top_n_percent / 100)))
        result.append(group_df.nlargest(k, score_col))

    return pd.concat(result, ignore_index=True)


def filter_interactions(
    df: pd.DataFrame,
    lifespan_col: str = "lifespan",
    geneset_col: str = "geneset",
    min_genes: int = 3,
    top_n_percent: float = 10.0,
) -> pd.DataFrame:
    """
    Full WGTDA topological filtering pipeline.

    Steps:
    1. Remove infinite persistence
    2. Filter by geneset size
    3. Extract top N% persistent holes per Betti number

    Parameters
    ----------
    df : pd.DataFrame
        Input interactions DataFrame.
    lifespan_col : str
        Column measuring persistence lifespan.
    geneset_col : str
        Column containing list of genes per topological feature.
    min_genes : int
        Minimum geneset size allowed.
    top_n_percent : float
        Percentage of top persistent features to keep.

    Returns
    -------
    pd.DataFrame
        Filtered topological interactions.
    """
    df_clean = remove_inf_rows(df, column=lifespan_col)
    df_clean = filter_by_list_length(df_clean, column=geneset_col, min_length=min_genes)
    df_clean = extract_top_n_by_column(
        df_clean,
        group_col="betti_number",
        score_col=lifespan_col,
        top_n_percent=top_n_percent,
    )

    return df_clean
