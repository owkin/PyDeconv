import pandas as pd
from scipy.stats import pearsonr


def correlation(x: pd.DataFrame, y: pd.DataFrame) -> pd.Series:
    """
    Calculate the correlation between two lists of numbers.

    Args:
        x (list): First list of numbers.
        y (list): Second list of numbers.

    Returns:
        float: The Pearson correlation coefficient between the two lists.
    """
    if not isinstance(x, pd.DataFrame) or not isinstance(y, pd.DataFrame):
        raise ValueError("Both x and y must be pandas DataFrames.")

    x, y = check_inputs(x, y)

    return pd.Series({col: pearsonr(x[col], y[col])[0] for col in x.columns})


def mse(x: pd.DataFrame, y: pd.DataFrame):
    """
    Calculate the Mean Square Error (MSE) between two lists of numbers.

    Args:
        x (list): First list of numbers.
        y (list): Second list of numbers.

    Returns:
        float: The RMSE between the two lists.
    """
    x, y = check_inputs(x, y)

    return ((x - y) ** 2).mean()


def rmse(x: pd.DataFrame, y: pd.DataFrame):
    """
    Calculate the Root Mean Square Error (RMSE) between two lists of numbers.

    Args:
        x (list): First list of numbers.
        y (list): Second list of numbers.

    Returns:
        float: The RMSE between the two lists.
    """
    x, y = check_inputs(x, y)

    return ((x - y) ** 2).mean() ** 0.5


def check_inputs(df1: pd.DataFrame, df2: pd.DataFrame):
    check_structure = (set(df1.columns) == set(df2.columns)) and (set(df1.index) == set(df2.index))
    if not check_structure:
        raise ValueError("DataFrames must have the same structure (columns and index).")

    # Format columns if not in the same order
    if not df1.columns.equals(df2.columns):
        df2 = df2[df1.columns]

    if not df1.index.equals(df2.index):
        df2 = df2.reindex(df1.index)
    return df1, df2
