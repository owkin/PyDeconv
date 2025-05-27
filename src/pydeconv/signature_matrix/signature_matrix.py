from pathlib import Path

import pandas as pd


class SignatureMatrix:
    """Signature matrix object.

    Attributes
    ----------
    values : pd.DataFrame
        Signature matrix.
    list_genes : list
        List of genes.
    list_cell_types : list
        List of cell types.

    Methods
    -------
    load(path, **kwargs)
        Load a signature matrix from a file.
    """

    def __init__(self, df: pd.DataFrame):
        """Signature matrix object.

        Parameters
        ----------
        df : pd.DataFrame
            Signature matrix, rows are genes, columns are cell types.
        """
        self._df = df

    @classmethod
    def load(clf, path: str, **kwargs):
        """Load a signature matrix from a file. Only support csv file for now."""
        df = load_signature_matrix(Path(path), **kwargs)
        return clf(df)

    @property
    def values(self):
        """Signature matrix."""
        return self._df

    @property
    def list_genes(self):
        """List of genes."""
        return self._df.index.to_list()

    @property
    def list_cell_types(self):
        """List of cell types."""
        return self._df.columns.to_list()

    def update_gene_list(self, new_gene_list):
        """Update the gene list."""
        self._df = self._df.loc[new_gene_list]


def load_signature_matrix(path: Path, **kwargs):
    """Load a signature matrix from a file.

    Parameters
    ----------
    path : Path
        Path to the signature matrix file.
    **kwargs
        Additional arguments to pass to the read_csv function.

    Returns
    -------
    pd.DataFrame
        Signature matrix.
    """
    if path.suffix == ".csv":
        df = pd.read_csv(path, **kwargs)
    else:
        raise ValueError("File format not supported")
    return df
