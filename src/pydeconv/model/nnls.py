import pandas as pd
from anndata import AnnData

from pydeconv import SignatureMatrix

from .base_model import SignatureBasedModel
from .solver import solver_nnls


class NNLS(SignatureBasedModel):
    """Non-Negative Least Squares (NNLS) model for deconvolution.
    Attributes
    ----------
    name : str
        Name of the model.
    signature_matrix : SignatureMatrix
        Signature matrix used for deconv.

    Methods
    -------
    transform(adata, layer, ratio=False)
        Fit the NNLS model to the data and return the deconvolution values.
    """

    def __init__(self, signature_matrix: SignatureMatrix):
        """
        Initialize the NNLS model.

        Parameters
        ----------
        signature_matrix : SignatureMatrix
            Signature matrix used for deconv.
        """
        super().__init__(name="NNLS", signature_matrix=signature_matrix)

    def transform(self, adata: AnnData, layer: str = "relative_counts", ratio=False) -> pd.DataFrame:
        """Fit the NNLS model to the data and return the deconvolution values.

        Parameters
        ----------
        adata : AnnData
            Anndata object containing the data.
        layer : str
            Layer of the data to use for deconvolution.
        ratio : bool, optional
            Transform raw output into cell proportions, by default False

        Returns
        -------
        pd.DataFrame
            Deconvolution values.
        """

        adata = self.valid_input(adata, layer=layer)
        values = solver_nnls(adata, self._signature_matrix, layer)

        values = self.format_output(
            values, columns=self._signature_matrix.list_cell_types, index=adata.obs.index, ratio=ratio
        )

        return values
