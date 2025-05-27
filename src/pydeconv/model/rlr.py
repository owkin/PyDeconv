from anndata import AnnData

from pydeconv import SignatureMatrix

from .base_model import SignatureBasedModel
from .solver import solver_rlr


class RLR(SignatureBasedModel):
    """Robust Linear Regression (RLR) model for deconvolution.

    Attributes
    ----------
    name : str
        Name of the model.
    signature_matrix : SignatureMatrix
        Signature matrix used for deconv.

    Methods
    -------

    transform(adata, layer, ratio=True)
        Fit the RLR model to the data and return the deconvolution values.
    """

    def __init__(self, signature_matrix: SignatureMatrix):
        """Initialize the RLR model.

        Parameters
        ----------
        signature_matrix : SignatureMatrix
            Signature matrix used for deconv.
        """
        super().__init__(name="RLR", signature_matrix=signature_matrix)

    def transform(self, adata: AnnData, layer: str = "relative_counts", ratio=True):
        """Fit the RLR model to the data and return the deconvolution values.

        Parameters
        ----------
        adata : AnnData
            Anndata object containing the data.
        layer : str
            Layer of the data to use for deconvolution.
        ratio : bool, optional
            Transform raw output into cell proportions, by default True

        Returns
        -------
        pd.DataFrame
            Deconvolution values.
        """

        adata = self.valid_input(adata, layer=layer)
        values = solver_rlr(adata, self._signature_matrix, layer)

        values = self.format_output(
            values, columns=self._signature_matrix.list_cell_types, index=adata.obs.index, ratio=ratio
        )

        return values
