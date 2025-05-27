from anndata import AnnData

from pydeconv import SignatureMatrix

from .base_model import SignatureBasedModel
from .solver import solver_ols


class OLS(SignatureBasedModel):
    """Ordinary Least Squares (OLS) model for deconvolution.

    Attributes
    ----------
    name : str
        Name of the model.
    signature_matrix : SignatureMatrix
        Signature matrix used for deconv.

    Methods
    -------

    transform(adata, layer, ratio=True)
        Fit the OLS model to the data and return the deconvolution values.
    """

    def __init__(self, signature_matrix: SignatureMatrix):
        """Initialize the OLS model.

        Parameters
        ----------
        signature_matrix : SignatureMatrix
            Signature matrix used for deconv.
        """
        super().__init__(name="OLS", signature_matrix=signature_matrix)

    def transform(self, adata: AnnData, layer: str = "relative_counts", ratio=True):
        """Fit the OLS model to the data and return the deconvolution values.

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
        if layer is None:
            layer = self.default_layer

        adata = self.valid_input(adata, layer=layer)
        values = solver_ols(adata, self._signature_matrix, layer)

        values = self.format_output(
            values, columns=self._signature_matrix.list_cell_types, index=adata.obs.index, ratio=ratio
        )

        return values
