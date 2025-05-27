from anndata import AnnData

from pydeconv import SignatureMatrix

from .base_model import SignatureBasedModel
from .solver import solver_nusvr, solver_svr


class SVR(SignatureBasedModel):
    def __init__(self, signature_matrix: SignatureMatrix):
        """Support Vector Regression (SVR) model for deconvolution.

        Parameters
        ----------
        signature_matrix : SignatureMatrix
            Signature matrix used for deconv.
        """
        super().__init__(name="SVR", signature_matrix=signature_matrix)

    def transform(self, adata: AnnData, layer: str = "relative_counts", ratio=False):
        """Fit the SVR model to the data and return the deconvolution values.

        Parameters
        ----------
        adata : AnnData
            Input anndata object.
        layer : str
            Layer of the data to use for deconvolution.
        ratio : bool, optional
            If True return a proportion of celltype instead raw values, by default False

        Returns
        -------
        pd.DataFrame
            Deconvolution values.
        """
        if layer is None:
            layer = self.default_layer

        adata = self.valid_input(adata, layer=layer)
        values = solver_svr(adata, self._signature_matrix, layer)

        values = self.format_output(
            values, columns=self._signature_matrix.list_cell_types, index=adata.obs.index, ratio=ratio
        )

        return values


class NuSVR(SignatureBasedModel):
    def __init__(self, signature_matrix: SignatureMatrix, norm=True, scale=False):
        """Support Vector Regression (SVR) model for deconvolution. NuSVR is similar to SVR but with a parameter to
        control the number of support vectors.

        Parameters
        ----------
        signature_matrix : SignatureMatrix
            Signature matrix used for deconv.
        norm : bool, optional
            if True, normalize the data, by default True
        scale : bool, optional
            if True, scale the data, by default False
        """
        super().__init__(name="NuSVR", signature_matrix=signature_matrix)
        self.norm = norm
        self.scale = scale

    def transform(self, adata: AnnData, layer: str = "relative_counts", ratio=False):
        """Fit the NuSVR model to the data and return the deconvolution values.

        Parameters
        ----------
        adata : AnnData
            Input anndata object.
        layer : str
            Layer of the data to use for deconvolution.
        ratio : bool, optional
            If True return a proportion of celltype instead raw values, by default False

        Returns
        -------
        pd.DataFrame
            Deconvolution values.
        """
        if layer is None:
            layer = self.default_layer

        adata = self.valid_input(adata, layer=layer)
        values = solver_nusvr(adata, self._signature_matrix, layer, norm=self.norm, scale=self.scale)

        values = self.format_output(
            values,
            columns=self._signature_matrix.list_cell_types,
            index=adata.obs.index,
            ratio=ratio,
            only_positive=True,
        )

        return values
