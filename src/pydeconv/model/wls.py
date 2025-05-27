from typing import Literal, Union

from anndata import AnnData

from pydeconv import SignatureMatrix

from .base_model import SignatureBasedModel
from .solver import solver_nnls, solver_ols, solver_svr, solver_wls


SOLVER_METHODS = {"nnls": solver_nnls, "ols": solver_ols, "svr": solver_svr}


class WLS(SignatureBasedModel):
    def __init__(self, signature_matrix: SignatureMatrix, solver_method: Literal["ols", "nnls", "svr"] = "nnls"):
        """Weighted Least Squares (WLS) model for deconvolution.

        Parameters
        ----------
        signature_matrix : SignatureMatrix
            Signature matrix used for deconv.
        solver_method : Literal[ols;, nnls, svr], optional


        Raises
        ------
        ValueError
            Raised when solver_method is not found.
        """
        super().__init__(name=f"wls_{solver_method}", signature_matrix=signature_matrix)
        self._baseline_solver = SOLVER_METHODS.get(solver_method)
        if self._baseline_solver is None:
            raise ValueError(f"Solver {solver_method} not found")

    def transform(
        self,
        adata: AnnData,
        layer: str = "relative_counts",
        ratio: bool = False,
        max_iter: int = 1000,
        tol: float = 1e-3,
    ):
        """Fit the WLS model to the data and return the deconvolution values.

        Parameters
        ----------
        adata : AnnData
            Anndata object containing the data.
        layer : str
            Layer of the data to use for deconvolution.
        ratio : bool, optional
            Transform raw output into cell proportions, by default False
        max_iter : int, optional
            Maximum number of iterations, by default 1000
        tol : float, optional
            Tolerance, by default 1e-3

        Returns
        -------
        pd.DataFrame
            Deconvolution values.
        """

        adata = self.valid_input(adata, layer=layer)
        values = solver_wls(adata, self._signature_matrix, layer=layer, max_iter=max_iter, tol=tol, dampened=None)

        values = self.format_output(
            values, columns=self._signature_matrix.list_cell_types, index=adata.obs.index, ratio=ratio
        )

        return values


class WNNLS(SignatureBasedModel):
    default_layer = "relative_counts"

    def __init__(self, signature_matrix: SignatureMatrix):
        """Weighted Non-Negative Least Squares (WNNLS) model for deconvolution.

        Parameters
        ----------
        signature_matrix : SignatureMatrix
            Signature matrix used for deconv.

        Raises
        ------
        ValueError
            Raised when solver_method is not found.
        """
        super().__init__(name="WNNLS", signature_matrix=signature_matrix)
        self._baseline_solver = SOLVER_METHODS.get("nnls")

    def transform(
        self,
        adata: AnnData,
        layer: str = "relative_counts",
        ratio: bool = False,
        max_iter: int = 1000,
        tol: float = 1e-3,
    ):
        """_summary_

        Parameters
        ----------
        adata : AnnData
            Input anndata object.
        layer : str
            Layer of the data to use for deconvolution.
        ratio : bool, optional
            if True return a proportion of celltype instead raw values, by default False
        max_iter : int, optional
            maximum number of iterations, by default 1000
        tol : _type_, optional
            tolerance, by default 1e-3

        Returns
        -------
        pd.DataFrame
            Deconvolution values.
        """

        adata = self.valid_input(adata, layer=layer)
        values = solver_wls(adata, self._signature_matrix, layer=layer, max_iter=max_iter, tol=tol, dampened=None)

        values = self.format_output(
            values, columns=self._signature_matrix.list_cell_types, index=adata.obs.index, ratio=ratio
        )

        return values


class DWLS(SignatureBasedModel):
    def __init__(self, signature_matrix: SignatureMatrix, solver_method: Literal["ols", "nnls", "svr"] = "nnls"):
        """Dumpened Weighted Least Squares (DWLS) model for deconvolution.

        Parameters
        ----------
        signature_matrix : SignatureMatrix
            Signature matrix used for deconv.
        solver_method : Literal[ols, nnls, svr], optional
            solver method to use, by default "nnls"

        Raises
        ------

        """
        super().__init__(name="DWLS", signature_matrix=signature_matrix)
        self._baseline_solver = SOLVER_METHODS.get(solver_method)
        if self._baseline_solver is None:
            raise ValueError(f"Solver {solver_method} not found")

    def transform(
        self,
        adata: AnnData,
        layer: str = "relative_counts",
        ratio: bool = False,
        max_iter: int = 1000,
        tol: float = 1e-2,
        dampened: Union[list, Literal["auto"], None] = "auto",
        parallel: bool = False,
    ):
        """Fit the DWLS model to the data and return the deconvolution values.

        Parameters
        ----------
        adata : AnnData
            Input anndata object.
        layer : str
            Layer of the data to use for deconvolution.
        ratio : bool, optional
            If True return a proportion of celltype instead raw values, by default False
        max_iter : int, optional
            Maximum number of iterations, by default 1000
        tol : float, optional
            Tolerance, by default 1
        dampened : Union[list, Literal[auto], None], optional
            Dampened values, by default "auto"
        parallel : bool, optional
            If True, use parallel computation, by default False

        Returns
        -------
        pd.DataFrame
            Deconvolution values.
        """
        if layer is None:
            layer = self.default_layer

        adata = self.valid_input(adata, layer=layer)
        values = solver_wls(
            adata,
            self._signature_matrix,
            layer=layer,
            max_iter=max_iter,
            tol=tol,
            dampened=dampened,
            solver_func=self._baseline_solver,
            parallel=parallel,
        )

        values = self.format_output(
            values, columns=self._signature_matrix.list_cell_types, index=adata.obs.index, ratio=ratio
        )
        return values
