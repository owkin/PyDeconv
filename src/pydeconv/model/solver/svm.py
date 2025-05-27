import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn import svm
from sklearn.preprocessing import StandardScaler

from pydeconv import SignatureMatrix


def solver_nusvr(
    adata: AnnData,
    signature_matrix: SignatureMatrix,
    layer: str,
    norm: bool = True,
    scale=False,
) -> np.ndarray:
    """Solve the linear regression using ordinary least squares

    Parameters
    ----------
    adata : AnnData
        Anndata object containing the data
    signature_matrix : SignatureMatrix
        Signature matrice object
    layer : str
        Layer to use for the regression
    norm : bool
        If True, the data will be normalized to the range [-1, 1]
    scale : bool
        If True, the data will be scaled to have zero mean and unit variance.

    Returns
    -------

    coef_ : np.ndarray
        Coefficients of the linear regression
    """
    model = svm.NuSVR(nu=0.5, gamma="auto", kernel="linear", C=1.0, cache_size=200)

    X = signature_matrix.values  # row: gene, column: cell type
    y = adata.to_df(layer=layer).T  # row: gene, column: samples

    # Normalize X and y to the range [-1, 1]
    if norm:
        # errited from R code
        X_max = np.max(X)
        X_min = np.min(X)
        upper_bounds = [np.max([X_max, y_max]) for y_max in y.max(axis=0).values]
        lower_bounds = [np.min([X_min, y_min]) for y_min in y.min(axis=0).values]

        X_norm = [
            ((X - lower_bound) / upper_bound) * 2 - 1
            for lower_bound, upper_bound in zip(lower_bounds, upper_bounds, strict=True)
        ]
        y_norm = [
            ((y[col] - lower_bound) / upper_bound) * 2 - 1
            for lower_bound, upper_bound, col in zip(lower_bounds, upper_bounds, y.columns, strict=True)
        ]

        if scale:
            X_norm_scale = [
                pd.DataFrame(StandardScaler().transform(X_norm_col), columns=X.columns, index=X.index)
                for X_norm_col in X_norm
            ]
            y_norm_scale = [StandardScaler().transform(y_norm_col.to_frame()) for y_norm_col in y_norm]
            coeffs = np.array(
                [
                    model.fit(X_norm_col, y_norm_col).coef_
                    for X_norm_col, y_norm_col in zip(X_norm_scale, y_norm_scale, strict=True)
                ]
            ).squeeze()
        else:
            coeffs = np.array(
                [model.fit(X_norm_col, y_norm_col).coef_ for X_norm_col, y_norm_col in zip(X_norm, y_norm, strict=True)]
            ).squeeze()

    else:
        if scale:
            X = pd.DataFrame(StandardScaler().transform(X), columns=X.columns, index=X.index)
            y = pd.DataFrame(StandardScaler().transform(y), columns=y.columns, index=y.index)

        coeffs = np.array([model.fit(X, y[col]).coef_ for col in y.columns]).squeeze()

    return coeffs


def solver_svr(adata: AnnData, signature_matrix: SignatureMatrix, layer: str) -> np.ndarray:
    """Solve the linear regression using ordinary least squares

    Parameters
    ----------

    adata : AnnData
        Anndata object containing the data
    signature_matrix : SignatureMatrix
        Signature matrice object
    layer : str
        Layer to use for the regression

    Returns
    -------

    coef_ : np.ndarray
        Coefficients of the linear regression
    """
    # model = make_pipeline(StandardScaler(), svm.NuSVR(nu=0.5, gamma="auto", kernel="linear", C=1.0, cache_size=40))
    model = svm.SVR(kernel="linear")

    X = signature_matrix.values  # row: gene, column: cell type
    y = adata.to_df(layer=layer).T  # row: gene, column: samples

    coeffs = np.array([model.fit(X, y[col]).coef_ for col in y.columns]).squeeze()

    return coeffs
