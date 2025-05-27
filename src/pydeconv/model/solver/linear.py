import logging
from typing import Literal, Union

import numpy as np
import pandas as pd
import statsmodels.api as sm
from anndata import AnnData
from quadprog import solve_qp
from scipy import linalg, optimize
from tqdm import tqdm

from pydeconv import SignatureMatrix


def solver_ols(adata: AnnData, signature_matrix: SignatureMatrix, layer: str):
    """
    Solve the linear regression using ordinary least squares.

    Parameters
    ----------
    adata : AnnData
        Anndata object containing the data.
    signature_matrix : SignatureMatrix
        Signature matrix object.
    layer : str
        Layer to use for the regression.

    Returns
    -------
    coef_ : np.ndarray
        Coefficients of the linear regression.
    """
    X = signature_matrix.values
    Y = adata.to_df(layer=layer).T

    coef_, _, _, _ = linalg.lstsq(X, Y)
    coef_ = coef_.T

    return coef_


def solver_nnls(adata: AnnData, signature_matrix: SignatureMatrix, layer: str):
    """
    Solve the linear regression using non-negative least squares.

    Parameters
    ----------
    adata : AnnData
        Anndata object containing the data.
    signature_matrix : SignatureMatrix
        Signature matrix object.
    layer : str
        Layer to use for the regressions.

    Returns
    -------
    coef_ : np.ndarray
        Coefficients of the linear regression.
    """
    X = signature_matrix.values
    Y = adata.to_df(layer=layer).T

    coef_ = np.vstack([optimize.nnls(X, yi)[0] for yi in Y.T.values])
    return coef_


def solver_rlr(adata: AnnData, signature_matrix: SignatureMatrix, layer: str):
    """
    Solve the linear regression using robust linear regression. and huber loss function.

    Parameters
    ----------
    adata : AnnData
        Anndata object containing the data.
    signature_matrix : SignatureMatrix
        Signature matrix object.
    layer : str
        Layer to use for the regression.

    Returns
    -------
    coef_ : np.ndarray
        Coefficients of the linear regression.
    """
    X = signature_matrix.values
    Y = adata.to_df(layer=layer).T

    coef_ = np.vstack([sm.RLM(y, X, M=sm.robust.norms.HuberT()).fit().params for y in Y.T.values])
    return coef_


def solver_wls(
    adata: AnnData,
    signature_matrix: SignatureMatrix,
    layer: str,
    max_iter: int = 1000,
    tol: float = 0.01,
    dampened: Union[list, Literal["auto"], None] = None,
    solver_func=solver_nnls,
    parallel: bool = False,
):
    """
    Solve the linear regression using Weighted Least Squares.

    Parameters
    ----------
    adata : AnnData
        Anndata object containing the data.
    signature_matrix : SignatureMatrix
        Signature matrix object.
    layer : str
        Layer to use for the regression.
    max_iter : int, optional
        Maximum number of iterations, by default 1000.
    tol : float, optional
        Tolerance for the convergence, by default 0.01.
    dampened : Union[list, Literal["auto"], None], optional
        Dampening constant for the weights. If None, no dampening is applied.
        If "auto", the dampening constant is calculated automatically, by default None.
    solver_func : Callable, optional
        Function to solve the linear regression, by default solver_nnls.
    parallel : bool, optional
        Whether to run the solver in parallel, by default False.

    Returns
    -------
    coef_ : np.ndarray
        Coefficients of the linear regression.
    """

    X = signature_matrix.values.copy()  # row: gene, column: cell type
    Y = adata.to_df(layer=layer).T.copy()  # row: genes, column: samples

    coef_init = solver_func(adata, signature_matrix, layer).T  # row: celltypes, column: samples

    # Define dampening constants
    if dampened == "auto":
        dampened = [find_dampening_constant(X, y, coef) for (_, y), coef in zip(Y.items(), coef_init.T, strict=True)]
    elif isinstance(dampened, float):
        dampened = np.array([dampened] * Y.shape[1])
    elif dampened is None:
        dampened = [None] * Y.shape[1]
    if len(dampened) != Y.shape[1]:
        raise ValueError("Dampening constant must be a float or a list of floats with the same length as the samples")

    # Solve the weighted linear regression
    if parallel:
        coef_ = solver_wls_parallel(X, Y, coef_init, max_iter=max_iter, tol=tol, dampened=dampened)
    else:
        coef_ = [
            solver_wls_per_patient(X, y, c, d, max_iter, tol)
            for (_, y), c, d in tqdm(
                zip(Y.items(), coef_init.T, dampened, strict=True), total=len(Y), position=3, leave=False
            )
        ]
    coef_ = np.array(coef_)

    return coef_  # output: row: samples, column: cell types


def solver_wls_per_patient(
    X: pd.DataFrame,
    y: np.ndarray,
    cell_pred: np.ndarray,
    dampened: Union[None, float],
    max_iter: int,
    tol: float,
):
    """
    Solve the linear regression using Weighted Least Squares for a single patient.

    Parameters
    ----------
    X : pd.DataFrame
        Signature matrix.
    y : np.ndarray
        Bulk data.
    cell_pred : np.ndarray
        Initial cell predictions.
    dampened : Union[None, float]
        Dampening constant for the weights. If None, no dampening is applied.
    max_iter : int
        Maximum number of iterations.
    tol : float
        Tolerance for the convergence.

    Returns
    -------
    cell_pred_new : np.ndarray
        Coefficients of the linear regression.
    """
    for iteration in range(max_iter):
        # solver dampeled weighted least squares j

        y_pred = X @ cell_pred

        y_pred = y_pred.clip(1e-4)  # Avoid division by zero (<1e-8) for y_pred**2
        weights = 1 / (y_pred**2)

        if dampened is not None:
            weights = weights / weights.min()
            weights = np.minimum(weights, dampened)

        cell_pred_new = solver_nnls_norm(X.values, y, weights)
        cell_pred_new = np.array(cell_pred_new)

        # reduce step for stability
        cell_pred_new = reduce_step(cell_pred, cell_pred_new)

        if np.linalg.norm(cell_pred_new - cell_pred, ord=1) <= tol:  # R implementation uses L1 norm so we do the same
            logging.info(f"Converged after {iteration+1} iterations")
            break
        else:
            cell_pred = cell_pred_new
        if iteration == max_iter - 1:
            logging.warning(f"Did not converge after {max_iter} iterations")
    return cell_pred


def solver_nnls_norm(X, y, wsDampened=None):
    """
    Solve the linear regression using non-negative least squares.

    Parameters
    ----------
    X : np.ndarray
        Signature matrix.
    y : np.ndarray
        Bulk matrix.
    wsDampened : np.ndarray, optional
        Dampened weights, by default None.

    Returns
    -------
    solution : np.ndarray
        Coefficients of the linear regression.
    """

    W = np.diag(wsDampened)

    # Step 2: Compute the weighted matrix D and vector d
    D = X.T @ W @ X
    d = X.T @ W @ y

    # Step 3: Normalize the matrix D and vector d
    sc = np.linalg.norm(D, 2)
    D = D / sc
    d = d / sc

    # Step 4: Set up the constraints for quadratic programming
    A = np.eye(X.shape[1])  # Identity matrix for non-negativity constraints
    bzero = np.zeros(X.shape[1])  # Vector of zeros for the constraints

    solution = solve_qp(D, d, A, bzero)[0]
    return solution


def reduce_step(solution, new_solution, method: str = "weighted_average", w=4):
    """
    Reduce the step size for the new solution.

    Parameters
    ----------
    solution : np.ndarray
        Original solution.
    new_solution : np.ndarray
        New solution.
    method : str, optional
        Method to reduce the step size, by default "weighted_average".
    w : int, optional
        Weight for the original solution, by default 4.

    Returns
    -------
    reduced_solution : np.ndarray
        Reduced solution.
    """
    if method == "weighted_average":
        # Weighted average of the two solutions with the original solution having a weight of 4
        solution_tiled = np.array([solution for _ in range(w)])
        combined = np.concatenate((new_solution[np.newaxis, :], solution_tiled))
        return combined.mean(axis=0)
    elif method == "average":
        return (solution + new_solution) / 2
    else:
        raise ValueError(f"Method {method} not found")


def solver_wls_parallel(
    X: pd.DataFrame,
    Y: pd.DataFrame,
    cell_pred: np.ndarray,
    dampened: list,
    max_iter: int,
    tol: float,
):
    """
    Solve the linear regression using Weighted Least Squares in parallel.

    Parameters
    ----------
    X : pd.DataFrame
        Signature matrix.
    Y : pd.DataFrame
        Bulk data.
    cell_pred : np.ndarray
        Initial cell predictions.
    dampened : list
        Dampening constant for the weights.
    max_iter : int
        Maximum number of iterations.
    tol : float
        Tolerance for the convergence.

    Returns
    -------
    np.ndarray
        Coefficients of the linear regression.
    """

    for iteration in tqdm(range(max_iter), desc="DWLS", unit="iterations"):
        # solver dampeled weighted least squares j

        bulk_pred = X @ cell_pred
        bulk_pred.columns = Y.columns

        bulk_pred = bulk_pred.clip(1e-4)  # Avoid division by zero (<1e-8) for y_pred**2
        weights = 1 / (bulk_pred**2)

        if dampened is not None:
            weights = weights / weights.min()
            weights = np.minimum(weights, dampened)

        cell_pred_new = np.array(
            [
                solver_nnls_norm(X.values, y.values, w.values)
                for (_, y), (_, w) in zip(Y.items(), weights.items(), strict=True)
            ]
        ).T

        # reduce step for stability
        cell_pred_new = reduce_step(cell_pred, cell_pred_new)

        if np.linalg.norm(cell_pred_new - cell_pred, ord=1) < tol:  # R implementation uses L1 norm so we do the same
            logging.info(f"Converged after {iteration+1} iterations")
            break
        else:
            cell_pred = cell_pred_new
        if iteration == max_iter - 1:
            logging.warning(f"Did not converge after {max_iter} iterations")
    return cell_pred.T  # output: row: samples, column: cell types


def find_dampening_constant(X: pd.DataFrame, y: pd.DataFrame, init_solution: np.ndarray):
    """
    Find the dampening constant for the weighted least squares regression.
    The dampening constant is chosen such that the cross-validation variance is minimized.

    Parameters
    ----------
    X : pd.DataFrame
        Signature matrix.
    y : pd.DataFrame
        Bulk data.
    init_solution : np.ndarray
        Initial cell predictions.

    Returns
    -------
    dampening_constant : int
        Dampening constant for the weighted least squares regression.
    """
    ws = (1 / (X @ init_solution)) ** 2
    ws_scaled = ws / np.min(ws)

    solutions_sd = []

    # Try multiple values of the dampening constant (multiplier)
    for j in range(int(np.ceil(np.log2(np.max(ws_scaled[np.isfinite(ws_scaled)]))))):
        multiplier = 2**j
        ws_dampened = np.minimum(ws_scaled, multiplier)

        solutions = []

        for i in range(100):
            np.random.seed(i)
            subset = np.random.choice(len(ws), size=int(len(ws) * 0.5), replace=False)

            # Solve dampened weighted least squares for subset
            model = sm.WLS(y.iloc[subset], X.iloc[subset], weights=ws_dampened.iloc[subset])
            fit = model.fit()
            sol = fit.params * np.sum(init_solution) / float(fit.params.sum())

            solutions.append(sol)

        solutions_sd.append(np.std(solutions, axis=0))

    # Choose dampening constant that results in least cross-validation variance
    mean_sd_solution = np.mean(np.array(solutions_sd).squeeze() ** 2, axis=1)
    j = int(np.argmin(mean_sd_solution))
    return 2**j
