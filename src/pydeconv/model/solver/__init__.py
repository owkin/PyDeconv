from .linear import solver_nnls, solver_ols, solver_rlr, solver_wls
from .nn import solver_mixupvi_module, solver_torch_module
from .svm import solver_nusvr, solver_svr


__all__ = [
    "solver_nnls",
    "solver_nusvr",
    "solver_ols",
    "solver_rlr",
    "solver_svr",
    "solver_torch_module",
    "solver_mixupvi_module",
    "solver_wls",
]
