from anndata import AnnData

from pydeconv import SignatureMatrix
from pydeconv.utils import is_torch_available

from .linear import solver_nnls


if is_torch_available():
    import torch


def solver_torch_module(model: torch.nn.Module, adata: AnnData, layer: str, device: str = "cpu"):
    """Solve the deconvolution problem using a torch model.

    Parameters
    ----------
    model : torch.nn.Module
        The model to use for deconvolution.
    adata : AnnData
        Anndata object containing the data.
    layer : str
        Layer of the data to use for deconvolution.
    device : str, optional
        Device to use for the computation, by default "cpu"
    """
    model.eval()
    x = adata.to_df(layer=layer)
    x = torch.tensor(x.values, dtype=torch.float32).to(device)
    model = model.to(device)
    # can't be use with torch no grad context, because of optional adaptative step
    cell_prop_ = model(x)
    cell_prop_ = cell_prop_.cpu().detach().numpy()
    return cell_prop_


def solver_mixupvi_module(
    model: torch.nn.Module, signature_latent_space: SignatureMatrix, adata: AnnData, layer: str, device: str = "cpu"
):
    """Placeholder for MixupVI module solver."""
    x = adata.to_df(layer=layer)
    x = torch.tensor(x.values, dtype=torch.float32).to(device)
    model.eval()
    latent_features = model.extract_latent(x)
    anndata_latent = AnnData(latent_features, layers={"latent_features": latent_features})
    cell_prop_ = solver_nnls(anndata_latent, signature_latent_space, layer="latent_features")
    return cell_prop_
