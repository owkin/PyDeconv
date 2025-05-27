# for deconvolution

from pydeconv.utils import is_scvi_available, requires_scvi


if is_scvi_available():
    import torch
    import torch.nn as nn

    torch.backends.cudnn.benchmark = True
    from pydeconv.model.nn.scvi.module_mixupvae import get_latent_representation

    from .scvi.module_mixupvae import _MixUpVAE


@requires_scvi
class MixupVIModule(nn.Module):
    """MixupVAE is a simple implementation of the Mixup Variational Autoencoder model.
    paper: https://arxiv.org/abs/1803.11096
    """

    def __init__(self, input_dim: int = 3000, latent_dim: int = 30):
        """Initialize the MixupVAE model.

        Parameters
        ----------
        input_dim : int
            Input dimension.
        latent_dim : int
            Latent dimension.
        """
        super().__init__()
        self.model = _MixUpVAE(
            n_input=input_dim,
            n_latent=latent_dim,
        )

    def extract_latent(self, x):
        """Forward pass of the MixupVAE model.

        Parameters
        ----------
        x : torch.Tensor
            Input tensor.

        Returns
        -------
        torch.Tensor
            Output tensor.
        """
        # This is a placeholder for the actual implementation
        self.model.n_pseudobulks = x.shape[0]
        latent_adata = get_latent_representation(x, self.model)
        return latent_adata

    def eval(self):
        """Set the model to evaluation mode."""
        self.model.eval()
        return self
