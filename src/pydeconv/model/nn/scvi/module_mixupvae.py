"""MixUpVAE module."""

from pydeconv.utils import is_scvi_available, requires_scvi


if is_scvi_available():
    import torch
    from scvi.module import VAE
    from scvi.module.base import auto_move_data

    torch.backends.cudnn.benchmark = True


@requires_scvi
class _MixUpVAE(VAE):
    """Variational auto-encoder model with linearity constraint within batches.

    The linearity constraint is inspired by the MixUp method
    (https://arxiv.org/abs/1710.09412v2).

    Parameters
    ----------
    n_input
        Number of input genes
    n_latent
        Dimensionality of the latent space
    seed
        The desired seed.
    library_log_means
        1 x n_batch array of means of the log library sizes. Parameterizes prior on library size if
        not using observed library size.
    library_log_vars
        1 x n_batch array of variances of the log library sizes. Parameterizes prior on library size if
        not using observed library size.
    """

    def __init__(
        self,
        # VAE arguments
        n_input: int,
        n_latent: int = 10,
        seed: int = 0,
    ):
        torch.manual_seed(seed)

        super().__init__(
            n_input=n_input,
            n_batch=1,
            n_labels=0,
            n_hidden=512,
            n_latent=n_latent,
            n_layers=1,
            n_continuous_cov=0,
            n_cats_per_cov=None,
            dropout_rate=0.1,
            dispersion="gene",
            log_variational=True,
            gene_likelihood="zinb",
            latent_distribution="normal",
            encode_covariates=False,
            deeply_inject_covariates=True,
            use_batch_norm="none",
            use_layer_norm="none",
            use_size_factor_key=False,
            use_observed_lib_size=True,
            library_log_means=None,
            library_log_vars=None,
            var_activation=None,
            extra_encoder_kwargs=None,
            extra_decoder_kwargs=None,
        )

    @auto_move_data
    def inference(
        self,
        x,
    ):
        """High level inference method for pseudobulks of single cells.

        Runs the inference (encoder) model.
        """
        x_ = x
        encoder_input = torch.log(1 + x_)

        # regular encoding
        qz, z = self.z_encoder(
            encoder_input,
            torch.ones((x_.shape[0], 1)),
        )

        outputs = {
            # encodings
            "z": z,
            "qz": qz,
        }

        return outputs


@torch.inference_mode()
def get_latent_representation(
    x: torch.Tensor,
    module: VAE,
):
    """Return the latent representation for each cell or pseudobulk.

    This is typically denoted as :math:`z_n`.

    Parameters
    ----------
    x: torch.Tensor
        The single-cells or pseudobulks to infer from.
    module: scvi.module.VAE
        The scVI module to infer from, typically MixUpVAE.

    Returns
    -------
    Low-dimensional representation for each cell or a tuple containing its mean and variance.
    """
    latent = []
    outputs = module.inference(x)
    qz = outputs["qz"]
    z = qz.loc
    latent += [z.cpu()]
    return torch.cat(latent).numpy()
