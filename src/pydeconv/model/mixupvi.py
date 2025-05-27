from functools import partial
from typing import Callable, Optional

from anndata import AnnData

from pydeconv.signature_matrix.registry import (
    signature_mixupvi_latent_space_1st_granularity,
    signature_mixupvi_latent_space_2nd_granularity,
)

from .base_model import BasePreprocessing, NeuralNetworkModel
from .nn.registry import mixupvi_module_cti
from .solver import solver_mixupvi_module


class MixupVI(NeuralNetworkModel):
    model_registry: dict[str, Callable] = {
        "cti_1st_level_granularity": partial(mixupvi_module_cti, granularity="1st_granularity"),
        "cti_2nd_level_granularity": partial(mixupvi_module_cti, granularity="2nd_granularity"),
    }
    signature_registry: dict[str, Callable] = {
        "cti_1st_level_granularity": signature_mixupvi_latent_space_1st_granularity,
        "cti_2nd_level_granularity": signature_mixupvi_latent_space_2nd_granularity,
    }

    def __init__(
        self, weights_version: str = "cti_1st_level_granularity", preprocessing: Optional[BasePreprocessing] = None
    ):
        """MixupVI model for deconvolution."""
        super().__init__(name="MixupVI", weights_version=weights_version)
        self.preprocessing = preprocessing
        self.signature_latent_space = self.signature_registry[weights_version]()

    def transform(self, adata: AnnData, layer: str = "raw_counts", ratio=False):
        """Fit the Scaden model to the data and return the deconvolution values.

        Parameters
        ----------
        adata : AnnData
            Anndata object containing the data.
        layer : str
            Layer of the data to use for deconvolution.
        ratio : bool, optional
            Transform raw output into cell proportions, by default False

        Returns
        -------
        pd.DataFrame
            Deconvolution values.
        """

        adata = self.valid_input(adata, layer=layer)
        if self.preprocessing is not None:
            adata_preprocessed = self.preprocessing(adata, layer, list_ordered_genes=self.params["input"]["gene_names"])
        else:
            adata_preprocessed = adata
        values = solver_mixupvi_module(self.module, self.signature_latent_space, adata_preprocessed, layer)

        values = self.format_output(
            values,
            columns=self.list_cell_types,
            index=adata.obs.index,
            ratio=ratio,
            only_positive=False,
        )

        return values


class PreprocessingMixupVI(BasePreprocessing):
    def preprocess(self, adata: AnnData, layer: str, list_ordered_genes: list):
        """Preprocess the data for Scaden model

        Parameters
        ----------
        adata : AnnData
            input anndata object.
        layer : str
            layer of the data to use for deconvolution.
        list_ordered_genes : list
            list of genes to use for deconvolution.

        Returns
        -------
        AnnData
            Preprocessed data
        """
        return adata
