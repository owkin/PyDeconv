from functools import partial
from typing import Callable

import numpy as np
from anndata import AnnData
from sklearn.preprocessing import MinMaxScaler

from .base_model import BasePreprocessing, NeuralNetworkModel
from .nn.registry import scaden_module_cti
from .solver import solver_torch_module


class Scaden(NeuralNetworkModel):
    model_registry: dict[str, Callable] = {
        "cti_1st_level_granularity": partial(scaden_module_cti, granularity="1st_granularity"),
        "cti_2nd_level_granularity": partial(scaden_module_cti, granularity="2nd_granularity"),
    }

    def __init__(self, weights_version: str = "cti_1st_level_granularity"):
        """Scaden model for deconvolution."""
        super().__init__(name="Scaden", weights_version=weights_version)
        self.preprocessing = PreprocessingScaden()

    def transform(self, adata: AnnData, layer: str = "counts_sum", ratio=False):
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
        adata_preprocessed = self.preprocessing(adata, layer, list_ordered_genes=self.params["input"]["gene_names"])
        values = solver_torch_module(self.module, adata_preprocessed, layer)

        values = self.format_output(
            values,
            columns=self.list_cell_types,
            index=adata.obs.index,
            ratio=ratio,
            only_positive=False,
        )

        return values


class PreprocessingScaden(BasePreprocessing):
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
        df = adata.to_df(layer)[list_ordered_genes]
        adata = AnnData(df, layers={layer: df})

        adata.layers[layer] = np.log2(adata.layers[layer] + 1)
        adata.layers[layer] = MinMaxScaler().fit_transform(adata.layers[layer].T).T
        return adata
