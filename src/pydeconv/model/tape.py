from functools import partial
from typing import Callable

import numpy as np
from anndata import AnnData
from sklearn.preprocessing import MinMaxScaler

from .base_model import BasePreprocessing, NeuralNetworkModel
from .nn.registry import tape_module_cti
from .solver import solver_torch_module


class Tape(NeuralNetworkModel):
    """Tape model for deconvolution.

    Attributes
    ----------
    name : str
        Name of the model.

    Methods
    -------
    transform(adata, layer, ratio=False)
        Fit the Tape model to the data and return the deconvolution values.
    """

    model_registry: dict[str, Callable] = {
        "cti_1st_level_granularity": partial(tape_module_cti, granularity="1st_granularity"),
        "cti_2nd_level_granularity": partial(tape_module_cti, granularity="2nd_granularity"),
    }

    def __init__(self, weights_version: str = "cti_1st_level_granularity", adaptative=True):
        """

        Parameters
        ----------
        adaptative : bool, optional
            _description_, by default True
        """
        super().__init__(name="TAPE", weights_version=weights_version, adaptative=adaptative)
        self.preprocessing = PreprocessingTape()

    def transform(self, adata: AnnData, layer: str = "counts_sum", ratio=False):
        """Fit the Tape model to the data and return the deconvolution values.

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
            only_positive=True,
        )

        return values


class PreprocessingTape(BasePreprocessing):
    def preprocess(self, adata: AnnData, layer: str, list_ordered_genes: list):
        """Preprocess the data for Tape model

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
        """
        df = adata.to_df(layer)[list_ordered_genes]
        adata = AnnData(df, layers={layer: df})

        adata.layers[layer] = np.log2(adata.layers[layer] + 1)
        adata.layers[layer] = MinMaxScaler().fit_transform(adata.layers[layer].T).T
        return adata
