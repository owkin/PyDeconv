from abc import ABC, abstractmethod
from typing import Callable

import numpy as np
import pandas as pd
from anndata import AnnData

from pydeconv.utils import valid_anndata


class BaseModel(ABC):
    def __init__(self, name: str):
        """Base model for deconvolution.

        Parameters
        ----------
        name : str
            Name of the model.
        signature_matrix : SignatureMatrix, optional
            Signature matrix used for deconvolution, by default None
        """
        self._name = name

    @abstractmethod
    def transform(self, adata: AnnData, layer: str, ratio=False):
        """Fit the model and return the deconvolution values.

        Parameters
        ----------
        adata : AnnData
            input anndata object.
        layer : str
            layer of the data to use for deconvolution.
        ratio : bool, optional
            If True return a proportion of celltype instead raw values, by default False
        """
        pass

    @staticmethod
    def format_output(
        values: np.ndarray, columns: list, index: list, ratio: bool, only_positive: bool = False
    ) -> pd.DataFrame:
        """Format the output of the model.

        Parameters
        ----------
        values : np.ndarray
            values to format

        columns : list
            name of the columns corresponding to the cell types
        index : list
            name of the index corresponding to the patients
        ratio : bool
            If True return a proportion of celltype instead raw values
        only_positive : bool, optional
            If True clip the values to be positive, by default False

        Returns
        -------
        pd.DataFrame
            Formatted values
        """
        if only_positive:
            values = values.clip(min=0)

        if ratio:
            values = values / values.sum(axis=1)[:, None]

        values = pd.DataFrame(values, columns=columns, index=index)

        return values

    @abstractmethod
    def valid_input(self, input_anndata: AnnData, layer="tpm", tol=0.3):
        pass


class SignatureBasedModel(BaseModel):
    def __init__(self, name: str, signature_matrix):
        """Base model for signature-based deconvolution.

        Parameters
        ----------
        name : str
            Name of the model.
        signature_matrix : SignatureMatrix
            Signature matrix used for deconvolution.
        """
        super().__init__(name)
        self._signature_matrix = signature_matrix
        self.list_cell_types = self._signature_matrix.list_cell_types

    def valid_input(self, input_anndata: AnnData, layer="tpm", tol=0.3):
        """Check if the input anndata is valid.

        Parameters
        ----------
        input_anndata : AnnData
            data to check
        layer : str, optional
            layer of the data to check, by default "tpm"
        tol : float, optional
            tolerance for the gene names, by default 0.3

        Raises
        ------
        ValueError
            If the input anndata is not valid.
        """

        adata, new_gene_list = valid_anndata(input_anndata, self._signature_matrix.list_genes, tol=tol)
        if new_gene_list is not None:
            self._signature_matrix.update_gene_list(new_gene_list)
        return adata


class NeuralNetworkModel(BaseModel):
    model_registry: dict[str, Callable] = {}

    def __init__(self, name: str, weights_version: str, **kwargs):
        """Base model for neural network-based deconvolution.

        Parameters
        ----------
        name : str
            Name of the model.
        params : dict
            Parameters of the model.
        """
        super().__init__(name)
        if weights_version not in self.model_registry:
            raise ValueError(
                f"Invalid weights version: {weights_version}. Available versions: {list(self.model_registry.keys())}"
            )
        self.module, self.params = self.model_registry[weights_version](**kwargs)
        self.list_cell_types = self.params["output"]["cell_types"]

    def valid_input(self, input_anndata: AnnData, layer="tpm", tol=0.3):
        """Check if the input anndata is valid.

        Parameters
        ----------
        input_anndata : AnnData
            data to check
        layer : str, optional
            layer of the data to check, by default "tpm"
        tol : float, optional
            tolerance for the gene names, by default 0.3

        Raises
        ------
        ValueError
            If the input anndata is not valid.
        """

        adata, new_gene_list = valid_anndata(input_anndata, self.params["input"]["gene_names"], tol=tol)
        if new_gene_list is not None:
            # we can't edit the architecture of neural network so we fill the adata with 0 instead
            empty_genes = list(set(self.params["input"]["gene_names"]) - set(new_gene_list))
            df = adata.to_df(layer=layer)
            df[empty_genes] = 0.0
            df = df[self.params["input"]["gene_names"]]
            adata = AnnData(df, layers={layer: df})

        return adata


class BasePreprocessing(ABC):
    @abstractmethod
    def preprocess(self, adata: AnnData, layer: str, list_ordered_genes: list):
        """Base class for preprocessing.

        Parameters
        ----------
        adata : AnnData
            _description_
        layer : str
            _description_
        list_ordered_genes : list
            _description_
        """
        pass

    def __call__(self, adata: AnnData, layer: str, list_ordered_genes, inplace=False):
        """Apply the preprocessing.

        Parameters
        ----------
        adata : AnnData
            Input anndata
        layer : str
            Anndata layer to use.
        list_ordered_genes : list
            List of genes to use.
        inplace : bool, optional
            If yes apply the modification inplace, by default False

        Returns
        -------
        AnnData
            Preprocessed anndata
        """
        if inplace:
            output = adata
        else:
            output = adata.copy()
        return self.preprocess(output, layer, list_ordered_genes)
