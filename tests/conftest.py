"""Fixtures module for the tests."""

from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr
import pytest
from anndata import AnnData


@pytest.fixture()
def fake_bulkrnaseq_anndata(request):
    num_genes = request.param.get("num_genes", 1000)
    num_samples = request.param.get("num_samples", 100)
    gene_names = request.param.get("gene_names", None)

    np.random.seed(42)
    if gene_names is None:
        gene_names = [f"Gene_{i}" for i in range(1, num_genes + 1)]
    sample_labels = [f"Sample_{j}" for j in range(1, num_samples + 1)]
    expression_matrix = np.random.randint(0, 1000, size=(num_samples, len(gene_names)))
    expression_df = pd.DataFrame(expression_matrix, index=sample_labels, columns=gene_names)
    adata = AnnData(expression_df)
    adata.layers["bulk_rnaseq"] = expression_df
    return adata


@pytest.fixture()
def fake_signature_matrix(request):
    num_genes = request.param.get("num_genes", 1000)
    num_cell_types = request.param.get("num_cell_types", 10)

    np.random.seed(42)
    signature_matrix = np.random.rand(num_genes, num_cell_types)
    gene_names = [f"Gene_{i}" for i in range(1, num_genes + 1)]
    cell_type_labels = [f"CellType_{j}" for j in range(1, num_cell_types + 1)]
    signature_df = pd.DataFrame(signature_matrix, index=gene_names, columns=cell_type_labels)
    return signature_df


def dwls_mca_signature_matrix():
    df_signature_matrix = pyreadr.read_r(Path(__file__).parent / "data/dwls/mca/Sig.RData")["Sig"]
    return df_signature_matrix


def dwls_mca_bulk():
    names = pyreadr.read_r(Path(__file__).parent / "data/dwls/mca/list_genes.RData")["names"]
    list_df = [
        pyreadr.read_r(Path(__file__).parent / f"data/dwls/mca/bulkData{i}.RData")[f"bulkData{i}"]
        for i in [1, 2, 3, 4, 9, 10, 11, 12]
    ]
    df = pd.concat(list_df, axis=1)
    df.index = names.names

    return df


@pytest.fixture()
def dwls_mca_input():
    df_signature_matrix = dwls_mca_signature_matrix()
    df_bulk = dwls_mca_bulk()

    genes = df_signature_matrix.index.intersection(df_bulk.index)

    df_signature_matrix = df_signature_matrix.loc[genes]
    df_bulk = df_bulk.loc[genes]
    df_bulk = df_bulk[~df_bulk.index.duplicated(keep="first")]
    df_bulk = df_bulk.T  # Transpose to have samples as rows

    return df_signature_matrix, df_bulk


@pytest.fixture()
def dwls_mca_dwls_dampened_values():
    # not stable:
    # generated from R with seed c(100:200)
    # return [36, 38, 15, 47, 14, 27, 14, 27]
    # generated from R with seed c(500:600)
    # return [36, 41, 34, 41, 14, 27, 14, 28]

    # generated from R with seed c(0:100)
    j = [16, 52, 15, 46, 14, 27, 14, 27]
    j = [2 ** (val - 1) for val in j]
    return j


def dwls_mca_results(method, options=""):
    df_results = pyreadr.read_r(Path(__file__).parent / f"data/dwls/mca/allCounts_{method}{options}.RData")[
        f"allCounts_{method}"
    ]
    df_results.columns = [f"bulkData{i}" for i in [1, 2, 3, 4, 9, 10, 11, 12]]
    df_results = df_results.T
    return df_results


@pytest.fixture()
def dwls_mca_results_nnls():
    df_results = dwls_mca_results("OLS")  # named OLS but actually NNLS
    return df_results


@pytest.fixture()
def dwls_mca_results_nusvr():
    df_results = dwls_mca_results("SVR")
    return df_results


@pytest.fixture()
def dwls_mca_results_nusvr_norm_no_scale():
    df_results = dwls_mca_results("SVR", "_norm_no_scale")
    return df_results


@pytest.fixture()
def dwls_mca_results_nusvr_no_norm_no_scale():
    df_results = dwls_mca_results("SVR", "_no_norm_no_scale")
    return df_results


@pytest.fixture()
def dwls_mca_results_nusvr_no_norm_scale():
    df_results = dwls_mca_results("SVR", "_no_norm_scale")
    return df_results


@pytest.fixture()
def dwls_mca_results_dwls():
    df_results = dwls_mca_results("DWLS")
    return df_results


def tape_input(method: str):
    path_folder = Path(__file__).parent / "data/tape"
    df_input = pd.read_csv(path_folder / f"pseudobulk_test_{method}.csv", index_col=0)
    return df_input


@pytest.fixture()
def tape_input_tape():
    return tape_input("tape")


@pytest.fixture()
def tape_input_scaden():
    return tape_input("scaden")


def tape_results(method: str, suffix: str = ""):
    path_folder = Path(__file__).parent / "data/tape"
    df_results = pd.read_csv(path_folder / f"deconv_results_{method}{suffix}.csv", index_col=0)
    return df_results


@pytest.fixture()
def tape_results_tape():
    return tape_results("tape", "_no_adaptative")


@pytest.fixture()
def tape_results_scaden():
    return tape_results("scaden")
