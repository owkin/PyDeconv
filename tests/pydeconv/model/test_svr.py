import numpy as np
import pytest
from anndata import AnnData

from pydeconv import SignatureMatrix
from pydeconv.model import SVR, NuSVR


@pytest.mark.parametrize(
    "fake_bulkrnaseq_anndata, fake_signature_matrix",
    [
        ({"num_genes": 1000, "num_samples": 5}, {"num_genes": 1000, "num_cell_types": 10}),
        ({"num_genes": 1000, "num_samples": 5}, {"num_genes": 995, "num_cell_types": 10}),
    ],
    indirect=True,
)
def test_nusvr(fake_signature_matrix, fake_bulkrnaseq_anndata):
    signature_matrix = SignatureMatrix(fake_signature_matrix)

    solver = NuSVR(signature_matrix)

    cell_scores = solver.transform(fake_bulkrnaseq_anndata, layer="bulk_rnaseq", ratio=False)
    assert cell_scores.shape == (len(fake_bulkrnaseq_anndata.obs), len(signature_matrix.list_cell_types))


def test_nusvr_equivalence_r_norm_no_scale(dwls_mca_input, dwls_mca_results_nusvr_norm_no_scale):
    df_signature_matrix, df_bulk = dwls_mca_input
    solver = NuSVR(SignatureMatrix(df_signature_matrix), norm=True)
    adata = AnnData(df_bulk, layers={"bulk_rnaseq": df_bulk})
    cell_scores = solver.transform(adata, layer="bulk_rnaseq", ratio=True)
    assert np.isclose(cell_scores, dwls_mca_results_nusvr_norm_no_scale, atol=1e-4).all()


def test_nusvr_equivalence_r_no_norm_no_scale(dwls_mca_input, dwls_mca_results_nusvr_no_norm_no_scale):
    df_signature_matrix, df_bulk = dwls_mca_input
    solver = NuSVR(SignatureMatrix(df_signature_matrix), norm=False)
    adata = AnnData(df_bulk, layers={"bulk_rnaseq": df_bulk})
    cell_scores = solver.transform(adata, layer="bulk_rnaseq", ratio=True)
    assert np.isclose(cell_scores, dwls_mca_results_nusvr_no_norm_no_scale, atol=1e-4).all()


# TODO: Fix the following tests

# def test_nusvr_equivalence_r(dwls_mca_input, dwls_mca_results_nusvr):
#     df_signature_matrix, df_bulk = dwls_mca_input
#     solver = NuSVR(SignatureMatrix(df_signature_matrix), norm=True, scale=True)
#     adata = AnnData(df_bulk, layers={"bulk_rnaseq": df_bulk})
#     cell_scores = solver.transform(adata, layer="bulk_rnaseq", ratio=True)
#     assert np.isclose(cell_scores, dwls_mca_results_nusvr, atol=1e-4).all()


# def test_nusvr_equivalence_r_no_norm_scale(dwls_mca_input, dwls_mca_results_nusvr_no_norm_scale):
#     df_signature_matrix, df_bulk = dwls_mca_input
#     solver = NuSVR(SignatureMatrix(df_signature_matrix), norm=False, scale=True)
#     adata = AnnData(df_bulk, layers={"bulk_rnaseq": df_bulk})
#     cell_scores = solver.transform(adata, layer="bulk_rnaseq", ratio=True)
#     assert np.isclose(cell_scores, dwls_mca_results_nusvr_no_norm_scale, atol=1e-4).all()


@pytest.mark.parametrize(
    "fake_bulkrnaseq_anndata, fake_signature_matrix",
    [
        ({"num_genes": 1000, "num_samples": 5}, {"num_genes": 1000, "num_cell_types": 10}),
    ],
    indirect=True,
)
def test_svr(fake_signature_matrix, fake_bulkrnaseq_anndata):
    signature_matrix = SignatureMatrix(fake_signature_matrix)

    solver = SVR(signature_matrix)

    cell_scores = solver.transform(fake_bulkrnaseq_anndata, layer="bulk_rnaseq", ratio=False)
    assert cell_scores.shape == (len(fake_bulkrnaseq_anndata.obs), len(signature_matrix.list_cell_types))
