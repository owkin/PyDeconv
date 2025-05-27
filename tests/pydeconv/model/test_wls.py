import numpy as np
import pytest
from anndata import AnnData

from pydeconv import SignatureMatrix
from pydeconv.model.wls import DWLS, WLS, WNNLS


@pytest.mark.parametrize("parallel", [False, True])
@pytest.mark.parametrize(
    "fake_bulkrnaseq_anndata, fake_signature_matrix",
    [
        ({"num_genes": 1000, "num_samples": 5}, {"num_genes": 1000, "num_cell_types": 10}),
        ({"num_genes": 1000, "num_samples": 5}, {"num_genes": 995, "num_cell_types": 10}),
    ],
    indirect=["fake_bulkrnaseq_anndata", "fake_signature_matrix"],
)
def test_dwls(fake_signature_matrix, fake_bulkrnaseq_anndata, parallel):
    signature_matrix = SignatureMatrix(fake_signature_matrix)
    solver = DWLS(signature_matrix)

    cell_scores = solver.transform(fake_bulkrnaseq_anndata, layer="bulk_rnaseq", ratio=False, parallel=parallel)

    assert cell_scores.shape == (len(fake_bulkrnaseq_anndata.obs), len(signature_matrix.list_cell_types))


# TODO: This text is failing because the results are not the same as the ones from the R package
# because when finding the dampening constant, the R package as not the same seed as the python
# package. We should find a way to fix this.
# def test_dwls_equivalence_r(dwls_mca_input, dwls_mca_results_dwls):
#     df_signature_matrix, df_bulk = dwls_mca_input
#     solver = DWLS(SignatureMatrix(df_signature_matrix))
#     adata = AnnData(df_bulk, layers={"tpm": df_bulk})

#     cell_scores = solver.transform(adata, layer="bulk_rnaseq", ratio=True, dampened="auto")
#     assert np.isclose(cell_scores, dwls_mca_results_dwls, atol=1e-4).all()


def test_dwls_equivalence_r_fixed_dampened(dwls_mca_input, dwls_mca_results_dwls, dwls_mca_dwls_dampened_values):
    df_signature_matrix, df_bulk = dwls_mca_input
    solver = DWLS(SignatureMatrix(df_signature_matrix))
    adata = AnnData(df_bulk, layers={"tpm": df_bulk})

    dampened = dwls_mca_dwls_dampened_values
    cell_scores = solver.transform(adata, layer="tpm", ratio=True, dampened=dampened)
    assert np.isclose(
        cell_scores, dwls_mca_results_dwls, atol=1e-4
    ).all()  # in R round 1e-5 but but 1e-4 to have exactly the same round


# TODO: very unstable because of the seed: need to find a way to fix this
# def test_dwls_find_dampening_constant_equivalence_r(dwls_mca_input, dwls_mca_dwls_dampened_values):
#     R_output = dwls_mca_dwls_dampened_values
#     df_signature_matrix, df_bulk = dwls_mca_input

#     signature_matrix = SignatureMatrix(df_signature_matrix)
#     adata = AnnData(df_bulk, layers={"tpm": df_bulk})

#     signature = signature_matrix.values  # row: gene, column: cell type
#     bulk = adata.to_df(layer="tpm").T  # row: genes, column: samples

#     cell_pred = solver_nnls(adata, signature_matrix, "tpm").T  # row: celltypes, column: samples

#     dampened = [find_dampening_constant(signature, b, c) for (_, b), c in zip(bulk.items(), cell_pred.T, strict=True)]

#     assert dampened == R_output


@pytest.mark.parametrize(
    "fake_bulkrnaseq_anndata, fake_signature_matrix",
    [
        ({"num_genes": 1000, "num_samples": 5}, {"num_genes": 1000, "num_cell_types": 10}),
        ({"num_genes": 1000, "num_samples": 5}, {"num_genes": 995, "num_cell_types": 10}),
    ],
    indirect=True,
)
def test_wls(fake_signature_matrix, fake_bulkrnaseq_anndata):
    signature_matrix = SignatureMatrix(fake_signature_matrix)
    solver = WLS(signature_matrix)

    cell_scores = solver.transform(fake_bulkrnaseq_anndata, layer="bulk_rnaseq", ratio=False)
    assert cell_scores.shape == (len(fake_bulkrnaseq_anndata.obs), len(signature_matrix.list_cell_types))


@pytest.mark.parametrize(
    "fake_bulkrnaseq_anndata, fake_signature_matrix",
    [
        ({"num_genes": 1000, "num_samples": 5}, {"num_genes": 1000, "num_cell_types": 10}),
        ({"num_genes": 1000, "num_samples": 5}, {"num_genes": 995, "num_cell_types": 10}),
    ],
    indirect=True,
)
def test_wnnls(fake_signature_matrix, fake_bulkrnaseq_anndata):
    signature_matrix = SignatureMatrix(fake_signature_matrix)
    solver = WNNLS(signature_matrix)

    cell_scores = solver.transform(fake_bulkrnaseq_anndata, layer="bulk_rnaseq", ratio=False)
    assert cell_scores.shape == (len(fake_bulkrnaseq_anndata.obs), len(signature_matrix.list_cell_types))
