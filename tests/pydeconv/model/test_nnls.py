import numpy as np
import pytest
from anndata import AnnData

from pydeconv import SignatureMatrix
from pydeconv.model import NNLS


@pytest.mark.parametrize(
    "fake_bulkrnaseq_anndata, fake_signature_matrix",
    [
        ({"num_genes": 1000, "num_samples": 5}, {"num_genes": 1000, "num_cell_types": 10}),
        ({"num_genes": 1000, "num_samples": 5}, {"num_genes": 995, "num_cell_types": 10}),
    ],
    indirect=True,
)
def test_nnls(fake_signature_matrix, fake_bulkrnaseq_anndata):
    signature_matrix = SignatureMatrix(fake_signature_matrix)
    solver = NNLS(signature_matrix)

    cell_scores = solver.transform(fake_bulkrnaseq_anndata, layer="bulk_rnaseq", ratio=False)

    assert (cell_scores > 0).all(None)
    assert cell_scores.shape == (len(fake_bulkrnaseq_anndata.obs), len(signature_matrix.list_cell_types))


def test_nnls_equivalence_r(dwls_mca_input, dwls_mca_results_nnls):
    df_signature_matrix, df_bulk = dwls_mca_input
    solver = NNLS(SignatureMatrix(df_signature_matrix))
    adata = AnnData(df_bulk, layers={"bulk_rnaseq": df_bulk})
    cell_scores = solver.transform(adata, layer="bulk_rnaseq", ratio=True)
    assert np.isclose(cell_scores, dwls_mca_results_nnls, atol=1e-4).all()
