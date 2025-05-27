import pytest

from pydeconv import SignatureMatrix
from pydeconv.model import OLS


@pytest.mark.parametrize(
    "fake_bulkrnaseq_anndata, fake_signature_matrix",
    [
        ({"num_genes": 1000, "num_samples": 5}, {"num_genes": 995, "num_cell_types": 10}),  # with more input genes
        ({"num_genes": 995, "num_samples": 5}, {"num_genes": 1000, "num_cell_types": 10}),  # with missing input genes
    ],
    indirect=True,
)
def test_ols(fake_signature_matrix, fake_bulkrnaseq_anndata):
    signature_matrix = SignatureMatrix(fake_signature_matrix)
    solver = OLS(signature_matrix)

    cell_scores = solver.transform(fake_bulkrnaseq_anndata, layer="bulk_rnaseq", ratio=False)
    assert cell_scores.shape == (len(fake_bulkrnaseq_anndata.obs), len(signature_matrix.list_cell_types))
