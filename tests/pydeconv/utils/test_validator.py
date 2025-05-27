import pytest

from pydeconv.signature_matrix import SignatureMatrix
from pydeconv.utils import valid_anndata


@pytest.mark.parametrize(
    "fake_bulkrnaseq_anndata, fake_signature_matrix",
    [
        ({"num_genes": 1000, "num_samples": 100}, {"num_genes": 1000, "num_cell_types": 10}),
        ({"num_genes": 1000, "num_samples": 20}, {"num_genes": 1000, "num_cell_types": 30}),
        ({"num_genes": 1000, "num_samples": 100}, {"num_genes": 50, "num_cell_types": 10}),
    ],
    indirect=True,
)
def test_check_valid_anndata(fake_bulkrnaseq_anndata, fake_signature_matrix):
    signature_matrix = SignatureMatrix(fake_signature_matrix)
    valid_anndata(fake_bulkrnaseq_anndata, signature_matrix.list_genes)


@pytest.mark.parametrize(
    "fake_bulkrnaseq_anndata, fake_signature_matrix",
    [
        ({"num_genes": 100, "num_samples": 100}, {"num_genes": 1000, "num_cell_types": 10}),
    ],
    indirect=True,
)
def test_check_non_valid_anndata(fake_bulkrnaseq_anndata, fake_signature_matrix):
    signature_matrix = SignatureMatrix(fake_signature_matrix)

    with pytest.raises(ValueError):
        valid_anndata(fake_bulkrnaseq_anndata, signature_matrix.list_genes)
