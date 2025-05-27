import pytest

from pydeconv.signature_matrix import SignatureMatrix


@pytest.mark.parametrize(
    "fake_signature_matrix",
    [
        ({"num_genes": 1000, "num_cell_types": 10}),
        ({"num_genes": 1000, "num_cell_types": 30}),
        ({"num_genes": 50, "num_cell_types": 10}),
    ],
    indirect=True,
)
def test_check_valid_anndata(fake_signature_matrix):
    signature_matrix = SignatureMatrix(fake_signature_matrix)
    assert len(signature_matrix.list_genes) == len(fake_signature_matrix.index)
    assert len(signature_matrix.list_cell_types) == len(fake_signature_matrix.columns)

    # Update gene list
    new_genes = [f"Gene_{i}" for i in range(1, 10)]
    signature_matrix.update_gene_list(new_genes)
    assert len(signature_matrix.list_genes) == len(new_genes)
