import pytest
from constant import GENE_LIST_SCADEN

from pydeconv.model import MixupVI


@pytest.mark.parametrize(
    "fake_bulkrnaseq_anndata",
    [
        (
            {
                "gene_names": GENE_LIST_SCADEN,
                "num_samples": 5,
            }
        ),
        (
            {
                "gene_names": GENE_LIST_SCADEN[:-10],
                "num_samples": 5,
            }
        ),
    ],
    indirect=True,
)
def test_mixupvi(fake_bulkrnaseq_anndata):
    solver = MixupVI()

    cell_scores = solver.transform(fake_bulkrnaseq_anndata, layer="bulk_rnaseq", ratio=False)
    assert cell_scores.shape == (len(fake_bulkrnaseq_anndata.obs), len(solver.list_cell_types))
