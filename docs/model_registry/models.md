# Model registry

## Models

---

### Scaden

Paper: [Deep learning–based cell composition analysis from tissue expression profiles](https://www.science.org/doi/10.1126/sciadv.aba2619)

```python
from pydeconv.model import Scaden
model = Scaden()
model.transform(adata, layer="tpm", ratio=True)
```

#### Registered models

* `cti_dirichlet_2nd_granularity`

    ... details dataset and model ...

    ??? info "Output Celltypes"
        `B`, `CD4T`, `CD8T`, `DC`, `Mast`, `Mono`, `NK`, `Plasma`, `Tregs`

### TAPE

Paper: [Deep autoencoder for interpretable tissue-adaptive deconvolution and cell-type-specific gene analysis](https://www.nature.com/articles/s41467-022-34550-9)

```python
from pydeconv.model import Tape
model = Tape()
model.transform(adata, layer="tpm", ratio=True)
```

#### Registered models

* `cti_dirichlet_2nd_granularity`

    ... details dataset and model ...

    ??? info "Output Celltypes"
        `B`, `CD4T`, `CD8T`, `DC`, `Mast`, `Mono`, `NK`, `Plasma`, `Tregs`

## Signature based models

---

### OLS

Paper:

```python
from pydeconv.model import OLS
model = OLS(signature_matrix)
model.transform(adata, layer="tpm", ratio=True)
```

### RLR

Paper:

```python
from pydeconv.model import RLR
model = RLR(signature_matrix)
model.transform(adata, layer="tpm", ratio=True)
```

### NNLS

Paper: [Deconvolution of blood microarray data identifies cellular activation patterns in systemic lupus erythematosus](https://pubmed.ncbi.nlm.nih.gov/19568420/)

```python
from pydeconv.model import NNLS
model = NNLS(signature_matrix)
model.transform(adata, layer="tpm", ratio=True)
```

### DWLS

Paper: [Accurate estimation of cell-type composition from gene expression data](https://www.nature.com/articles/s41467-019-10802-z)

```python
from pydeconv.model import DWLS
model = DWLS(signature_matrix)
model.transform(adata, layer="tpm", ratio=True)
```

### WNNLS (MuSiC)

Paper: [Bulk tissue cell type deconvolution with multi-subject single-cell expression reference](https://pubmed.ncbi.nlm.nih.gov/30670690/)

```python
from pydeconv.model import WNNLS
model = WNNLS(signature_matrix)
model.transform(adata, layer="tpm", ratio=True)
```
