# Quickstart

The goal of this project is to have a robust python implementation of deconvolution methods for bulk RNA-seq data.

## Concepts

* [Introduction to deconvolution](<https://www.sc-best-practices.org/deconvolution/bulk_deconvolution.html>)

## 1. Load an already registered signature matrix

```python
from pydeconv.signature_matrix.registry import sig_matrix_laughney_lung_cancer, sig_matrix_laughney_lung_cancer
signature_matrix = sig_matrix_laughney_lung_cancer()
```

Checkout [here](https://github.com/owkin/PyDeconv/blob/main/pydeconv/signature_matrix/registry.py) for more description of other registered signature matrix.

## 2. Load a custom signature matrix

```python
from pydeconv import SignatureMatrix
signature_matrix = SignatureMatrix.load("path/to/signature_matrix.csv") #index: gene names, column: cell types
```

!!! note

    For the moment only `.csv` format is supported. You can add any kwargs arguments from `pd.read_csv` after the path.

## 3. Predict

```python
from pydeconv.model import Tape, Scaden

adata = AnnData("path/to/adata.h5ad") # index: sample_id, columns: gene_names
adata.layers["relative_counts"] = ...

solver = Scaden()
cell_prop = solver.transform(adata, layer="relative_counts", ratio=True)
```

!!! note

    The model will check that you have the corresponding gene names in your input data.

## 4. Predict (signature based method)

```python
from pydeconv.model import OLS, NNLS, DWLS

signature_matrix = ...
adata = AnnData("path/to/adata.h5ad")
adata.layers["relative_counts"] = ...

solver = DWLS(signature_matrix)
cell_prop = solver.transform(adata, layer="relative_counts", ratio=True)
```
