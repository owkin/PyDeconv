import warnings

from anndata import AnnData


def valid_anndata(adata: AnnData, list_genes, tol=0.3):
    """Check if genes in signature matrix are in anndata object.

    Parameters
    ----------
    adata : AnnData
        Anndata object containing the data. Should have the same genes as the signature matrix.
        Rows are samples, columns are genes.
    list_genes : list
        List of genes to check.
    tol : float, optional
        Tolerance for missing genes, by default 0.3
    """
    diff_genes = set(list_genes) - set(adata.var_names)

    ratio_missing = len(diff_genes) / len(list_genes)

    if len(diff_genes) == 0:
        return adata[:, list_genes], None
    elif ratio_missing <= tol:
        warnings.warn(
            f"Genes in signature matrix/model not found in anndata object, missing: {diff_genes}", stacklevel=2
        )
        inter_genes = [gene for gene in list_genes if gene in adata.var_names]  # do not use set to keep ordered genes
        return adata[:, inter_genes], inter_genes
    else:
        raise ValueError(
            rf"More than {tol*10:.0f}\% of genes in signature matrix/model not found in anndata object, "
            " check that you are using the right indexes"
        )
