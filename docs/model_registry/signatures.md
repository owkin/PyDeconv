# Signature registry

The signature matrix is a matrix that contains the expression of genes for each cell type. It is used to predict the cell type composition of a sample. PyDeconv provides a registry of signature matrices that can be used to predict the cell type composition of a sample.

The signature matrices are stored in CSV files and can be loaded using the `pydeconv.signature_matrix.registry` module.

!!! success "Commercialy available signature matrices"
    The following signature matrices are commercially available. You can use them without any restrictions.

    ### Laugney Lung Cancer

    * indication: `lung`
    * dataset used:
    * n patients:
    * n cells:
    * extraction comments:

    ```python
    from pydeconv.signature_matrix.registry import sig_matrix_laughney_lung_cancer
    ```

    ??? info "Output Celltypes"

        `TNK`, `B`, `MonoMacro`, `Mast`, `Endothelial`, `Malignant`, `Stroma`, `DC`, `Epithelial`

    ### Cross tissue immune main CP

    * indication:
    * dataset used:
    * n patients:
    * n cells:
    * extraction comments:

    ```python
    from pydeconv.signature_matrix.registry import sig_matrix_cti_granularity_1
    ```

    ??? info "Output Celltypes"

        `TNK`, `B`, `MonoMacro`, `DC`, `Mast`

    ### Cross tissue immune main 2nd granularity

    * indication:
    * dataset used:
    * n patients:
    * n cells:
    * extraction comments:

    ```python
    from pydeconv.signature_matrix.registry import sig_matrix_cti_granularity_2
    ```

    ??? info "Output Celltypes"

        `CD4T`, `CD8T`, `B`, `Plasma`, `DC`, `NK`, `Mono`, `Mast`, `Tregs`

    ### Cross tissue immune main 3rd granularity

    * indication:
    * dataset used:
    * n patients:
    * n cells:
    * extraction comments:

---

!!! warning "Research only signature matrices"
    The following signature matrices are for research purposes only.

---

!!! danger "Private only signature matrices"
    The following signature matrices are private and can only be used by
