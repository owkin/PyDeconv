import random
from pathlib import Path
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import requests
import scanpy as sc
from joblib import Memory
from tqdm import tqdm


memory = Memory(location=".cache", verbose=0)


def load_pseudobulk_dataset(n_sample: int, n_cells: int, granularity: str, folder_data: Path):
    anndata = load_cti(folder_data=folder_data)
    anndata, (_, index_test) = add_cell_types_grouped(anndata, group=granularity, folder_data=folder_data)
    pseudo_bulk_test, labels = create_dirichlet_pseudobulk_dataset(
        anndata[index_test], n_cells=n_cells, n_sample=n_sample
    )
    return pseudo_bulk_test, labels


@memory.cache
def load_cti(folder_data: Path, **kwargs):
    """Load and preprocess the CTI scRNAseq dataset.

    Parameters
    ----------
    folder_data: Path
        The folder where the CTI dataset is stored or will be downloaded.
    **kwargs: dict
        Additional keyword arguments for preprocessing.

    Return
    ------
    ad.AnnData
        The preprocessed CTI dataset as an AnnData object.
    """
    if not (folder_data / "cti_adata_global.h5ad").exists():
        download_cti(folder_data=folder_data)
    adata = sc.read(folder_data / "cti_adata_global.h5ad")
    adata.X = adata.raw.X
    adata = preprocess_scrna(adata, batch_key="donor_id")
    return adata


def download_cti(folder_data: Path):
    url = "https://datasets.cellxgene.cziscience.com/51007192-dc1e-4718-b48c-dbf35e8216f7.h5ad"
    output_path = folder_data / "cti_adata_global.h5ad"

    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        total_size = int(r.headers.get("content-length", 0))

        with (
            open(output_path, "wb") as f,
            tqdm(total=total_size, unit="B", unit_scale=True, desc="Downloading CTI dataset") as pbar,
        ):
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:  # filter out keep-alive chunks
                    f.write(chunk)
                    pbar.update(len(chunk))


def preprocess_scrna(adata: ad.AnnData, keep_genes: int = 2000, batch_key: Optional[str] = None):
    """Preprocess single-cell RNA data for deconvolution benchmarking.

    Parameters
    ----------
    adata: ad.AnnData
        in adata.X, the normalized log1p counts are saved
        in adata.layers["raw_counts"], raw counts are saved
        in adata.layers["relative_counts"], the relative counts are saved
        => The highly variable genes can be found in adata.var["highly_variable"]
    keep_genes: int | None
        The number of most variable genes to keep. If None, all genes are kept.
    batch_key: str | None
        The key in `adata.obs` to use for batch correction. If None, no batch correction is applied.
    Returns
    -------
    ad.AnnData
        The preprocessed AnnData object.
    """
    sc.pp.filter_genes(adata, min_counts=3)
    adata.layers["raw_counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    adata.layers["relative_counts"] = adata.X.copy()
    sc.pp.log1p(adata)
    adata.raw = adata

    return adata


def add_cell_types_grouped(adata: ad.AnnData, group: str, folder_data: Path) -> ad.AnnData:
    """Add the cell types grouped columns in Anndata according to the grouping choice.
    It uses and returns the train_test_index csv file created for the signature matrix.
    """
    if group == "1st_level_granularity":
        idx_split_group = pd.read_csv(folder_data / "train_test_index_matrix_common.csv", index_col=0)
    elif group == "2nd_level_granularity":
        idx_split_group = pd.read_csv(folder_data / "train_test_index_matrix_granular_updated.csv", index_col=0)
    adata.obs["cell_types_grouped"] = idx_split_group["group"]

    train_index = list(idx_split_group[idx_split_group.train_index].index)
    test_index = list(idx_split_group[idx_split_group.test_index].index)
    return adata, (train_index, test_index)


def create_dirichlet_pseudobulk_dataset(
    adata: ad.AnnData,
    prior_alphas: np.array = None,
    n_sample: int = 400,
    cell_type_group: str = "cell_types_grouped",
    aggregation_method: str = "mean",
    n_cells: int = 100,
    is_n_cells_random: bool = False,
):
    """Create pseudobulk dataset from single-cell RNA data, sampled from a dirichlet
    distribution. If a prior belief on the cell fractions (e.g. prior knowledge from
    specific tissue), then it can be incorporated. Otherwise, it will just be a non-
    informative prior. Then, compute dirichlet posteriors to sample cells - dirichlet is
    conjugate to the multinomial distribution, thus giving an easy posterior
    calculation.
    """
    seed = random.randint(0, 1000)
    random_state = np.random.RandomState(seed=seed)
    cell_types = adata.obs[cell_type_group].value_counts()
    if prior_alphas is None:
        prior_alphas = np.ones(len(cell_types))  # non-informative prior
    likelihood_alphas = cell_types / adata.n_obs  # multinomial likelihood
    alpha_posterior = prior_alphas + likelihood_alphas
    posterior_dirichlet = random_state.dirichlet(alpha_posterior, n_sample)
    if is_n_cells_random:
        n_cells = np.random.randint(50, 1001, size=posterior_dirichlet.shape[0])
        posterior_dirichlet = np.round(np.multiply(posterior_dirichlet, n_cells))
    else:
        posterior_dirichlet = np.round(posterior_dirichlet * n_cells)
    posterior_dirichlet = posterior_dirichlet.astype(np.int64)  # number of cells to sample
    groundtruth_fractions = posterior_dirichlet / posterior_dirichlet.sum(axis=1, keepdims=True)

    random.seed(seed)
    averaged_data, _ = {"relative_counts": [], "raw_counts": [], "counts_sum": []}, []
    all_adata_samples = []
    for i in range(n_sample):
        sample_data = []
        for j, cell_type in enumerate(likelihood_alphas.index):
            # If sample larger than cell population, sample with replacement
            if posterior_dirichlet[i][j] > cell_types[cell_type]:
                cell_sample = random.choices(
                    list(adata.obs.loc[adata.obs.cell_types_grouped == cell_type].index),
                    k=posterior_dirichlet[i][j],
                )
            else:
                cell_sample = random.sample(
                    list(adata.obs.loc[adata.obs.cell_types_grouped == cell_type].index),
                    posterior_dirichlet[i][j],
                )
            sample_data.extend(cell_sample)
        adata_sample = adata[sample_data]
        if aggregation_method == "mean":
            averaged_data["relative_counts"].append(adata_sample.layers["relative_counts"].mean(axis=0).tolist()[0])
            X = np.array(adata_sample.layers["raw_counts"].mean(axis=0).tolist()[0])
            X_sum = np.array(adata_sample.layers["raw_counts"].sum(axis=0).tolist()[0])

            averaged_data["raw_counts"].append(X)
            averaged_data["counts_sum"].append(X_sum)
        all_adata_samples.append(adata_sample)

    # pseudobulk dataset
    adata_pseudobulk = create_anndata_pseudobulk(adata.obs, adata.var_names, np.array(averaged_data["raw_counts"]))
    adata_pseudobulk_rc = create_anndata_pseudobulk(
        adata.obs, adata.var_names, np.array(averaged_data["relative_counts"])
    )
    adata_pseudobulk_cs = create_anndata_pseudobulk(adata.obs, adata.var_names, np.array(averaged_data["counts_sum"]))
    adata_pseudobulk.layers["relative_counts"] = adata_pseudobulk_rc.to_df()
    adata_pseudobulk.layers["counts_sum"] = adata_pseudobulk_cs.to_df()

    # ground truth fractions
    groundtruth_fractions = pd.DataFrame(
        groundtruth_fractions, index=adata_pseudobulk.obs_names, columns=list(cell_types.index)
    )
    groundtruth_fractions = groundtruth_fractions.fillna(0)  # The Nan are cells not sampled

    return adata_pseudobulk, groundtruth_fractions


def create_anndata_pseudobulk(adata_obs: pd.DataFrame, adata_var_names: list, x: np.array) -> ad.AnnData:
    """Creates an anndata object from a pseudobulk sample.

    Parameters
    ----------
    adata_obs: pd.DataFrame
        Obs dataframe from anndata object storing training set
    adata_var_names: list
        Gene names from the anndata object
    x: np.array
        pseudobulk sample

    Return
    ------
    ad.AnnData
        Anndata object storing the pseudobulk array
    """
    df_obs = pd.DataFrame.from_dict([{col: adata_obs[col].value_counts().index[0] for col in adata_obs.columns}])
    if len(x.shape) > 1 and x.shape[0] > 1:
        # several pseudobulks, so duplicate df_obs row
        df_obs = df_obs.loc[df_obs.index.repeat(x.shape[0])].reset_index(drop=True)
        df_obs.index = [f"sample_{idx}" for idx in df_obs.index]
    adata_pseudobulk = ad.AnnData(X=x, obs=df_obs)
    adata_pseudobulk.var_names = adata_var_names
    adata_pseudobulk.layers["raw_counts"] = np.copy(x)
    adata_pseudobulk.raw = adata_pseudobulk

    return adata_pseudobulk
