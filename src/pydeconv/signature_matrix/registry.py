from pydeconv.model.nn.registry import ORIGIN as ORIGIN_MODEL
from pydeconv.utils.hub import download_file_to_cache

from .signature_matrix import SignatureMatrix


REPO_URL = "https://raw.githubusercontent.com/owkin/PyDeconv/"
ORIGIN = "main/hub/signature_matrix/"


# CROSS TISSUE IMMUNE
def sig_matrix_cti_granularity_1():
    """Load the cross tissue immune signature matrix (1st granualrity level)."""
    path = download_file_to_cache(
        repo_url=REPO_URL, origin=ORIGIN, relative_path="sig_matrix_Cross_Tissue_Immune_granularity_1.csv"
    )
    return SignatureMatrix.load(path, index_col=0)


def sig_matrix_cti_granularity_2():
    """Load the cross tissue immune signature matrix (2nd granualrity level)."""
    path = download_file_to_cache(
        repo_url=REPO_URL, origin=ORIGIN, relative_path="sig_matrix_Cross_Tissue_Immune_granularity_2.csv"
    )
    return SignatureMatrix.load(path, index_col=0)


# LAUGHNEY LUNG CANCER
def sig_matrix_laughney_lung_cancer():
    """Load the Laugney lung cancer signature matrix."""
    path = download_file_to_cache(repo_url=REPO_URL, origin=ORIGIN, relative_path="sig_matrix_laughney_lung_cancer.csv")
    return SignatureMatrix.load(path, index_col=0)


# MIXUPVI LATENT SPACE
def signature_mixupvi_latent_space_1st_granularity():
    path = download_file_to_cache(
        repo_url=REPO_URL,
        origin=ORIGIN_MODEL,
        relative_path="mixupvi/cti_dirichlet_1st_granularity/sig_matrix_mixupvi_latent_space.csv",
    )
    return SignatureMatrix.load(path)


def signature_mixupvi_latent_space_2nd_granularity():
    path = download_file_to_cache(
        repo_url=REPO_URL,
        origin=ORIGIN_MODEL,
        relative_path="mixupvi/cti_dirichlet_2nd_granularity/sig_matrix_mixupvi_latent_space.csv",
    )
    return SignatureMatrix.load(path)
