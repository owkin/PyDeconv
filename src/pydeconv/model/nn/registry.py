import json
from typing import Literal, Tuple

from pydeconv.utils import is_torch_available
from pydeconv.utils.hub import download_file_to_cache

from .mixupvi import MixupVIModule
from .scaden import DEFAULT_SCADEN_ARCHITECTURES, ScadenModule
from .tape import TapeModule


if is_torch_available():
    import torch
    import torch.nn as nn

REPO_URL = "https://raw.githubusercontent.com/owkin/PyDeconv/"
ORIGIN = "main/hub/model/"


def scaden_module_cti(
    granularity: Literal["1st_granularity", "2nd_granularity"], device="cpu"
) -> Tuple[nn.Module, dict]:
    """Factory function to create a ScadenModule with predefined parameters.

    Parameters
    ----------
    version :str
        The version of the ScadenModule to create.

    Returns
    ----------
    nn.Module
        An instance of ScadenModule with the specified version.
    """
    path_parameters = download_file_to_cache(
        repo_url=REPO_URL,
        origin=ORIGIN,
        relative_path=f"scaden/cti_dirichlet_{granularity}/scaden_params.json",
    )
    with open(path_parameters) as f:
        params = json.load(f)

    model = ScadenModule(
        input_dim=params["input"]["input_dim"],
        output_dim=params["output"]["output_dim"],
        model_params_dict=DEFAULT_SCADEN_ARCHITECTURES,
    )
    for shape in ["256", "512", "1024"]:
        path_model = download_file_to_cache(
            repo_url=REPO_URL,
            origin=ORIGIN,
            relative_path=f"scaden/cti_dirichlet_{granularity}/scaden_weights_{shape}.pth",
        )
        state_dict = torch.load(
            path_model,
            map_location=device,
            weights_only=True,
        )
        for key in list(state_dict.keys()):
            state_dict[key.replace("model.", "mlp.")] = state_dict.pop(key)

        getattr(model, f"model_{shape}").load_state_dict(state_dict)
    return model, params


def tape_module_cti(
    granularity: Literal["1st_granularity", "2nd_granularity"], adaptative, device="cpu"
) -> Tuple[nn.Module, dict]:
    """Factory function to create a ScadenModule with predefined parameters.

    Parameters
    ----------
    version :str
        The version of the ScadenModule to create.

    Returns
    ----------
    nn.Module: An instance of ScadenModule with the specified version.
    """
    path_parameters = download_file_to_cache(
        repo_url=REPO_URL,
        origin=ORIGIN,
        relative_path=f"tape/cti_dirichlet_{granularity}/tape_params.json",
    )
    with open(path_parameters) as f:
        params = json.load(f)

    model = TapeModule(
        input_dim=params["input"]["input_dim"], output_dim=params["output"]["output_dim"], adaptative=adaptative
    )

    path_model = download_file_to_cache(
        repo_url=REPO_URL,
        origin=ORIGIN,
        relative_path=f"tape/cti_dirichlet_{granularity}/tape_weights.pth",
    )
    state_dict = torch.load(
        path_model,
        map_location=device,
        weights_only=True,
    )
    model.load_state_dict(state_dict)
    return model, params


def mixupvi_module_cti(
    granularity: Literal["1st_granularity", "2nd_granularity"], device="cpu"
) -> Tuple[nn.Module, dict]:
    """Factory function to create a MixupVIModule with predefined parameters.
    Parameters
    ----------
    granularity : Literal["1st_granularity", "2nd_granularity"]
        The granularity level for the MixupVIModule.
    device : str, optional
        The device to load the model on, by default "cpu".
    Returns
    -------
    Tuple[nn.Module, dict]
        A tuple containing the MixupVIModule instance and its parameters.
    """
    path_parameters = download_file_to_cache(
        repo_url=REPO_URL,
        origin=ORIGIN,
        relative_path=f"mixupvi/cti_dirichlet_{granularity}/mixupvi_params.json",
    )
    with open(path_parameters) as f:
        params = json.load(f)

    model = MixupVIModule(
        input_dim=params["input"]["input_dim"],
        latent_dim=params["output"]["latent_dim"],
    )

    path_model = download_file_to_cache(
        repo_url=REPO_URL,
        origin=ORIGIN,
        relative_path=f"mixupvi/cti_dirichlet_{granularity}/mixupvi_weights.pt",
    )

    model_dict = torch.load(
        path_model,
        map_location=device,
        weights_only=False,
    )

    model.model.load_state_dict(model_dict["model_state_dict"])
    return model, params
