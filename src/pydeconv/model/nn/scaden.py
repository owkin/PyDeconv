from typing import List

from pydeconv.utils import is_torch_available, requires_torch


if is_torch_available():
    import torch.nn as nn

DEFAULT_SCADEN_ARCHITECTURES = {
    "model_256": {
        "hidden_units": [256, 128, 64, 32],
        "dropout_rates": [0, 0, 0, 0],
    },
    "model_512": {
        "hidden_units": [512, 256, 128, 64],
        "dropout_rates": [0, 0.3, 0.2, 0.1],
    },
    "model_1024": {
        "hidden_units": [1024, 512, 256, 128],
        "dropout_rates": [0, 0.6, 0.3, 0.1],
    },
}


@requires_torch
class _ScadenBaseModule(nn.Module):
    """ScadenBaseModule is a simple implementation of the Scaden model from the original implementation.
    paper: https://www.nature.com/articles/s41467-022-34550-9
    """

    def __init__(self, input_dim: int, output_dim: int, hidden_units: List[int], dropout_rates: List[int]):
        """ScadenBaseModule is a simple implementation of the Scaden model from the original implementation.

        Parameters
        ----------
        input_dim : int
            Input dimension.
        output_dim : int
            Output dimension.
        hidden_units : List[int]
            Hidden units.
        dropout_rates : List[int]
            Dropout rates.
        """
        super().__init__()

        self.mlp = nn.Sequential(
            nn.Linear(input_dim, hidden_units[0]),
            nn.Dropout(dropout_rates[0]),
            nn.ReLU(),
            nn.Linear(hidden_units[0], hidden_units[1]),
            nn.Dropout(dropout_rates[1]),
            nn.ReLU(),
            nn.Linear(hidden_units[1], hidden_units[2]),
            nn.Dropout(dropout_rates[2]),
            nn.ReLU(),
            nn.Linear(hidden_units[2], hidden_units[3]),
            nn.Dropout(dropout_rates[3]),
            nn.ReLU(),
            nn.Linear(hidden_units[3], output_dim),
            nn.Softmax(dim=1),
        )

    def forward(self, x):
        """Forward pass of the model."""
        output = self.mlp(x)
        return output


@requires_torch
class ScadenModule(nn.Module):
    """ScadenModule is a simple implementation of the Scaden model from the original implementation.
    paper: https://www.nature.com/articles/s41467-022-34550-9

    This module is a combination of three ScadenBaseModule with different parameters.

    Attributes
    ----------
    model_256 : ScadenBaseModule
        ScadenBaseModule with 256 hidden units.
    model_512 : ScadenBaseModule
        ScadenBaseModule with 512 hidden units.
    model_1024 : ScadenBaseModule
        ScadenBaseModule with 1024 hidden

    Methods
    -------
    forward(x)
        Forward pass of the model.
    """

    def __init__(self, input_dim: int, output_dim: int, model_params_dict: dict):
        """ScadenModule is a simple implementation of the Scaden model from the original implementation.

        Parameters
        ----------
        input_dim : int
            Input dimension.
        output_dim : int
            Output dimension.
        model_params_dict : dict
            model parameters.
        """
        super().__init__()

        # Initialize the models with different parameters
        self.model_256 = _ScadenBaseModule(input_dim=input_dim, output_dim=output_dim, **model_params_dict["model_256"])
        self.model_512 = _ScadenBaseModule(input_dim=input_dim, output_dim=output_dim, **model_params_dict["model_512"])
        self.model_1024 = _ScadenBaseModule(
            input_dim=input_dim, output_dim=output_dim, **model_params_dict["model_1024"]
        )

    def forward(self, x):
        """Forward pass of the model.

        Parameters
        ----------
        x : torch.Tensor
            Input tensor.

        Returns
        -------
        torch.Tensor
            Output
        """
        # Get predictions from each model
        output_1 = self.model_256(x)
        output_2 = self.model_512(x)
        output_3 = self.model_1024(x)

        # Combine the outputs (e.g., averaging)
        output = (output_1 + output_2 + output_3) / 3
        return output
