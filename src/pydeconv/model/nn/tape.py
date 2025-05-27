import copy
from typing import Tuple, Union

from pydeconv.utils import is_torch_available, requires_torch, seed_all


if is_torch_available():
    import torch
    import torch.nn as nn
    import torch.nn.functional as F


@requires_torch
class TapeModule(nn.Module):
    """TapeModule is a simple implementation of the TAPE model from the original implementation.

    paper: https://www.nature.com/articles/s41467-022-34550-9
    """

    def __init__(self, input_dim: int, output_dim: int, adaptative: bool = False):
        """TapeModule is a simple implementation of the TAPE model from the original implementation.

        Parameters
        ----------
        input_dim : int
            Input dimension.
        output_dim : int
            Output dimension.
        adaptative : bool, optional
            Adaptative step for predicting the model, by default True.
        """
        super().__init__()
        self.adaptative_step: bool = adaptative

        self.encoder = nn.Sequential(
            nn.Dropout(),
            nn.Linear(input_dim, 512),
            nn.CELU(),
            nn.Dropout(),
            nn.Linear(512, 256),
            nn.CELU(),
            nn.Dropout(),
            nn.Linear(256, 128),
            nn.CELU(),
            nn.Dropout(),
            nn.Linear(128, 64),
            nn.CELU(),
            nn.Linear(64, output_dim),
        )

        # not used in the predit method
        self.decoder = nn.Sequential(
            nn.Linear(output_dim, 64, bias=False),
            nn.Linear(64, 128, bias=False),
            nn.Linear(128, 256, bias=False),
            nn.Linear(256, 512, bias=False),
            nn.Linear(512, input_dim, bias=False),
            nn.ReLU(),
        )
        self.activation = nn.ReLU()

    @property
    def signature_matrix(self) -> torch.Tensor:
        """Pseudo property to get the signature matrix from the model.

        Returns
        -------
        torch.Tensor
            Signature matrix from the model.
        """
        sm = torch.mm(self.decoder[0].weight.T, self.decoder[1].weight.T)
        sm = torch.mm(sm, self.decoder[2].weight.T)
        sm = torch.mm(sm, self.decoder[3].weight.T)
        sm = torch.mm(sm, self.decoder[4].weight.T)
        sm = F.relu(sm)
        return sm

    def forward(self, x: torch.Tensor) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """Forward pass of the model.

        Parameters
        ----------
        x : torch.Tensor
            Input tensor.

        Returns
        -------
        Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]
            Output tensor or tuple of output tensor and reconstructed
            tensor.
        """
        z = self.encoder(x)
        if self.training:
            x_reconstruct = torch.mm(z, self.signature_matrix)
            output = (z, x_reconstruct)
        else:
            if self.adaptative_step:
                model_adapted = adaptative_step(self, x, device=x.device)
                model_adapted.eval()
                output = model_adapted(x)
            else:
                output = self.activation(z)
        return output


def adaptative_step(
    model: TapeModule,
    data: torch.Tensor,
    n_steps: int = 10,
    n_iter: int = 5,
    device: Union[str, torch.device] = "cpu",
    lr: float = 1e-4,
) -> nn.Module:
    """Adaptive stage for training the model from the original implementation.
    https://github.com/poseidonchan/TAPE/blob/8ffb2f4600e1cbc689c6b1b1f428e1ddac773c6e/TAPE/train.py#L42

    Parameters
    ----------
    model : TapeModule
        The model to adapt.
    data : torch.Tensor
        Input data.
    n_steps : int, optional
        Number of steps to perform, by default 10
    n_iter : int, optional
        Number of iterations to perform, by default 5
    device : Union[str, torch.device], optional
        Device to use for the computation, by default "cpu"
    lr : float, optional
        Learning rate, by default 1e-4
    """
    data = data.to(device)
    # Because the model is trained in place, we need to copy the model to not modify the original model weights
    # if we perform multiple predicts
    model_copy = copy.deepcopy(model)
    model_copy.eval()
    model_copy.adaptative_step = False

    opt_encoder = torch.optim.Adam(model_copy.encoder.parameters(), lr=lr)
    opt_decoder = torch.optim.Adam(model_copy.decoder.parameters(), lr=lr)

    origin_pred = model_copy(data)
    origin_sigmatrix = model_copy.signature_matrix.detach()
    origin_pred = origin_pred.detach()

    for _ in range(n_iter):
        model_copy.train()
        for _ in range(n_steps):
            seed_all(seed=0)
            opt_decoder.zero_grad()
            _, x_recon = model_copy(data)
            batch_loss = F.l1_loss(x_recon, data) + F.l1_loss(model_copy.signature_matrix, origin_sigmatrix)
            batch_loss.backward()
            opt_decoder.step()

        for _ in range(n_steps):
            seed_all(seed=0)
            opt_encoder.zero_grad()
            pred, x_recon = model_copy(data)
            batch_loss = F.l1_loss(origin_pred, pred) + F.l1_loss(x_recon, data)
            batch_loss.backward()
            opt_encoder.step()

    return model_copy
