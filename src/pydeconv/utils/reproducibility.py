import random

import numpy as np

from .optional_imports import is_torch_available


if is_torch_available():
    import torch


def seed_all(seed=0):
    random.seed(seed)
    np.random.seed(seed)
    if is_torch_available():
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
            torch.backends.cudnn.deterministic = True
