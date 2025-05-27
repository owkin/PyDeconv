from functools import wraps


try:
    import torch

    torch_available = True
except ImportError:
    torch_available = False


def is_torch_available():
    return torch_available


def get_torch():
    if not torch_available:
        raise ImportError(
            'torch is required for this functionality but not installed. (pip install "owkin-pydeconv[torch]")'
        )
    return torch


def requires_torch(cls):
    @wraps(cls)
    def wrapper(*args, **kwargs):
        get_torch()
        return cls(*args, **kwargs)

    return wrapper


try:
    import scvi

    scvi_available = True
except ImportError:
    scvi_available = False


def is_scvi_available():
    return scvi_available


def get_scvi():
    if not scvi_available:
        raise ImportError(
            'scvi is required for this functionality but not installed. (pip install "owkin-pydeconv[torch]")'
        )
    return scvi


def requires_scvi(cls):
    @wraps(cls)
    def wrapper(*args, **kwargs):
        get_scvi()
        return cls(*args, **kwargs)

    return wrapper
