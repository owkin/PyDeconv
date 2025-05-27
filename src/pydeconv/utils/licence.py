import functools
import warnings


def non_commercial_only(func):
    """
    Decorator that raises a warning about the license when the function is called.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        warnings.warn(
            f"Warning: The {func.__name__} is only for non-commercial use. "
            "Make sure you are not working on a commercial project.",
            category=UserWarning,
            stacklevel=2,
        )
        return func(*args, **kwargs)

    return wrapper


def internal_only(func):
    """
    Decorator that raises a warning about the license when the function is called.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        warnings.warn(
            f"Warning: The {func.__name__} is only for private use. Make sure you are not sharing the results.",
            category=UserWarning,
            stacklevel=2,
        )
        return func(*args, **kwargs)

    return wrapper
