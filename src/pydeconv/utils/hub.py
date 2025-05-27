import logging
import os
import posixpath
import shutil
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from pydeconv import __version__


def get_cache_path(relative_path: str, version: str = __version__, libname: str = "pydeconv"):
    cache_base = os.getenv("XDG_CACHE_HOME", os.path.expanduser("~/.cache"))
    return Path(cache_base) / libname / version / relative_path


def download_file_to_cache(repo_url: str, origin: str, relative_path: str, version=__version__, libname="pydeconv"):
    """Download a file from a GitHub repository and cache it locally.
    Parameters
    ----------
    repo_url : str
        The URL of the GitHub repository (e.g., "
    origin : str
        The origin of the file, e.g., "pydeconv.hub.model".
    relative_path : str
        The relative path to the file in the repository.
    version : str, optional
        The version of the library, by default __version__ from pydeconv.
    libname : str, optional
        The name of the library, by default "pydeconv".
    Returns
    -------
    str
        The local path to the cached file.
    Raises
    ------
    RuntimeError
        If the download fails (non-200 status code).
    """
    cache_path = get_cache_path(relative_path, version, libname)
    if cache_path.exists():
        logging.log(msg=f"[✓] Using cached file: {cache_path}", level=logging.DEBUG)
        return str(cache_path)

    url = posixpath.join(repo_url, origin, relative_path)
    logging.log(msg=f"[↓] Downloading from: {url}", level=logging.DEBUG)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urlopen(req) as response, open(cache_path, "wb") as out_file:
            shutil.copyfileobj(response, out_file)
    except HTTPError as e:
        raise RuntimeError(f"Failed to download: {url} (HTTP {e.code})") from e
    except URLError as e:
        raise RuntimeError(f"Failed to download: {url} (URL Error: {e.reason})") from e

    return str(cache_path)
