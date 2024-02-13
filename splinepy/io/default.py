import pathlib

from splinepy.io import gismo as _gismo
from splinepy.io import iges as _iges
from splinepy.io import irit as _irit
from splinepy.io import json as _json
from splinepy.io import mfem as _mfem
from splinepy.io import npz as _npz
from splinepy.io import svg as _svg

extension_to_io = {
    ".xml": _gismo,
    ".iges": _iges,
    ".igs": _iges,
    ".itd": _irit,
    ".irit": _irit,
    ".json": _json,
    ".mesh": _mfem,
    ".npz": _npz,
    ".svg": _svg,
}


def export(fname, splinepy_obj):
    """
    Exports a splinepy object to a file.

    Tries to find out the correct exporter from the file extension. If the
    extension is not supported, a ValueError is raised.

    Parameters
    ----------
    fname : str
      Export Filename
    splinepy_obj: Splinepy object
      Splinepy object to be exported.

    Returns
    -------
    output_dict : dict
      Dictornay data written into file

    """
    # Determine file extension
    extension = pathlib.Path(fname).suffix
    if extension not in extension_to_io:
        raise ValueError("Extension not supported")

    # Export
    extension_to_io[extension].export(fname, splinepy_obj)


def load(fname):
    """
    Loads a splinepy object from a file.

    Tries to find out the correct importer from the file extension. If the
    extension is not supported, a ValueError is raised.

    Parameters
    ----------
    fname : str
      Filename to load

    Returns
    -------
    splinepy_obj : Splinepy object
      Splinepy object loaded from file

    """
    # Determine file extension
    extension = pathlib.Path(fname).suffix
    if extension not in extension_to_io:
        raise ValueError("Extension not supported")

    # Load
    return extension_to_io[extension].load(fname)
