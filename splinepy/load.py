"""splinepy/splinepy/load.py

Single function file containing `load_splines`.
"""

import os as _os

from splinepy import io as _io
from splinepy.io.ioutils import abs_fname as _abs_fname
from splinepy.io.ioutils import dict_to_spline as _dict_to_spline
from splinepy.nurbs import NURBS as _NURBS
from splinepy.splinepy_core import read_iges as _read_iges


def load_splines(fname, as_dict=False):
    """
    Loads spline files of extension
      - `.iges`
      - `.xml`
      - `.itd`
      - `.json`
      - `.npz`
      - `.mesh` #only 2D-single-patch.

    Parameters
    -----------
    fname: str
    as_dict: bool
      Default is False. If True, returns spline properties as dict.

    Returns
    --------
    splines: list of BSpline/NURBS or dict
    """
    fname = str(fname)
    fname = _abs_fname(fname)

    ext = _os.path.splitext(fname)[1]

    if ext in {".iges", ".igs"}:
        loaded_splines = _read_iges(fname)

    elif ext == ".itd":
        loaded_splines = _io.load.irit(fname)

    elif ext == ".npz":
        # it should only have one spline, but
        # put it in a list to keep output format consistent
        loaded_splines = [_io.npz.load(fname)]

    elif ext == ".mesh":
        loaded_splines = [_io.mfem.load(fname)]

    elif ext == ".json":
        json_load_raw = _io.json.load(fname)
        # json load returns splines organized by their types in dict
        # pretty neat. but for now, let's keep the format consistent
        loaded_splines = []
        for splines in json_load_raw.values():
            loaded_splines.extend(splines)

    else:
        raise NotImplementedError(
            "We can only import "
            "< .iges | .xml | .itd | .npz | .mesh | .json > "
            "spline files."
        )

    # exit early for as_dict
    if as_dict:
        return loaded_splines

    # turn dicts into real splines
    splines = _dict_to_spline(loaded_splines)

    return splines


def load_solution(fname, reference_spline):
    """
    Reads a solution file and returns a solution spline.
    If output is spline, use load_splines.

    Parameters
    -----------
    fname: str
    reference_spline: Spline

    Returns
    --------
    solution_spline: Spline
    """
    fname = str(fname)
    fname = _abs_fname(fname)

    ext = _os.path.splitext(fname)[1]

    if ext == ".gf":
        return _NURBS(**_io.mfem.read_solution(fname, reference_spline))

    else:
        raise NotImplementedError("We can only import < .gf >  solution files")
