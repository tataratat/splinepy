"""splinepy/splinepy/load.py

Single function file containing `load_splines`.
"""

import os

from splinepy._splinepy import read_iges, read_xml, read_irit
from splinepy.bspline import BSpline
from splinepy.nurbs import NURBS
from splinepy.utils import abs_fname
from splinepy import io

def load_splines(fname, as_dict=False):
    """
    Loads spline files of extension 
      - `.iges`
      - `.xml`
      - `.itd`
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
    fname = abs_fname(fname)

    ext = os.path.splitext(fname)[1]

    if ext == ".iges":
        loaded_splines = read_iges(fname)

    elif ext == ".xml":
        loaded_splines = read_xml(fname)

    elif ext == ".itd":
        loaded_splines = read_irit(fname)

    elif ext == ".npz":
        # it should only have one spline, but
        # put it in a list to keep output format consistent
        loaded_splines = [io.npz.read_npz(fname)]

    elif ext == ".mesh":
        loaded_splines = [io.mfem.read_mfem(fname)]

    else:
        raise NotImplementedError(
            "We can only import < .iges | .xml | .itd | .npz | .mesh > "
            "spline files."
        )

    # exit early for as_dict
    if as_dict:
        return loaded_splines

    # turn dicts into real splines
    splines = []
    for s in loaded_splines:
        # are you nurbs?
        is_nurbs = s.get("weights", None)

        if is_nurbs is not None:
            splines.append(NURBS(**s))
        else:
            splines.append(BSpline(**s))
            
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
    fname = abs_fname(fname)

    ext = os.path.splitext(fname)[1]
    
    if ext == ".gf":
        return NURBS(**io.mfem.read_solution(fname, reference_spline))

    else:
        raise NotImplementedError(
            "We can only import < .gf >  solution files"
        )

