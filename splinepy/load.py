"""splinepy/splinepy/load.py

Single function file containing `load_splines`.
"""

import os

from splinepy._splinepy import Reader
from splinepy.bspline import BSpline
from splinepy.nurbs import NURBS
from splinepy.utils import abs_fname
from splinepy import io

def load_splines(fname):
    """
    Loads spline files of extension 
      - `.iges`
      - `.xml`
      - `.itd`
      - `.npz`

    Parameters
    -----------
    fname: str

    Returns
    --------
    splines: list of BSpline or/and NURBS
    """
    fname = str(fname)
    fname = abs_fname(fname)

    sr = Reader()

    ext = os.path.splitext(fname)[1]
    
    if ext == ".iges":
        loaded_splines = sr.read_iges(fname)
    elif ext == ".xml":
        loaded_splines = sr.read_xml(fname)
    elif ext == ".itd":
        loaded_splines = sr.read_irit(fname)
    elif ext == ".npz":
        npz = np.load(fname)
        whatami = npz["whatami"][0]
        if whatami.startswith("NURBS"):
            spline = NURBS()
            spline.weights = npz["weights"]
        else:
            spline = BSpline()

        spline.control_points = npz["control_points"]
        spline.degrees = npz["degrees"]
        spline.knot_vectors = eval(npz["knot_vectors"][0])

        # For now, npz only has 1 spline. However, to keep the output format
        # consistent, return a list.
        return [spline]

    elif ext == ".mesh":
        return [NURBS(**io.mfem.read_mfem(fname))]

    else:
        raise NotImplementedError(
            "We can only import < .iges | .xml | .itd | .npz | .mesh > "
            + "spline files."
        )

    splines = []
    # Format s => [weights, degrees, knot_vectors, control_points]
    for s in loaded_splines:
        if s[0] is None:
            # BSpline.
            tmp_spline = BSpline()
            tmp_spline.degrees = s[1]
            tmp_spline.knot_vectors = s[2]
            tmp_spline.control_points = s[3]
            splines.append(tmp_spline)
 
        else:
            # Make nurbs
            tmp_spline = NURBS()
            tmp_spline.weights = s[0]
            tmp_spline.degrees = s[1]
            tmp_spline.knot_vectors = s[2]
            tmp_spline.control_points = s[3]
            splines.append(tmp_spline)

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

    sr = Reader()

    ext = os.path.splitext(fname)[1]
    
    if ext == ".gf":
        return NURBS(**io.mfem.read_solution(fname, reference_spline))

    else:
        raise NotImplementedError(
            "We can only import < .gf >  solution files"
        )

