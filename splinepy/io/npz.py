"""
npz io.

keys in raw files are:
  - control_points : (ncpts, dim) np.ndarray, float64
  - degrees : (para_dim,) np.ndarray, int32
  - knot_vectors : (1,) np.ndarray, str
  - whatami : (1,) np.ndarray, str
"""

import numpy as np
import splinepy as sp


def load(
    fname,
):
    """
    Read spline in `.npz` form.

    Parameters
    -----------
    fname: str

    Returns
    --------
    dict_spline: dict
    """
    loaded = np.load(fname)
    whatami = loaded["whatami"][0]

    # create NURBS
    if whatami.startswith("NURBS"):
        loaded_spline = sp.NURBS(
        degrees=loaded["degrees"],
        knot_vectors=loaded["knot_vectors"],
        control_points=loaded["control_points"],
        weights=loaded["weights"],
        )

    # create BSpline
    elif whatami.startswith("BSpline"):
        loaded_spline = sp.BSpline(
        degrees=loaded["degrees"],
        knot_vectors=loaded["knot_vectors"],
        control_points=loaded["control_points"],
        )

    # create Bezier
    elif whatami.startswith("Bezier"):
        loaded_spline = sp.Bezier(
        degrees=loaded["degrees"],
        control_points=loaded["control_points"],
        )

    # create Rational Bezier
    elif whatami.startswith("RationalBezier"):
        loaded_spline = sp.RationalBezier(
        degrees=loaded["degrees"],
        control_points=loaded["control_points"],
        weights=loaded["weights"],
        )

    else:
        raise TypeError("No proper spline type detected.")

    return loaded_spline


def export(fname, spline):
    """
    Save spline as `.npz`.

    Parameters
    -----------
    fname: str
    spline: Bezier or BSpline or NURBS

    Returns
    --------
    None
    """
    property_dicts = {
        "degrees": spline.degrees,
        "knot_vectors": np.array(spline.knot_vectors),
        "control_points": spline.control_points,
        "whatami": np.array([spline.whatami]),
    }

    if spline.whatami.startswith("NURBS"):
        property_dicts.update(weights=spline.weights)

    np.savez(
        fname,
        **property_dicts,
    )
