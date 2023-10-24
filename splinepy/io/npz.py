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
    list_spline: list
    """
    loaded = np.load(fname)

    dict_spline = {}
    dict_spline.update(control_points = loaded["control_points"], degrees = loaded["degrees"])

    if "knot_vectors_0" in loaded:
        kvs = []
        for i in range(dict_spline["degrees"].size):
            kvs.append(loaded[f"knot_vectors_{i}"])
        dict_spline["knot_vectors"] = kvs

    if "weights" in loaded:
        dict_spline.update(weights=loaded["weights"])

    list_spline = [dict_spline]

    return sp.io.ioutils.dict_to_spline(list_spline)


def export(fname, spline):
    """
    Save spline as `.npz`.

    Parameters
    -----------
    fname: str
    spline: NURBS or BSpline or Bezier or Rational Bezier

    Returns
    --------
    None
    """
    property_dicts = {
        "degrees": spline.degrees,
        "control_points": spline.control_points,
        "whatami": np.array([spline.whatami])
    }

    if spline.is_rational:
        property_dicts.update(weights=spline.weights)

    if spline.has_knot_vectors:
        for i in range(spline.degrees.size):
            property_dicts.update({f"knot_vectors_{i}": spline.knot_vectors[i]})

    np.savez(
        fname,
        **property_dicts,
    )