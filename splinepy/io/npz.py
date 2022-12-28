"""
npz io.

keys in raw files are:
  - control_points : (ncpts, dim) np.ndarray, float64
  - degrees : (para_dim,) np.ndarray, int32
  - knot_vectors : (1,) np.ndarray, str
  - whatami : (1,) np.ndarray, str
"""

import numpy as np


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

    dict_spline = dict()
    # update weights
    if whatami.startswith("NURBS"):
        dict_spline.update(weights=loaded["weights"])

    # update the rest
    dict_spline.update(control_points=loaded["control_points"])
    dict_spline.update(degrees=loaded["degrees"])
    dict_spline.update(knot_vectors=eval(loaded["knot_vectors"][0]))

    return dict_spline


def export(fname, spline):
    """
    Save spline as `.npz`.

    Parameters
    -----------
    fname: str
    spline: BSpline or NURBS

    Returns
    --------
    None
    """
    property_dicts = dict(
        degrees=spline.degrees,
        knot_vectors=np.array([str(spline.knot_vectors)]),
        control_points=spline.control_points,
        whatami=np.array([spline.whatami]),
    )

    if spline.whatami.startswith("NURBS"):
        property_dicts.update(weights=spline.weights)

    np.savez(
        fname,
        **property_dicts,
    )
