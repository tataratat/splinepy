"""
npz io.

keys in raw files are:
  - control_points : (ncpts, dim) np.ndarray, float64
  - degrees : (para_dim,) np.ndarray, int32
  - knot_vectors : (1,) np.ndarray, str
  - whatami : (1,) np.ndarray, str
"""

import numpy as np

from splinepy import io


def load(
    fname,
):
    """
    Read a list of splines in `.npz` form.

    Parameters
    -----------
    fname: str

    Returns
    --------
    list_spline: list of NURBS / BSpline / Bezier / Rational Bezier
    """
    loaded = np.load(fname)

    # Initialize an empty list to store the splines
    list_of_spline_dicts = []

    processed_keys = set()

    # Loop through keys in loaded dictionary and find the ones starting with 'spline_'
    for key in loaded:
        if key.startswith("spline_") and key not in processed_keys:
            # Extract the index of the spline from the key
            i = int(key.split("_")[1])

            # Initialize an empty dictionary to store the properties of the spline
            dict_spline = {}

            # Add the common properties of all splines
            dict_spline["control_points"] = loaded[
                f"spline_{i}_control_points"
            ]
            dict_spline["degrees"] = loaded[f"spline_{i}_degrees"]

            # Add the weights if the spline is rational
            if f"spline_{i}_weights" in loaded:
                dict_spline["weights"] = loaded[f"spline_{i}_weights"]

            # Add the knot vectors if the spline has them
            if f"spline_{i}_knot_vectors_0" in loaded:
                kvs = []
                for j in range(dict_spline["degrees"].size):
                    kvs.append(loaded[f"spline_{i}_knot_vectors_{j}"])
                dict_spline["knot_vectors"] = kvs

            # Append dictionary of spline to the list
            list_of_spline_dicts.append(dict_spline)

            # Add keys to processed_keys set
            processed_keys.add(f"spline_{i}_control_points")
            processed_keys.add(f"spline_{i}_degrees")
            processed_keys.add(f"spline_{i}_weights")
            for j in range(dict_spline["degrees"].size):
                processed_keys.add(f"spline_{i}_knot_vectors_{j}")

    return io.ioutils.dict_to_spline(list_of_spline_dicts)


def export(fname, list_of_splines):
    """
    Save a list of splines as `.npz`.

    Parameters
    -----------
    fname: str
    splines: list of NURBS / BSpline / Bezier / Rational Bezier

    Returns
    --------
    None
    """
    # Initialize an empty dictionary to store the properties of each spline
    property_dicts = {}

    # Loop through the list of splines and add their properties to the dictionary
    for i, spline in enumerate(list_of_splines):
        # Use a prefix to distinguish different splines
        prefix = f"spline_{i}_"

        # Add the common properties of all splines
        property_dicts[f"{prefix}degrees"] = spline.degrees
        property_dicts[f"{prefix}control_points"] = spline.control_points

        # Add the weights if the spline is rational
        if spline.is_rational:
            property_dicts[f"{prefix}weights"] = spline.weights

        # Add the knot vectors if the spline has them
        if spline.has_knot_vectors:
            for j in range(spline.degrees.size):
                property_dicts[
                    f"{prefix}knot_vectors_{j}"
                ] = spline.knot_vectors[j]

    # Save the dictionary as `.npz`
    np.savez(
        fname,
        **property_dicts,
    )
