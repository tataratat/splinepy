"""
vtk io.
This one samples and exports a vtk compatible geometries.
Requires sample resolutions.
Doesn't load.
"""
from splinepy import splinepy_core


def export(fname, splines, resolutions_per_spline):
    """
    Save spline as `.`.

    Parameters
    -----------
    fname: str
    splines: (n) list
    resolutions_per_spline: (n) list
      list of array-like defining

    Returns
    --------
    None
    """
    return splinepy_core.export_vtk(fname, splines, resolutions_per_spline)
