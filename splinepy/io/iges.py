"""
iges io.
"""
from splinepy import splinepy_core


def load(fname):
    """
    Read spline in `.iges` form.

    Parameters
    -----------
    fname: str

    Returns
    --------
    dict_splines: list
    """
    return splinepy_core.read_iges(fname)


def export(fname, splines):
    """
    Save splines as `.iges`.

    Parameters
    -----------
    fname: str
    spline: list
      list of Splines

    Returns
    --------
    None
    """
    return splinepy_core.export_iges(fname, splines)
