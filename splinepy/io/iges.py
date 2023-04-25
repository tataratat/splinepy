"""
iges io.
"""
from splinepy import splinepy_core
from splinepy.io import ioutils


def load(fname):
    """
    Read spline in `.iges` form.

    Parameters
    -----------
    fname: str

    Returns
    --------
    splines: list
      Spline Type defined in NAME_TO_TYPE
    """
    return ioutils.dict_to_spline(splinepy_core.read_iges(fname))


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
