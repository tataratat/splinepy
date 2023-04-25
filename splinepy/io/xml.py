"""
xml io.
Keyword follows RWTH CATS's spline format.
Possiblely will be turned into pure python io.
"""
from splinepy import splinepy_core
from splinepy.io import ioutils


def load(fname):
    """
    Read spline in `.xml` form.
    RWTH CATS spline format.

    Parameters
    -----------
    fname: str

    Returns
    --------
    splines: list
      Spline Type defined in NAME_TO_TYPE
    """
    return ioutils.dict_to_spline(splinepy_core.read_xml(fname))


def export(fname, splines):
    """
    Save spline as `.xml`.

    Parameters
    -----------
    fname: str
    splines: list

    Returns
    --------
    None
    """
    return splinepy_core.export_xml(fname, splines)
