"""
irit io.
"""
from splinepy import splinepy_core


def load(fname):
    """
    Read spline in `.itd` form.

    Parameters
    -----------
    fname: str

    Returns
    --------
    splines: list
    """
    return splinepy_core.read_irit(fname)


def export(fname, splines):
    """
    Save splines as `.itd`.

    Parameters
    -----------
    fname: str
    splines: list

    Returns
    --------
    None
    """
    return splinepy_core.export_irit(fname, splines)
