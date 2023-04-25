"""
irit io.
"""
from splinepy import splinepy_core
from splinepy.io import ioutils


def load(fname, save_replace=True):
    """
    Read spline in `.itd` form.

    Parameters
    -----------
    fname: str

      Path to the irit file to read in.
    save_replace: bool
      Replace all tabs in the irit file. Defaults to True.

    Returns
    --------
    splines: list
      Spline Type defined in NAME_TO_TYPE
    """
    if save_replace:
        ioutils.strip_tabs(fname)
    return ioutils.dict_to_spline(splinepy_core.read_irit(fname))


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
