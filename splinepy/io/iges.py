"""
iges io.
"""
from splinepy import splinepy_core
from splinepy.utils import log


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
    # create empty valid spline array
    valid_splines = []
    # check if any spline contains volume elements
    for ispl, spline in enumerate(splines):
        if spline.para_dim > 2 or spline.dim > 3:
            log.warning(
                f"Skipping volumetric {type(spline)} found at splines[{ispl}]"
                "which cannot be handled by the .iges file format \n"
                "To include this spline, use a boundary representation "
                "instead, by calling e.g. spline.extract.boundaries()"
            )
        elif spline.para_dim > spline.dim:
            log.warning(
                f"Skipping {type(spline)} at splines[{ispl}]\n"
                f"Parametric dimension ({spline.para_dim}) is higher than the"
                f" physical dimension ({spline.dim}) at splines[{ispl}]\n"
                " Parametric dimension must be smaller or equal to the"
                " physical dimension."
            )
        else:
            valid_splines.append(spline)
    if not valid_splines:
        log.warning("No valid splines found in the given input")
    return splinepy_core.export_iges(fname, valid_splines)
