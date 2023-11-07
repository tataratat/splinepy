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
    try:
        return splinepy_core.export_iges(fname, splines)

    except RuntimeError:
        # create empty valid spline array
        valid_splines = []
        # check if any spline contains volume elements
        for ispl, spline in enumerate(splines):
            if spline.para_dim > 2:
                log.warning(
                    f"Skipping volumetric {type(spline)} at splines[{ispl}]"
                    "which cannot be handled by the .iges file format \n"
                    "To include this spline, use a boundary representation "
                    "instead, by calling e.g. spline.extract.boundaries()"
                )
            elif spline.dim > 3:
                log.warning(
                    f"Skipping high dim {type(spline)} at splines[{ispl}]"
                    "which cannot be handled by the .iges file format."
                    "iges supports up to dim=3."
                )
            else:
                valid_splines.append(spline)

        if len(valid_splines) == 0:
            raise ValueError("No valid splines found in the given input")
