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
    # check if any spline contains volume elements
    for ispl, spline in enumerate(splines):
        if spline.para_dim > 2 or spline.dim > 3:
            raise ValueError(
                f"Volumetric {type(spline)} found at splines[{ispl}], which cannot be handled by the .iges file format \n"
                "Use a boundary representation instead by calling e.g. multipatch.boundary_multipatch()"
            )
    return splinepy_core.export_iges(fname, splines)
