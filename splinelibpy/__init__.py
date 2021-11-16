from splinelibpy import _splinelibpy as _s
from splinelibpy.bspline import BSpline
from splinelibpy.nurbs import NURBS
from splinelibpy._spline import Spline
from splinelibpy.load import load_splines

# Alias for "Legacy" classes
BSplineCurve2D = _s.BSpline1P2D
BSplineCurve3D = _s.BSpline1P3D
BSplineSurface2D = _s.BSpline2P2D
BSplineSurface3D = _s.BSpline2P3D
BSplineSolid = _s.BSpline3P3D

NurbsCurve2D = _s.NURBS1P2D
NurbsCurve3D = _s.NURBS1P3D
NurbsSurface2D = _s.NURBS2P2D
NurbsSurface3D = _s.NURBS2P3D
NurbsSolid = _s.NURBS3P3D

__all__ = [
    "Spline",
    "BSpline",
    "NURBS",
    "Reader",
    "BSplineCurve2D",
    "BSplineCurve3D",
    "BSplineSurface2D",
    "BSplineSurface3D",
    "BSplineSolid",
    "NurbsCurve2D",
    "NurbsCurve3D",
    "NurbsSurface2D",
    "NurbsSurface3D",
    "NurbsSolid",
]
