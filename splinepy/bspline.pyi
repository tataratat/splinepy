import numpy.typing as _npt
from _typeshed import Incomplete

from splinepy import spline as _spline

has_scipy: bool

class BSplineBase(_spline.Spline):
    def __init__(self, *args, **kwargs) -> None: ...
    def insert_knots(self, parametric_dimension, knots): ...
    def knot_insertion_matrix(
        self,
        parametric_dimension: Incomplete | None = None,
        knots: Incomplete | None = None,
        beziers: bool = False,
    ): ...
    def remove_knots(
        self, parametric_dimension, knots, tolerance: Incomplete | None = None
    ): ...
    def normalize_knot_vectors(self) -> None: ...
    def extract_bezier_patches(self): ...

class BSpline(BSplineBase):
    def __init__(
        self,
        degrees: _npt.ArrayLike | None = None,
        knot_vectors: _npt.ArrayLike | None = None,
        control_points: _npt.ArrayLike | None = None,
        spline: Incomplete | None = None,
    ) -> None: ...
    @property
    def bspline(self): ...
    @property
    def nurbs(self): ...
