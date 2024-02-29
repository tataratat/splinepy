from _typeshed import Incomplete

from splinepy.bspline import BSplineBase as _BSplineBase

class NURBS(_BSplineBase):
    def __init__(
        self,
        degrees: Incomplete | None = None,
        knot_vectors: Incomplete | None = None,
        control_points: Incomplete | None = None,
        weights: Incomplete | None = None,
        spline: Incomplete | None = None,
    ) -> None: ...
    @property
    def nurbs(self): ...
