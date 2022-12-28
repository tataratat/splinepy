from splinepy.bspline import BSplineBase


class NURBS(BSplineBase):
    """
    Non-Uniform Rational B-Spline.
    """

    __slots__ = ()

    def __init__(
        self,
        degrees=None,
        knot_vectors=None,
        control_points=None,
        weights=None,
        spline=None,
    ):
        """
        NURBS.

        Parameters
        -----------
        degrees: (para_dim,) list-like
        knot_vectors: (para_dim, n) list
        control_points: (m, dim) list-like
        weights: (m,) list-like

        Returns
        --------
        None
        """
        super().__init__(
            spline=spline,
            degrees=degrees,
            knot_vectors=knot_vectors,
            control_points=control_points,
            weights=weights,
        )
