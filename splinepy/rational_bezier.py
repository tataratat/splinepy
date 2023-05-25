from splinepy.bezier import BezierBase


class RationalBezier(BezierBase):
    r"""
    RationalBezier (Spline).

    Passes all arguments to :code:`super.__init__()`, see :class:`.BezierBase`.

    Rational Beziers are an extension to classical Beziers, similar like NURBS
    are an extension to B-Splines. By introducing a scalar weighting function
    :math:`W^B:\mathbb{R}^{N_{param}}\to\mathbb{R}^{+}_{*}`, we can construct
    modified (rational) basis functions, that allow to define rational
    Beziers.
    For a one-dimensional parameter space (:math:`N_{param}=1`), this
    weighting function reads

    .. math::
            W^B(u) = \sum_{i=0}^{l} B_{i;p}(u) w_i

    with scalar weights :math:`w_i\in\mathbb{R}^{+}_{*}`.
    Consequently, for the two-dimensional case (:math:`N_{param}=2`) we have

    .. math::
            W^B(u,v) = \sum_{i=0}^{l} \sum_{j=0}^{m} B_{i;p}(u) B_{j;q}(v)
                w_{i,j}

    and for the three-dimensional case (:math:`N_{param}=3`)

    .. math::
            W^B(u,v,w) = \sum_{i=0}^{l} \sum_{j=0}^{m} \sum_{k=0}^{n}
                B_{i;p}(u) B_{j;q}(v) B_{k;r}(w) w_{i,j,k}

    We proceed with a description of how to construct rational Beziers:

    .. note::
        For simplicity, we restrict ourselves to the three most common types 
        of splines in the following, namely curves, surfaces and volumes,
        although :code:`splinepy` also supports higher dimensions, see
        the documentation of :class:`.Spline` for more information.

    1. A rational Bezier of degree :math:`p` with control points
    :math:`P_i\in\mathbb{R}^{N_{phys}}`, weights
    :math:`w_i\in\mathbb{R}^{+}_{*}`, and a one-dimensional parameter space
    (i.e., :math:`N_{param}=1`) corresponds to a line embedded into the
    physical space

    .. math::
            C^B(u) = \sum_{i=0}^{l} R^B_{i;p}(u) P_i

    with the modified (rational) basis functions

    .. math::
            R^B_{i;p}(u) = \frac{B_{i;p}(u) w_i}{W^B(u)}

    2. A rational Bezier of degrees :math:`p,q` with control points
    :math:`P_{i,j}\in\mathbb{R}^{N_{phys}}`, weights
    :math:`w_{i,j}\in\mathbb{R}^{+}_{*}`, and a two-dimensional parameter
    space (i.e., :math:`N_{param}=2`) corresponds to a surface embedded into
    the physical space

    .. math::
            S^B(u,v) = \sum_{i=0}^{l} \sum_{j=0}^{m} R^B_{i,j;p,q}(u,v) P_{i,j}

    with the modified (rational) basis functions

    .. math::
            R^B_{i,j;p,q}(u,v) =
                \frac{\tilde{N}_{i,j;p,q}(u,v) w_{i,j}}{W^B(u,v)}

    3. A rational Bezier of degrees :math:`p,q,r` with control points
    :math:`P_{i,j,k}\in\mathbb{R}^{N_{phys}}`, weights
    :math:`w_{i,j,k}\in\mathbb{R}^{+}_{*}`, and a three-dimensional parameter
    space (i.e., :math:`N_{param}=3`) corresponds to a volume embedded into
    the physical space

    .. math::
            V^B(u,v,w) = \sum_{i=0}^{l} \sum_{j=0}^{m} \sum_{k=0}^{n}
                R^B_{i,j,k;p,q,r}(u,v,w) P_{i,j,k}

    with the modified (rational) basis functions

    .. math::
            R^B_{i,j,k;p,q,r}(u,v,w) =
                \frac{
                    \tilde{N}_{i,j,k;p,q,r}(u,v,w) w_{i,j,k}
                }{
                    W^B(u,v,w)
                }

    Higher-dimensional instances are constructed accordingly.

    **Usage**:

    .. code-block:: python

        # Rational bezier curve (quarter arc)
        rational_bezier = splinepy.RationalBezier(
            degrees=[2],
            control_points=[
                [0.0, 0.0],
                [1.0, 0.0],
                [1.0, 1.0]
            ],
            weights=[1.0, 2.0**-0.5, 1.0]
        )
    
    Parameters
    -----------
    degrees: (para_dim,) array-like
    control_points: (m, dim) array-like
    weights: (m) array-like
    spline: Spline

    Returns
    --------
    None
    """

    __slots__ = ()

    def __init__(
        self, degrees=None, control_points=None, weights=None, spline=None
    ):
        super().__init__(
            spline=spline,
            degrees=degrees,
            control_points=control_points,
            weights=weights,
        )
