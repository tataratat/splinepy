from splinepy import settings, splinepy_core
from splinepy.bezier import BezierBase


class RationalBezier(BezierBase):
    r"""
    RationalBezier (Spline).

    Passes all arguments to :code:`super.__init__()`, see :class:`.BezierBase`.

    Rational Beziers are an extension to classical Beziers, similar like NURBS
    are an extension to B-Splines. By introducing a scalar weighting function
    :math:`W:\mathbb{R}^{N_{param}}\to\mathbb{R}^{+}_{*}`, we can construct
    modified (rational) basis functions, that allow to define rational
    Beziers.
    For a one-dimensional parameter space (:math:`N_{param}=1`), this
    weighting function reads

    .. math::
            W(u) = \sum_{i=0}^{l} B^{i;p}(u) w^i

    with scalar weights :math:`w^i\in\mathbb{R}^{+}_{*}`.
    Consequently, for the two-dimensional case (:math:`N_{param}=2`) we have

    .. math::
            W(u,v) = \sum_{i=0}^{l} \sum_{j=0}^{m} B^{i;p}(u) B^{j;q}(v)
                w^{i,j}

    and for the three-dimensional case (:math:`N_{param}=3`)

    .. math::
            W(u,v,w) = \sum_{i=0}^{l} \sum_{j=0}^{m} \sum_{k=0}^{n}
                B^{i;p}(u) B^{j;q}(v) B^{k;r}(w) w^{i,j,k}

    We proceed with a description of how to construct rational Beziers:

    .. note::
        For simplicity, we restrict ourselves to the three most common types
        of splines in the following, namely curves, surfaces and volumes,
        although :code:`splinepy` also supports higher dimensions, see
        the documentation of :class:`.Spline` for more information.

    1. A rational Bezier of degree :math:`p` with control points
    :math:`P^i\in\mathbb{R}^{N_{phys}}`, weights
    :math:`w^i\in\mathbb{R}^{+}_{*}`, and a one-dimensional parameter space
    (i.e., :math:`N_{param}=1`) corresponds to a line embedded into the
    physical space

    .. math::
            C(u) = \sum_{i=0}^{l} R^{i;p}(u) P^i

    with the modified (rational) basis functions

    .. math::
            R^{i;p}(u) = \frac{B^{i;p}(u) w^i}{W(u)}

    2. A rational Bezier of degrees :math:`p,q` with control points
    :math:`P^{i,j}\in\mathbb{R}^{N_{phys}}`, weights
    :math:`w^{i,j}\in\mathbb{R}^{+}_{*}`, and a two-dimensional parameter
    space (i.e., :math:`N_{param}=2`) corresponds to a surface embedded into
    the physical space

    .. math::
            S(u,v) = \sum_{i=0}^{l} \sum_{j=0}^{m} R^{i,j;p,q}(u,v) P^{i,j}

    with the modified (rational) basis functions

    .. math::
            R^{i,j;p,q}(u,v) =
                \frac{\tilde{N}^{i,j;p,q}(u,v) w^{i,j}}{W(u,v)}

    3. A rational Bezier of degrees :math:`p,q,r` with control points
    :math:`P^{i,j,k}\in\mathbb{R}^{N_{phys}}`, weights
    :math:`w^{i,j,k}\in\mathbb{R}^{+}_{*}`, and a three-dimensional parameter
    space (i.e., :math:`N_{param}=3`) corresponds to a volume embedded into
    the physical space

    .. math::
            V(u,v,w) = \sum_{i=0}^{l} \sum_{j=0}^{m} \sum_{k=0}^{n}
                R^{i,j,k;p,q,r}(u,v,w) P^{i,j,k}

    with the modified (rational) basis functions

    .. math::
            R^{i,j,k;p,q,r}(u,v,w) =
                \frac{
                    \tilde{N}^{i,j,k;p,q,r}(u,v,w) w^{i,j,k}
                }{
                    W(u,v,w)
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

    @property
    def rationalbezier(self):
        return self.copy(saved_data=False)

    @property
    def nurbs(self):
        return settings.NAME_TO_TYPE["NURBS"](
            spline=splinepy_core.same_spline_with_knot_vectors(self)
        )
