from splinepy.bspline import BSplineBase


class NURBS(BSplineBase):
    r"""
    Non-Uniform Rational B-Spline.

    Passes all arguments to :code:`super.__init__()`,
    see :class:`.BSplineBase`.

    NURBS are an extension of B-Splines overcoming their drawback of being
    unable to represent circular shapes. This is achieved by the introduction
    of a weighting function
    :math:`W:\mathbb{R}^{N_{param}}\to\mathbb{R}^{+}_{*}`, defined as (for a
    one-dimensional parameter space, i.e., :math:`N_{param}=1`)

    .. math::
            W(u) = \sum_{i=1}^{l} N^{i;p}(u) w^i

    with scalar weights :math:`w^i\in\mathbb{R}^{+}_{*}`. Note, that the same
    basis functions are used for both the weighting function and projection
    space. Consequently, for the two-dimensional case (:math:`N_{param}=2`) we
    have

    .. math::
            W(u,v) = \sum_{i=1}^{l} \sum_{j=1}^{m} N^{i;p}(u) N^{j;q}(v)
                w^{i,j}

    and for the three-dimensional case (:math:`N_{param}=3`)

    .. math::
            W(u,v,w) = \sum_{i=1}^{l} \sum_{j=1}^{m} \sum_{k=1}^{n}
                N^{i;p}(u) N^{j;q}(v) N^{k;r}(w) w^{i,j,k}

    We can now construct different spline embeddings similar to the
    description in :class:`.BSplineBase`, only including the additional
    weighting:

    .. note::
        For simplicity, we restrict ourselves to the three most common types
        of splines in the following, namely curves, surfaces and volumes,
        although :code:`splinepy` also supports higher dimensions, see
        the documentation of :class:`.Spline` for more information.

    1. A NURBS of degree :math:`p` with control points
    :math:`P^i\in\mathbb{R}^{N_{phys}}`, weights
    :math:`w^i\in\mathbb{R}^{+}_{*}`, and a one-dimensional parameter space
    (:math:`N_{param}=1`) corresponds to a line, embedded into the physical
    space:

    .. math::
            C(u) = \sum_{i=1}^{l} R^{i;p}(u) P^i

    with the modified (rational) basis functions

    .. math::
            R^{i;p}(u) = \frac{N^{i;p}(u)w^i}{W(u)}

    2. A NURBS of degrees :math:`p,q` with control points
    :math:`P^{i,j}\in\mathbb{R}^{N_{phys}}`, weights
    :math:`w^{i,j}\in\mathbb{R}^{+}_{*}`, and a two-dimensional parameter
    space (:math:`N_{param}=2`) corresponds to a surface, embedded into the
    physical space:

    .. math::
            S(u,v) = \sum_{i=1}^{l} \sum_{j=1}^{m} R^{i,j;p,q}(u,v) P^{i,j}

    with the modified (rational) basis functions

    .. math::
            R^{i,j;p,q}(u,v) =
                \frac{\tilde{N}^{i,j;p,q}(u,v) w^{i,j}}{W(u,v)}

    3. A NURBS of degrees :math:`p,q,r` with control points
    :math:`P^{i,j,k}\in\mathbb{R}^{N_{phys}}`, weights
    :math:`w^{i,j,k}\in\mathbb{R}^{+}_{*}`, and a three-dimensional parameter
    space (:math:`N_{param}=3`) corresponds to a volume, embedded into the
    physical space:

    .. math::
            V(u,v,w) = \sum_{i=1}^{l} \sum_{j=1}^{m} \sum_{k=1}^{n}
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

        # NURBS curve (circle)
        nurbs_curve = splinepy.NURBS(
            degrees=[2],
            knot_vectors=[0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, \
                1.0, 1.0, 1.0],
            control_points=[
                [ 1.0,  0.0],
                [ 1.0,  1.0],
                [ 0.0,  1.0],
                [-1.0,  1.0],
                [-1.0,  0.0],
                [-1.0, -1.0],
                [ 0.0, -1.0],
                [ 1.0, -1.0],
                [ 1.0,  0.0]
            ]
        )

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

    __slots__ = ()

    def __init__(
        self,
        degrees=None,
        knot_vectors=None,
        control_points=None,
        weights=None,
        spline=None,
    ):
        super().__init__(
            spline=spline,
            degrees=degrees,
            knot_vectors=knot_vectors,
            control_points=control_points,
            weights=weights,
        )

    @property
    def nurbs(self):
        return self.copy(saved_data=False)
