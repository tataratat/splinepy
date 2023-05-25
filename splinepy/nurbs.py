from splinepy.bspline import BSplineBase


class NURBS(BSplineBase):
    r"""
    Non-Uniform Rational B-Spline.

    NURBS are an extension of B-Splines overcoming their drawback of being
    unable to represent circular shapes. This is achieved by the introduction
    of a weighting function
    :math:`W^N:\mathbb{R}^{N_{param}}\to\mathbb{R}^{+}_{*}`, defined as (for
    a one-dimensional parameter space, i.e., :math:`N_{param}=1`)

    .. math::
            W^N(u) = \sum_{i=1}^{l} N_{i;p}(u) w_i

    with scalar weights :math:`w_i\in\mathbb{R}^{> 0}`. Note, that the same basis functions are used for both the weighting function and projection space.
    Consequently, for the two-dimensional case (:math:`N_{param}=2`) we have

    .. math::
            W^N(u,v) = \sum_{i=1}^{l} \sum_{j=1}^{m} N_{i;p}(u) N_{j;q}(v)
                w_{i,j}

    and for the three-dimensional case (:math:`N_{param}=3`)

    .. math::
            W^N(u,v,w) = \sum_{i=1}^{l} \sum_{j=1}^{m} \sum_{k=1}^{n}
                N_{i;p}(u) N_{j;q}(v) N_{k;r}(w) w_{i,j,k}

    We can now construct different spline embeddings similar to the
    description in :class:`.BSplineBase`, only including the additional
    weighting:

    1. A NURBS of degree :math:`p` with control points
    :math:`P_i\in\mathbb{R}^{N_{phys}}`, weights
    :math:`w_i\in\mathbb{R}^{> 0}`, and a one-dimensional parameter space
    (:math:`N_{param}=1`) corresponds to a line, embedded into the physical
    space:

    .. math::
            S^N(u) = \sum_{i=1}^{l} R^N_{i;p}(u) P_i

    with the modified (rational) basis functions

    .. math::
            R^N_{i;p}(u) = \frac{N_{i;p}(u)w_i}{W^N(u)}

    2. A NURBS of degrees :math:`p,q` with control points
    :math:`P_{i,j}\in\mathbb{R}^{N_{phys}}`, weights
    :math:`w_{i,j}\in\mathbb{R}^{> 0}`, and a two-dimensional parameter
    space (:math:`N_{param}=2`) corresponds to a surface, embedded into the
    physical space:

    .. math::
            S^N(u,v) = \sum_{i=1}^{l} \sum_{j=1}^{m} R^N_{i,j;p,q}(u,v) P_{i,j}

    with the modified (rational) basis functions

    .. math::
            R^N_{i,j;p,q}(u,v) =
                \frac{\tilde{N}_{i,j;p,q}(u,v) w_{i,j}}{W^N(u,v)}

    3. A NURBS of degrees :math:`p,q,r` with control points
    :math:`P_{i,j,k}\in\mathbb{R}^{N_{phys}}`, weights
    :math:`w_{i,j,k}\in\mathbb{R}^{> 0}`, and a three-dimensional parameter
    space (:math:`N_{param}=3`) corresponds to a volume, embedded into the
    physical space:

    .. math::
            S^N(u,v,w) = \sum_{i=1}^{l} \sum_{j=1}^{m} \sum_{k=1}^{n}
                R^N_{i,j,k;p,q,r}(u,v,w) P_{i,j,k}

    with the modified (rational) basis functions

    .. math::
            R^N_{i,j,k;p,q,r}(u,v,w) =
                \frac{
                    \tilde{N}_{i,j,k;p,q,r}(u,v,w) w_{i,j,k}
                }{
                    W^N(u,v,w)
                }
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
