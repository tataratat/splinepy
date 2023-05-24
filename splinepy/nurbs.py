from splinepy.bspline import BSplineBase


class NURBS(BSplineBase):
    r"""
    Non-Uniform Rational B-Spline.

    NURBS are an extension of B-Splines overcoming their drawback of being
    unable to represent circular shapes. This is achieved by the introduction
    of a weighting function :math:`W:\\mathbb{R}^{N_{param}}\\to\\mathbb{R}`,
    defined as (for a one-dimensional parameter space)

    .. math::
            W(u) = \\sum_{i=1}^{l} N_{i;p}(u) w_i

    with scalar weights :math:`w_i\\in\\mathbb{R}`.
    Consequently, for two-dimensional case, we have

    .. math::
            W(u,v) = \\sum_{i=1}^{l} \\sum_{j=1}^{m} N_{i;p}(u) N_{j;q}(v) 
                w_{i,j}

    and for the three-dimensional case

    .. math::
            W(u,v,w) = \\sum_{i=1}^{l} \\sum_{j=1}^{m} \\sum_{k=1}^{n} 
                N_{i;p}(u) N_{j;q}(v) N_{k;r}(w) w_{i,j,k}

    We can now construct different spline embeddings similar to the
    description in :class:`.BSplineBase`, only including the additional
    weighting:

    #. A NURBS of degree :math:`p` with control points
    :math:`P_i\\in\\mathbb{R}^{N_{phys}`, weights :math:`w_i\\in\\mathbb{R}`,
    and a one-dimensional parameter space corresponds to a line, embedded into
    the physical space:

    .. math::
            S(u) = \\sum_{i=1}^{l} R_{i;p}(u) P_i

    with the modified (rational) basis functions

    .. math::
            R_{i;p}(u) = \\frac{N_{i;p}(u)w_i}{W(u)}

    #. A NURBS of degrees :math:`p,q` with control points
    :math:`P_{i,j}\\in\\mathbb{R}^{N_{phys}`, weights
    :math:`w_{i,j}\\in\\mathbb{R}`, and a two-dimensional parameter space
    corresponds to a line, embedded into the physical space:

    .. math::
            S(u,v) = \\sum_{i=1}^{l} \\sum_{j=1}^{m} R_{i,j;p,q}(u,v) P_{i,j}

    with the modified (rational) basis functions

    .. math::
            R_{i,j;p,q}(u,v) =
                \\frac{\\tilde{N}_{i,j;p,q}(u,v) w_{i,j}}{W(u,v)}

    #. A NURBS of degrees :math:`p,q,r` with control points
    :math:`P_{i,j,k}\\in\\mathbb{R}^{N_{phys}`, weights
    :math:`w_{i,j,k}\\in\\mathbb{R}`, and a three-dimensional parameter space
    corresponds to a line, embedded into the physical space:

    .. math::
            S(u,v,w) = \\sum_{i=1}^{l} \\sum_{j=1}^{m} \\sum_{k=1}^{n} 
                R_{i,j,k;p,q,r}(u,v,w) P_{i,j,k}

    with the modified (rational) basis functions

    .. math::
            R_{i,j,k;p,q,r}(u,v,w) =
                \\frac{
                    \\tilde{N}_{i,j,k;p,q,r}(u,v,w) w_{i,j,k}
                }{
                    W(u,v,w)
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
