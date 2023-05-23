from splinepy.bspline import BSplineBase


class NURBS(BSplineBase):
    r"""
    Non-Uniform Rational B-Spline.

    NURBS are an extension of B-Splines overcoming their drawback of being 
    unable to represent circular shapes. This is achieved by the introduction 
    of a weighting function :math:`W:\\mathbb{R}^{N_{param}}\\to\\mathbb{R}`, 
    defined as (for a one-dimensional parameter space)
    
    .. math::
            W(u^1) = \\sum_{i=1}^{n^1} N_{i;p^1}^{\\Theta^1}(u^1) w_i

    with scalar weights :math:`w_i\\in\\mathbb{R}`.
    Consequently, for two-dimensional case, we have

    .. math::
            W(u^1,u^2) = \\sum_{i=1}^{n^1} \\sum_{j=1}^{n^2} 
                N_{i;p^1}^{\\Theta^1}(u^1) N_{j;p^2}^{\\Theta^2}(u^2) w_{i,j}

    and for the three-dimensional case

    .. math::
            W(u^1,u^2,u^3) = \\sum_{i=1}^{n^1} \\sum_{j=1}^{n^2} 
                \\sum_{k=1}^{n^3} N_{i;p^1}^{\\Theta^1}(u^1) 
                N_{j;p^2}^{\\Theta^2}(u^2) N_{j;p^3}^{\\Theta^3}(u^3) w_{i,j,k}
    
    We can now construct different spline embeddings similar to the
    description in :class:`.BSplineBase` just including the additional 
    weighting:

    #. A NURBS of degree :math:`p^1` with control points 
    :math:`P_i\\in\\mathbb{R}^{N_{phys}`, weights :math:`w_i\\in\\mathbb{R}`,
    and a one-dimensional parameter space corresponds to a line, embedded into 
    the physical space:

    .. math:: 
            S(u^1) = \\sum_{i=1}^{n^1} R_{i;p^1}(u^1) P_i

    with the modified (rational) basis functions

    .. math::
            R_{i;p^1}(u^1) = \\frac{N_{i;p^1}^{\\Theta^1}(u^1)w_i}{W(u^1)}

    #. A NURBS of degrees :math:`p^1,p^2` with control points 
    :math:`P_{i,j}\\in\\mathbb{R}^{N_{phys}`, weights 
    :math:`w_{i,j}\\in\\mathbb{R}`, and a two-dimensional parameter space 
    corresponds to a line, embedded into the physical space:

    .. math:: 
            S(u^1,u^2) = \\sum_{i=1}^{n^1} \\sum_{j=1}^{n^2} 
                R_{i,j;p^1,p^2}(u^1,u^2) P_{i,j}

    with the modified (rational) basis functions

    .. math::
            R_{i,j;p^1,p^2}(u^1,u^2) = 
                \\frac{\\tilde{N}_{i,j;p^1,p^2}(u^1,u^2)w_{i,j}}{W(u^1,u^2)}

    #. A NURBS of degrees :math:`p^1,p^2,p^3` with control points 
    :math:`P_{i,j,k}\\in\\mathbb{R}^{N_{phys}`, weights 
    :math:`w_{i,j,k}\\in\\mathbb{R}`, and a three-dimensional parameter space 
    corresponds to a line, embedded into the physical space:

    .. math:: 
            S(u^1,u^2,u^3) = \\sum_{i=1}^{n^1} \\sum_{j=1}^{n^2} 
                \\sum_{k=1}^{n^3} R_{i,j,k;p^1,p^2,p^3}(u^1,u^2,u^3) P_{i,j,k}

    with the modified (rational) basis functions

    .. math::
            R_{i,j,k;p^1,p^2,p^3}(u^1,u^2,u^3) = 
                \\frac{
                    \\tilde{N}_{i,j,k;p^1,p^2,p^3}(u^1,u^2,u^3)w_{i,j,k}
                }{
                    W(u^1,u^2,u^3)
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
