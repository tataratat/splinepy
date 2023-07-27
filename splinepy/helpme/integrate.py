import numpy as np

from splinepy.utils.data import cartesian_product


def volume(spline, orders=None):
    r"""Compute volume of a given spline

    Determinante has degree

    .. math::
        p_i^{det} = n_{dim} \cdot p_i - 1

    cf. [Mantzaflaris et al., 2017,
    DOI:http://dx.doi.org/10.1016/j.cma.2016.11.013]

    Parameters
    ----------
    spline : Spline to be integrated
      splinepy - spline type
    orders : array-like (optional)
      order for gauss quadrature

    Returns
    -------
    volume : float
      Integral of dim-dimensional object
    """
    from splinepy.spline import Spline

    # Check input type
    if not isinstance(spline, Spline):
        raise NotImplementedError("Extrude only works for splines")

    # Determine integration points
    positions = []
    weights = []

    # Check dimensionality
    if not (spline.dim == spline.para_dim):
        raise ValueError(
            "`Volume` of embedded spline depends on projection, "
            "integration is aborted"
        )

    # Determine integration orders
    if orders is None:
        expected_degree = spline.degrees * spline.para_dim + 1
        if spline.is_rational:
            expected_degree += 2
            spline._logw(
                "Integration on rational spline is only approximation"
            )

        # Gauss-Legendre is exact for polynomials 2*n-1
        quad_orders = np.ceil((expected_degree + 1) * 0.5).astype("int")
    else:
        quad_orders = np.ascontiguousarray(orders, dtype=int).flatten()
        if quad_orders.size != spline.para_dim:
            raise ValueError(
                "Integration order must be array of size para_dim"
            )

    for order in quad_orders:
        # Get legendre quadratuer points
        pos, ws = np.polynomial.legendre.leggauss(order)
        # from (-1,1) to (0,1)
        positions.append(pos * 0.5 + 0.5)
        weights.append(ws * 0.5)

    # summarize integration points
    positions = cartesian_product(positions)
    weights = cartesian_product(weights).prod(axis=1)

    # Calculate Volume
    if spline.has_knot_vectors:
        volume = np.sum(
            [
                np.sum(np.linalg.det(b.jacobian(positions)) * weights)
                for b in spline.extract.beziers()
            ]
        )
    else:
        volume = np.sum(np.linalg.det(spline.jacobian(positions)) * weights)
    return volume


def parametric_function(function, orders):
    """Integrate a function defined within the parametric domain"""
    raise NotImplementedError(
        "Function not implemented yet. Please feel free to write an issue, if "
        "you need it: github.com/tatarata/splinepy/issues"
    )


def physical_function(function, orders):
    """Integrate a function defined within the physical domain"""
    raise NotImplementedError(
        "Function not implemented yet. Please feel free to write an issue, if "
        "you need it: github.com/tatarata/splinepy/issues"
    )


class Integrator:
    """Helper class to integrate some values on a given spline geometry

    Examples
    --------

    .. code-block:: python

      splinepy.helpme.integrate.volume()
      # Equivalent to
      spline.integrate.volume()

    Parameters
    ----------
    spline : Spline
      Spline parent
    """

    def __init__(self, spl):
        self.spline = spl

    def volume(self, *args, **kwargs):
        return volume(self.spline, *args, **kwargs)
