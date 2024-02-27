from functools import wraps as _wraps

import numpy as _np

from splinepy._base import SplinepyBase as _SplinepyBase
from splinepy.utils.data import cartesian_product as _cartesian_product


def _get_integral_measure(spline):
    """
    Determines the appropriate measure to be used in integration

    If the spline dimension matches its parametric dimension it will return a
    Callable in the form

    .. math::
        \\mathcal{J}_S = det(\\mathbf(J))

    If the physical dimension is bigger then the paramtric dimension it will
    return

    .. math::
        \\mathcal{J}_S = \big( det(\\mathbf(J)^T \\mathbf(J)\big)^{0.5}

    Parameters
    ----------
    spline : Spline / Multipatch
      For parametric and physical dimension

    Returns
    -------
    measure : Callable
      single patch only
    """
    # Check dimensionality
    if spline.dim == spline.para_dim:

        def measure(spline_patch, positions):
            return _np.linalg.det(spline_patch.jacobian(positions))

        return measure

    elif spline.dim > spline.para_dim:

        def measure(spline_patch, positions):
            jacs = spline_patch.jacobian(positions)
            return _np.sqrt(
                _np.linalg.det(_np.einsum("oji,ojk->oik", jacs, jacs))
            )

        return measure

    else:
        raise ValueError("`Volume` not supported if para_dim > dim")


def _default_quadrature_orders(spline):
    expected_degree = spline.degrees * spline.para_dim + 1
    if spline.is_rational:
        expected_degree += 2
        spline._logd("Integration on rational spline is only approximation")

    # Gauss-Legendre is exact for polynomials 2*n-1
    return _np.ceil((expected_degree + 1) * 0.5).astype("int")


def _get_quadrature_information(spline, orders=None):
    """
    Select appropriate integration order (gauss-legendre)

    Determinante of a polynomial spline with para_dim==dim has degree

    .. math::
        p_i^{det} = n_{dim} \\cdot p_i - 1

    cf. [Mantzaflaris et al., 2017,
    DOI:http://dx.doi.org/10.1016/j.cma.2016.11.013]

    The same order approximation will also be used for the metric

    .. math::
        \\mathcal{J}_S = \big( det(\\mathbf(J)^T \\mathbf(J)\big)^{0.5}

    Parameters
    ----------
    spline : Spline
      Spline for integration
    orders : array-like (optional)
      Orders along every parametric dimension

    Returns
    -------
    positions : np.ndarray
      quadrature position in unit-square
    weights : np.ndarray
      quadrature weights
    """

    # Determine integration points
    positions = []
    weights = []

    # Determine integration orders
    if orders is None:
        quad_orders = _default_quadrature_orders(spline)

    else:
        quad_orders = _np.ascontiguousarray(orders, dtype=int).flatten()
        if quad_orders.size != spline.para_dim:
            raise ValueError(
                "Integration order must be array of size para_dim"
            )

    for order in quad_orders:
        # Get legendre quadratuer points
        pos, ws = _np.polynomial.legendre.leggauss(order)
        # from (-1,1) to (0,1)
        positions.append(pos * 0.5 + 0.5)
        weights.append(ws * 0.5)

    # summarize integration points
    return (
        _cartesian_product(positions),
        _cartesian_product(weights).prod(axis=1),
    )


def volume(spline, orders=None):
    r"""Compute volume of a given spline

    Parameters
    ----------
    spline : Spline
        (self if called via integrator)
      splinepy - spline type
    orders : array-like (optional)
      order for gauss quadrature

    Returns
    -------
    volume : float
      Integral of dim-dimensional object
    """
    from splinepy.spline import Spline as _Spline

    # Check i_nput type
    if not isinstance(spline, _Spline):
        raise NotImplementedError("integration only works for splines")

    # Retrieve aux info
    meas = _get_integral_measure(spline)
    positions, weights = _get_quadrature_information(spline, orders)

    # Calculate Volume
    if spline.has_knot_vectors:
        volume = _np.sum(
            [
                _np.sum(meas(b, positions) * weights)
                for b in spline.extract.beziers()
            ]
        )
    else:
        volume = _np.sum(meas(spline, positions) * weights)
    return volume


def parametric_function(
    spline,
    function,
    orders=None,
):
    """Integrate a function defined within the parametric domain

    Parameters
    ----------
    spline : Spline
        (self if called via integrator)
    function : Callable
    orders : optional

    Returns
    -------
    integral : np.ndarray
    """
    from splinepy.spline import Spline as _Spline

    # Check i_nput type
    if not isinstance(spline, _Spline):
        raise NotImplementedError("integration only works for splines")

    # Retrieve aux info
    meas = _get_integral_measure(spline)
    positions, weights = _get_quadrature_information(spline, orders)

    # Calculate Volume
    if spline.has_knot_vectors:
        # positions must be mapped into each knot-element
        para_view = spline.create.parametric_view(axes=False)

        # get initial shape
        initial = function([positions[0]])
        result = _np.zeros(initial.shape[1])
        for bezier_element in para_view.extract.beziers():
            quad_positions = bezier_element.evaluate(positions)
            result += _np.einsum(
                "i...,i,i->...",
                function(quad_positions),
                meas(spline, quad_positions),
                weights,
                optimize=True,
            )

    else:
        result = _np.einsum(
            "i...,i,i->...",
            function(positions),
            meas(spline, positions),
            weights,
            optimize=True,
        )
    return result


def physical_function(
    function,  # noqa ARG001
    orders,  # noqa ARG001
):
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

    __slots__ = ("_helpee",)

    def __init__(self, spl):
        self._helpee = spl

    @_wraps(volume)
    def volume(self, *args, **kwargs):
        return volume(self._helpee, *args, **kwargs)

    @_wraps(parametric_function)
    def parametric_function(self, *args, **kwargs):
        return parametric_function(self._helpee, *args, **kwargs)


class FieldIntegrator(_SplinepyBase):
    __slots__ = (
        "_positions",
        "_weights",
        "_helpee",
        "_orders",
        "_supports",
        "_bezier_patches",
        "_global_positions",
        "_jacobians",
        "_jacobian_weights",
    )

    def __init__(self, spline, orders=None):
        """ """
        self._helpee = spline

        self.reset(orders)

    def reset(self, orders=None):
        """ """
        if orders is None:
            self._orders = _default_quadrature_orders(self._helpee)
        else:
            self._orders = orders

        self._positions, self._weights = _get_quadrature_information(
            self._helpee, self._orders
        )

        self._bezier_patches = None
        self._global_positions = None
        self._supports = None
        self._jacobians = None
        self._jacobian_weights = None

    @property
    def positions(self):
        """
        Normalized Quadrature positions. Can use this value for
        """
        return self._positions

    @property
    def global_positions(self):
        """
        Quadrature points in global position
        """
        if self._global_positions is not None:
            return self._global_positions

        # TODO: clamped knot vector check once it's merged
        lower_bounds_per_dim = []
        span_scales_per_dim = []
        for ukv in self._helpee.unique_knots:
            lower_bounds_per_dim.append(ukv[:-1])
            span_scales_per_dim.append(_np.diff(ukv))
        lower_bounds = _cartesian_product(lower_bounds_per_dim, reverse=True)
        span_scales = _cartesian_product(span_scales_per_dim, reverse=True)

        # add lower bound as offsets using np.broadcast rules
        n_quads, dim = self.positions.shape
        n_elem = len(lower_bounds)

        # create normalized quad points for each element
        self._global_positions = _np.tile(self.positions, (n_elem, 1)).reshape(
            n_elem, n_quads, dim
        )
        # scale them
        self._global_positions *= span_scales.reshape(n_elem, 1, dim)
        # apply offset
        self._global_positions += lower_bounds.reshape(n_elem, 1, dim)

        return self._global_positions

    @property
    def supports(self):
        """ """
        if self._supports is not None:
            return self._supports

        self._supports = self._helpee.supports(self.global_positions)
        return self._supports

    @property
    def quadrature_weights(self):
        """ """
        return self._weights

    @property
    def bezier_patches(self):
        """Returns extracted bezier patches.
        These beziers can be used to compute jacobians and weights thereof.
        """
        if self._bezier_patches is not None:
            return self._bezier_patches

        self._bezier_patches = self._helpee.extract.beziers()
        return self._bezier_patches

    @property
    def jacobians(self):
        if self._jacobians is not None:
            return self._jacobians

        self._jacobians = [
            bp.jacobian(self.positions) for bp in self.bezier_patches
        ]

        return self._jacobians

    @property
    def jacobian_inverses(self):
        if self._jacobian_inverses is not None:
            return self._jacobian_inverses

        # will raise if this isn't square
        self._jacobian_inverses = [_np.linalg.inv(j) for j in self.jacobians]

        return self._jacobian_inverses

    @property
    def jacobian_weights(self):
        """
        Jacobian determinants if jacobian is a square matrix.
        Otherwise,
        """
        if self._jacobian_weights is not None:
            return self._jacobian_weights

        meas = _get_integral_measure(self._helpee)
        self._jacobian_weights = [
            meas(bp, self.positions) for bp in self.bezier_patches
        ]

        return self._jacobian_weights

    @property
    def assemble_matrix(self, function, dim, matrix=None):
        """
        """
        pass

    @property
    def assemble_vector(self, function, dim, matrix=None):
        """
        """
        pass