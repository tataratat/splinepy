import numpy as _np

from splinepy import settings as _settings
from splinepy.helpme.multi_index import MultiIndex
from splinepy.utils import log as _log
from splinepy.utils.data import has_scipy as _has_scipy
from splinepy.utils.data import make_matrix as _make_matrix

if _has_scipy:
    from scipy.sparse.linalg import spsolve as _spsolve


def parameterize(fitting_points, size, centripetal):
    """
    Parameterizes given fitting_points later used
    in interpolation/approximation.
    Values are then used to build adequate knot vectors.

    Parameters
    ----------
    fitting_points:(n, dim) array-like
      points to be interpolated/approximated
    size: array-like
      number of points per parametric dimension
    centripetal: bool
      if True -> centripetal parameterization will be used

    Returns
    -------
    parametric_coordinates: list<np.ndarray>
      array containing parameterization values per parametric dimension
    """

    def parametrize_to_line(reorganized_queries, axis, centripetal):
        """parameterize_line all at once.
        Expects with-MultiIndex-Reorganized points"""
        n_points = reorganized_queries.shape[axis]
        u_k = _np.empty(n_points)

        # norm from axis -1 as last axis is dim.
        chord_lengths = _np.linalg.norm(
            _np.diff(reorganized_queries, axis=axis), axis=-1
        )

        if centripetal:
            chord_lengths = _np.sqrt(chord_lengths)

        # indices for total chord_reshapes
        chord_lengths_shape = list(chord_lengths.shape)
        chord_lengths_shape[axis] = 1

        total_chord_length = _np.sum(chord_lengths, axis=axis).reshape(
            *chord_lengths_shape
        )
        chord_lengths /= total_chord_length

        # Compute cumulated sum
        if chord_lengths.ndim == 1:
            u_k[1:] = _np.cumsum(chord_lengths, axis=axis)
        else:
            # Indices for slices
            mean_slicer = list(range(fitting_points.ndim))
            mean_slicer.pop(axis)
            mean_slicer = tuple(mean_slicer)
            # Higher dimensional (surface ++)
            u_k[1:] = _np.mean(
                _np.cumsum(chord_lengths, axis=axis),
                axis=mean_slicer,
            )
        u_k[0] = 0.0
        # Numerical inaccuracies where an issue with the last value sometimes
        u_k[-1] = 1.0

        return u_k.reshape(-1, 1)

    mi = MultiIndex(size)

    # Loop over all dimensions and append each para_coords to list
    reorganized = fitting_points[mi._raveled_indices]
    parametric_coordinates = []
    for k in range(len(size)):
        parametric_coordinates.append(
            parametrize_to_line(reorganized, k, centripetal)
        )

    return parametric_coordinates


def compute_knot_vector(degree, n_control_points, u_k, n_fitting_points):
    """
    Computes the knot_vector for a spline.

    Parameters
    ----------
    degree: int
      degree of spline
    n_control_points: int
      number of control points
    u_k: (n, 1) array-like
      array containing parameterization values
    n_fitting_points: int
      number of fitting_points

    Returns
    -------
    knot_vector: (degree + n_control_points + 1, 1) np.ndarray
      array containing knots of spline
    """
    # initialize with m = n + p + 1
    knot_vector = _np.empty(degree + n_control_points + 1)

    knot_vector[: (degree + 1)] = 0.0
    knot_vector[-(degree + 1) :] = 1.0

    # check if interpolation or approximation
    if n_fitting_points == n_control_points:
        # interpolation
        # NURBS Book eq. 9.8 (n_fitting_points = m + 1,
        # n_control_points = n + 1)
        knot_vector[(degree + 1) : -(degree + 1)] = (
            _np.convolve(u_k.ravel(), _np.ones(degree), "valid")[1:-1] / degree
        )

    else:
        # approximation
        # NURBS Book eq. (9.68) (n_fitting_points = m + 1,
        # n_control_points = n + 1)
        d = (n_fitting_points) / (n_control_points - degree)

        u_k_flat = u_k.flat

        # NURBS Book eq. (9.69)
        for j in range(1, n_control_points - degree):
            i = int(j * d)
            alpha = (j * d) - i
            knot_vector[j + degree] = (1 - alpha) * u_k_flat[
                i - 1
            ] + alpha * u_k_flat[i]
    return knot_vector


def solve_for_control_points(
    fitting_points, fitting_spline, queries, interpolate_endpoints
):
    """
    Builds and solves the system to calculate control points.

    Parameters
    ----------
    fitting_points:(m + 1 x dim) array
        points to be interpolated/approximated
    fitting_spline:spline
        spline used for interpolation/approximation
    queries:(n x 1) array
        parametric values where the spline will be evaluated
    interpolate_endpoints:bool
        if True -> endpoints are interpolated

    Returns
    -------
    residual:float
        residual (coefficient_matrix @ control_points - fitting_points)
    """
    # build matrix
    coefficient_matrix = _make_matrix(
        *fitting_spline.basis_and_support(queries),
        fitting_spline.control_points.shape[0],
    )

    if (
        interpolate_endpoints
        and fitting_spline.control_points.shape != fitting_points.shape
    ):
        # reduction of matrix necessary for interpolation of endpoints
        # copy fitting_points to keep original fitting_points unchanged
        rhs = fitting_points.copy()

        fitting_spline.control_points[0] = rhs[0]
        fitting_spline.control_points[-1] = rhs[-1]

        if fitting_spline.control_points.shape[0] == 2:
            # control points equal to endpoints (straight line)
            residual = _np.linalg.norm(
                coefficient_matrix @ fitting_spline.control_points
                - fitting_points
            )
        elif _has_scipy:
            # move known values to RHS
            rhs[1:-1, :] -= (
                coefficient_matrix[1:-1, [0]] * fitting_points[0, :]
                + coefficient_matrix[1:-1, [-1]] * fitting_points[-1, :]
            ).todense()

            # solve system A^TAx = A^Tb
            # (tried using scipy.sparse.linalg.lsqr but had worse performance)
            at = coefficient_matrix[1:-1, 1:-1].transpose(copy=True)
            fitting_spline.control_points[1:-1] = _spsolve(
                at @ coefficient_matrix[1:-1, 1:-1], at @ rhs[1:-1]
            )
        else:
            # move known values to RHS
            rhs[1:-1, :] += (
                -coefficient_matrix[1:-1, [0]] * fitting_points[0, :]
                - coefficient_matrix[1:-1, [-1]] * fitting_points[-1, :]
            )

            # solve system using lstsq

            fitting_spline.control_points[1:-1] = _np.linalg.lstsq(
                coefficient_matrix[1:-1, 1:-1], rhs[1:-1]
            )[0]

    # interpolation or approximation with original system
    elif _has_scipy:
        at = coefficient_matrix.transpose(copy=True)
        fitting_spline.control_points = _spsolve(
            at @ coefficient_matrix, at @ fitting_points
        )

    else:
        fitting_spline.control_points = _np.linalg.lstsq(
            coefficient_matrix, fitting_points
        )[0]

    residual = _np.linalg.norm(
        coefficient_matrix @ fitting_spline.control_points - fitting_points
    )
    return residual


def _validate_specifications(n_query, degree, n_control_points, knot_vector):
    """Given main values of fitting, validates degree and n_control_points.
    If there's mismatch between expected value and given value, values will be
    overwritten."""
    # sanity checks for given values
    if degree is not None and knot_vector is not None:
        expected_ncps = len(knot_vector) - degree - 1
        if n_control_points is not None and n_control_points != expected_ncps:
            _log.error(
                f"n_control_points should be {expected_ncps}. Overwriting."
            )

        n_control_points = expected_ncps

    if knot_vector is not None and n_control_points is not None:
        expected_degree = len(knot_vector) - n_control_points - 1
        if degree is not None and degree != expected_degree:
            _log.error(f"degree should be {expected_degree}. Overwriting.")

        degree = expected_degree

    # we need degree
    if degree is None:
        raise ValueError(
            "Not enough input to determine degree. Please set degree."
        )

    # no n_control point -> same as query and this will be interpolation
    if n_control_points is None:
        _log.debug(
            "n_control_points not specified. "
            "Setting it the same as the number of queries"
        )
        n_control_points = n_query

    # degree check
    if degree >= n_control_points:
        raise ValueError("Given degree must be lower than n_control_points!")

    # check number of fitting_points and control points
    if n_control_points > n_query:
        raise ValueError(
            f"(n_fitting_points({n_query}) >= "
            f"n_control_points({n_control_points})) not satisfied!"
        )

    return degree, n_control_points


def _prepare_default_bspline1pd(
    n_queries, dim, degree, n_control_points, knot_vector, u_k
):
    # validate values. may raise if there's conflict.
    degree, n_control_points = _validate_specifications(
        n_queries, degree, n_control_points, knot_vector
    )

    # empty init of default fitting spline
    fitting_spline = _settings.NAME_TO_TYPE["BSpline"](
        [degree],
        [
            compute_knot_vector(
                degree=degree,
                n_control_points=n_control_points,
                u_k=u_k,
                n_fitting_points=n_queries,
            )
        ],
        _np.empty((n_control_points, dim)),
    )

    return fitting_spline


def curve(
    fitting_points,
    degree=None,
    n_control_points=None,
    knot_vector=None,
    fitting_spline=None,
    associated_queries=None,
    centripetal=True,
    interpolate_endpoints=True,
    verbose_output=False,
):
    """
    Fits a spline with given parameters through given fitting_points.
    Spline will be interpolated if n_control_points = n_fitting_points and
    approximated if not.

    Parameters
    ----------
    fitting_points: (m, dim) array
        points to be interpolated/approximated
    degree: int
        degree of spline
    n_control_points: int
        number of control points
    knot_vector: list
        desired knot vector of spline
    fitting_spline: Spline
        spline used for interpolation/approximation
    associated_queries: (n, 1) np.ndarray
        values where the spline will be evaluated
    centripetal: bool (default = True)
        if True -> centripetal parameterization will be used
    interpolate_endpoints: bool
        if True -> endpoints are interpolated
    verbose_outputs: dict
        returns additional information as dict

    Returns
    -------
    fitted_spline: Spline
        interpolated/approximated spline
    residual: float
        residual (coefficient_matrix @ control_points - fitting_points)
    """
    fitting_points = _np.asanyarray(fitting_points)

    # calculate n_fitting_points due to multiple usage
    n_fitting_points = fitting_points.shape[0]

    # determine evaluation points to build a linear system
    if associated_queries is not None:
        u_k = _np.asanyarray(associated_queries)
        if u_k.shape[1] != 1:
            raise ValueError("Associated queries need to have (-1, 1) shape")

    else:
        u_k = parameterize(
            fitting_points=fitting_points,
            size=[n_fitting_points],
            centripetal=centripetal,
        )[0]

    # create fitting_spline
    if fitting_spline is not None:
        # spline dimension check
        if fitting_spline.para_dim != 1:
            raise ValueError(
                "parametric dimension of fitting_spline must be 1!"
            )

        # some values maybe ignored.
        _log.debug(
            "degrees, n_control_points, or knot_vector maybe ignored as "
            "fitting_spline given."
        )

        # create copy, as we will set control points directly
        fitting_spline = fitting_spline.copy()

    else:
        # create bspline based on specifications
        fitting_spline = _prepare_default_bspline1pd(
            *fitting_points.shape, degree, n_control_points, knot_vector, u_k
        )

    # solve system for control points
    residual = solve_for_control_points(
        fitting_points=fitting_points,
        fitting_spline=fitting_spline,
        queries=u_k,
        interpolate_endpoints=interpolate_endpoints,
    )

    if verbose_output:
        return {
            "knot_vector": fitting_spline.knot_vector,
            "degree": fitting_spline.degrees,
            "control_points": fitting_spline.control_points,
            "fitting_points": fitting_points,
            "associated_queries": u_k,
            "residual": residual,
        }

    return fitting_spline, residual


def surface(
    fitting_points,
    size,
    degrees=None,
    n_control_points=None,
    knot_vectors=None,
    fitting_spline=None,
    associated_queries=None,
    centripetal=True,
    interpolate_endpoints=True,
    verbose_output=False,
):
    """
    Fits a surface spline with given parameters through given fitting_points.
    Spline will be interpolated if
    n_control_points[0]*n_control_points[1] = n_fitting_points and
    approximated if not.

    Parameters
    ----------
    fitting_points: (m + 1 x dim) array
        points to be interpolated/approximated
    degrees: array-like (int int)
        degree of fitted spline in both directions
    n_control_points: array-like (int, int)
        number of control points in each direction
    knot_vectors: list of lists [list list]
        knot vectors of fitted spline in both directions
    fitting_spline: spline
        spline used for interpolation/approximation
        must have parametric dimension of 2
    associated_queries: list of lists [list list]
        values where the splines will be evaluated (both directions!)
        will only be used if knot_vectors or fitting_spline is also given!
    centripetal: bool
        Default is True. if True -> centripetal parameterization will be used
    interpolate_endpoints: bool
        if True -> endpoints are interpolated
    verbose_outputs:
        returns additional information as dict

    Returns
    -------
    fitted_spline: spline
        interpolated/approximated spline
    residual: float
        residual (coefficient_matrix @ control_points - fitting_points)
    """
    fitting_points = _np.asanyarray(fitting_points)
    dim = fitting_points.shape[1]

    # get evaluation locations
    u_k = [None, None]
    if associated_queries is not None and (
        knot_vectors is not None or fitting_spline is not None
    ):
        # check dimensions of queries
        if (
            associated_queries[0].shape[1] != 1
            or associated_queries[1].shape[1] != 1
        ):
            raise ValueError(
                "Associated queries in each direction must have dimension 1!"
            )
        u_k[0] = _np.asanyarray(associated_queries[0])
        u_k[1] = _np.asanyarray(associated_queries[1])

    else:
        u_k = parameterize(
            fitting_points=fitting_points, size=size, centripetal=centripetal
        )

    # initialize list for fitting splines
    # structure is needed due to function calls later
    fitting_splines = [None, None]
    if fitting_spline is None:
        # initialize knot_vectors
        if knot_vectors is None:
            knot_vectors = [None, None]
        if degrees is None:
            degrees = [None, None]
        if n_control_points is None:
            n_control_points = [None, None]

        # create base spline for each direction
        spline_u = _prepare_default_bspline1pd(
            size[0],
            dim,
            degrees[0],
            n_control_points[0],
            knot_vectors[0],
            u_k[0],
        )
        spline_v = _prepare_default_bspline1pd(
            size[1],
            dim,
            degrees[1],
            n_control_points[1],
            knot_vectors[1],
            u_k[1],
        )
        fitting_splines = [spline_u, spline_v]
        fitted_spline = type(spline_u)(
            [*spline_u.ds, *spline_v.ds],
            [*spline_u.kvs, *spline_v.kvs],
            _np.empty((len(spline_u.cps) * len(spline_v.cps), dim)),
        )

    else:
        if fitting_spline.para_dim != 2:
            raise ValueError(
                "If fitting spline is given as one spline, "
                "parametric dimension must be 2!"
            )
        # extract spline in each direction (knot_vectors, weights etc.
        # of original spline), works for all type of splines!
        fitting_splines = fitting_spline.extract.boundaries([2, 0])
        fitted_spline = fitting_spline.copy()

    n_control_points = fitted_spline.control_mesh_resolutions

    # index helpers
    mi_pts = MultiIndex(size)
    mi_interim_cps = MultiIndex((n_control_points[0], size[1]))

    # loop first dim
    # curve fit for every j in n_points_v
    interim_control_points = _np.empty((n_control_points[0] * size[1], dim))
    residual_u = _np.empty(size[1])
    fitted_spline_u = fitting_splines[0]
    for v in range(size[1]):
        # previously fitted spline is used for the other fits
        fitted_spline_u, residual_u[v] = curve(
            fitting_points=fitting_points[mi_pts[:, v]],
            fitting_spline=fitted_spline_u,
            associated_queries=u_k[0],
            centripetal=centripetal,
            interpolate_endpoints=interpolate_endpoints,
        )
        # cps in u direction (later fitted in v direction)
        interim_control_points[mi_interim_cps[:, v]] = (
            fitted_spline_u.control_points
        )

    # loop second dim
    # curve fit for every k in n_control_points_u
    residual_v = _np.empty(n_control_points[0])
    fitted_spline_v = fitting_splines[1]
    for u in range(n_control_points[0]):
        # previously fitted spline is used for the other fits
        fitted_spline_v, residual_v[u] = curve(
            fitting_points=interim_control_points[mi_interim_cps[u, :]],
            fitting_spline=fitted_spline_v,
            associated_queries=u_k[1],
            centripetal=centripetal,
            interpolate_endpoints=interpolate_endpoints,
        )

        interim_control_points[mi_interim_cps[u, : n_control_points[1]]] = (
            fitted_spline_v.control_points
        )

    # copy
    fitted_spline.control_points = interim_control_points[
        mi_interim_cps[: n_control_points[0], : n_control_points[1]]
    ]

    residual = _np.linalg.norm(
        (_np.linalg.norm(residual_u), _np.linalg.norm(residual_v))
    )

    if verbose_output:
        return {
            "knot_vectors": fitted_spline.knot_vector,
            "degrees": fitted_spline.degrees,
            "control_points": fitted_spline.control_points,
            "fitting_points": fitting_points,
            "associated_queries": u_k,
            "residual": residual,
        }
    return fitted_spline, residual
