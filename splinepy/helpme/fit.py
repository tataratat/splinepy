import numpy as _np

from splinepy import settings as _settings
from splinepy.utils import log as _log
from splinepy.utils.data import has_scipy as _has_scipy
from splinepy.utils.data import make_matrix as _make_matrix

if _has_scipy:
    from scipy.sparse.linalg import spsolve as _spsolve

    # from scipy.sparse.linalg import lsqr as _lsqr


def parametrise(points, size, centripetal):
    """
    parametrises given points to be used in interpolation/approximation.
    Values are then used to build adequate knot vectors.

    Parameters
    ----------
    points: (m + 1 x dim) array
        points to be interpolated/approximated

    size: array-like
        number of points per parametric dimension

    centripetal: bool
        if True -> centripetal parametrisation will be used

    Returns
    -------
    parametric_coordinates: array-like
        array containing parametrisation values per parametric dimension
    """

    def parametrise_line(points, n_points, centripetal):
        u_k = _np.empty((n_points, 1))

        chord_lengths = _np.linalg.norm(_np.diff(points, axis=0), axis=1)

        if centripetal:
            chord_lengths = _np.sqrt(chord_lengths)

        total_chord_length = _np.sum(chord_lengths)

        chord_lengths = chord_lengths / total_chord_length
        u_k[0] = 0.0
        u_k[-1] = 1.0

        u_k[1:-1] = _np.cumsum(chord_lengths[:-1]).reshape(-1, 1)
        return u_k.reshape(-1, 1)

    from splinepy.helpme.multi_index import MultiIndex

    mi = MultiIndex(size)

    # Loop over all dimensions and append each para_coords to list
    n_dimensions = len(size)
    parametric_coordinates = [None] * n_dimensions
    for k in range(n_dimensions):
        u_k = _np.zeros((size[k], 1))
        entry_indices = [slice(None) for _ in range(n_dimensions)]
        if n_dimensions == 1:
            u_k += parametrise_line(
                points=points,
                n_points=size[k],
                centripetal=centripetal,
            )
            parametric_coordinates = [u_k]
        else:
            for i in range(size[k]):
                entry_indices[k] = i
                u_k += parametrise_line(
                    points=points[mi[tuple(entry_indices)]],
                    n_points=size[k],
                    centripetal=centripetal,
                )
            parametric_coordinates[k] = (u_k/size[k])

    return parametric_coordinates


def compute_knot_vector(degree, n_control_points, u_k, n_points):
    """
    Computes the knot_vector for a spline.

    Parameters
    ----------
    degree: int
        degree of spline

    n_control_points: int
        number of control points

    u_k: (n x 1) array-like
        array containing parametrisation values

    n_points: int
        number of points

    Returns
    -------
    knot_vector: ((degree + n_control_points + 1) x 1) array
        array containing knots of spline
    """
    # initialize with m = n + p + 1
    knot_vector = _np.empty(degree + n_control_points + 1)

    knot_vector[: (degree + 1)] = 0.0
    knot_vector[-(degree + 1) :] = 1.0

    # check if interpolation or approximation
    if n_points == n_control_points:
        # interpolation
        # NURBS Book eq. 9.8 (n_points = m + 1, n_control_points = n + 1)
        knot_vector[(degree + 1) : -(degree + 1)] = (
            _np.convolve(u_k.ravel(), _np.ones(degree), "valid")[1:-1] / degree
        )

    else:
        # approximation
        # NURBS Book eq. (9.68) (n_points = m + 1, n_control_points = n + 1)
        d = (n_points) / (n_control_points - degree)
        # NURBS Book eq. (9.69)
        for j in range(1, n_control_points - degree):
            i = int(j * d)
            alpha = (j * d) - i
            knot_vector[j + degree] = (1 - alpha) * u_k[i - 1] + alpha * u_k[i]
    return knot_vector


def solve_for_control_points(
    target_points, target_spline, queries, interpolate_endpoints
):
    """
    Builds and solves the system to calculate control points.

    Parameters
    ----------
    target_points: (m + 1 x dim) array
        points to be interpolated/approximated

    target_spline: spline
        spline used for interpolation/approximation

    queries: (n x 1) array
        values where the spline will be evaluated

    interpolate_endpoints: bool
        if True -> endpoints are interpolated

    Returns
    -------
    residual: float
        residual (coefficient_matrix @ control_points - points)
    """
    # build matrix
    coefficient_matrix = _make_matrix(
        *target_spline.basis_and_support(queries),
        target_spline.control_points.shape[0],
    )
    points = target_points.copy()

    if (
        interpolate_endpoints
        and target_spline.control_points.shape != points.shape
    ):
        # reduction of matrix necessary for interpolation of endpoints
        target_spline.control_points[0] = points[0]
        # [-1] index doesn't work -> therefore actual index is used
        target_spline.control_points[
            target_spline.control_points.shape[0] - 1
        ] = points[-1]

        if target_spline.control_points.shape[0] == 2:
            # control points = endpoints
            residual = _np.linalg.norm(
                coefficient_matrix @ target_spline.control_points - points
            )
        elif _has_scipy:
            # move known values to RHS
            points[1:-1, :] += (
                -coefficient_matrix[1:-1, [0]] * points[0, :]
                - coefficient_matrix[1:-1, [-1]] * points[-1, :]
            ).todense()

            # solve system A^TAx=A^Tb
            # (tried using scipy.sparse.linalg.lsqr but had worse performance)
            at = coefficient_matrix[1:-1, 1:-1].transpose(copy=True)
            target_spline.control_points[1:-1] = _spsolve(
                at @ coefficient_matrix[1:-1, 1:-1], at @ points[1:-1]
            )

            residual = _np.linalg.norm(
                coefficient_matrix @ target_spline.control_points - points
            )
        else:
            # move known values to RHS
            points[1:-1, :] += (
                -coefficient_matrix[1:-1, [0]] * points[0, :]
                - coefficient_matrix[1:-1, [-1]] * points[-1, :]
            )

            # solve system using lstsq

            (
                target_spline.control_points[1:-1],
                residual,
            ) = _np.linalg.lstsq(
                coefficient_matrix[1:-1, 1:-1], points[1:-1]
            )[:2]
            residual = _np.sqrt(_np.linalg.norm(residual))

    # interpolation or approximation with original system
    elif _has_scipy:
        at = coefficient_matrix.transpose(copy=True)
        target_spline.control_points = _spsolve(
            at @ coefficient_matrix, at @ points
        )

        residual = _np.linalg.norm(
            coefficient_matrix @ target_spline.control_points - points
        )

    else:
        target_spline.control_points, residual = _np.linalg.lstsq(
            coefficient_matrix, points
        )[:2]
        residual = _np.sqrt(_np.linalg.norm(residual))

    return residual


def curve(
    points,
    degree=None,
    n_control_points=None,
    knot_vector=None,
    target_spline=None,
    associated_queries=None,
    centripetal=True,
    interpolate_endpoints=True,
    verbose_output=False,
):
    """
    Fits a spline with given parameters through given points.
    Spline will be interpolated if n_control_points = n_points and
    approximated if not.

    Parameters
    ----------
    points: (m + 1 x dim) array
        points to be interpolated/approximated

    degree: int
        degree of spline

    n_control_points: int
        number of control points

    knot_vector: list
        desired knot vector of spline

    target_spline: spline
        spline used for interpolation/approximation

    associated_queries: (n x 1) array
        values where the spline will be evaluated
        will only be used if a knot_vector or target_spline is also given!

    centripetal: bool (default = True)
        if True -> centripetal parametrisation will be used

    interpolate_endpoints: bool
        if True -> endpoints are interpolated

    verbose_outputs:
        returns additional information as dict

    Returns
    -------
    fitted_spline: spline
        interpolated/approximated spline

    residual: float
        residual (coefficient_matrix @ control_points - points)
    """
    # calculate n_points due to multiple usage
    n_points = points.shape[0]

    if associated_queries is not None:
        u_k = associated_queries

    else:
        u_k = parametrise(
            points=points, size=[n_points], centripetal=centripetal
        )[0]

    # check dimension of associated queries
    if associated_queries is not None and associated_queries.shape[1] != 1:
        raise ValueError(
            "Parametric dimension of target_spline "
            "and associated_queries do not match!"
        )

    if target_spline is not None:
        # Sanity checks for target_spline mit raises

        # check compatibility of dimensions
        if points.shape[1] != target_spline.dim:
            raise ValueError(
                "Physical dimension of target_spline "
                "and points do not match!"
            )

        if degree is not None:
            _log.warning("Overwrite degrees (target_spline given).")
        degree = target_spline.degrees[0]

        if n_control_points is not None:
            _log.warning("Ignore n_control_points (target_spline given).")
        n_control_points = target_spline.control_points.shape[0]

        if knot_vector is not None:
            _log.warning("Overwrite knot_vector (target_spline given).")
        knot_vector = target_spline.knot_vectors[0]

    else:
        # sanity checks for given values
        if (
            degree is not None
            and n_control_points is not None
            and knot_vector is not None
        ):
            _log.warning(
                "Problem is over specified: n_control_points will be "
                "calculated with degree and knot_vector"
            )
            n_control_points = len(knot_vector) - degree - 1

        if degree is None:
            if knot_vector and n_control_points is not None:
                _log.info(
                    "Neither degree nor target_vector was given. Degree was "
                    "calculated with knot_vector and n_control_points."
                )
                degree = len(knot_vector) - n_control_points - 1

            else:
                raise ValueError(
                    "Neither degree nor target_vector was given "
                    "and n_control_points or knot_vector is None "
                    "-> unable to calculate degree."
                )

        if n_control_points is None:
            if knot_vector is not None:
                _log.info(
                    "n_control_points was not given and therefore calculated "
                    "with knot_vector and degree"
                )
                n_control_points = len(knot_vector) - degree - 1
            else:
                _log.warning(
                    "Neither n_control_points nor target_vector was given "
                    "and degree or knot_vector is None "
                    "-> set n_control_points = n_points ."
                )
                n_control_points = n_points

        if knot_vector is None:
            if degree >= n_control_points:
                raise ValueError(
                    "Given degree must be lower than n_control_points!"
                )

            knot_vector = compute_knot_vector(
                degree=degree,
                n_control_points=n_control_points,
                u_k=u_k,
                n_points=n_points,
            )

        # build target_spline
        control_points = _np.empty((n_control_points, points.shape[1]))

        target_spline = _settings.NAME_TO_TYPE["BSpline"](
            degrees=[degree],
            knot_vectors=[knot_vector],
            control_points=control_points,
        )

    # check number of points and control points
    if n_control_points > n_points:
        raise ValueError(
            "Requirement (n_points >= n_control_points) not satisfied!"
        )

    # solve system for control points
    residual = solve_for_control_points(
        target_points=points,
        target_spline=target_spline,
        queries=u_k,
        interpolate_endpoints=interpolate_endpoints,
    )

    if verbose_output:
        return {
            "knot_vector": target_spline.knot_vector,
            "degree": target_spline.degrees[0],
            "control_points": target_spline.control_points,
            "points": points,
            "associated_queries": u_k,
            "residual": residual,
        }
    return target_spline, residual


def surface(
    points,
    size,
    degrees=None,
    n_control_points=None,
    knot_vectors=None,
    target_splines=None,
    associated_queries=None,
    centripetal=True,
    interpolate_endpoints=True,
    verbose_output=False,
):
    """
    Fits a surface spline with given parameters through given points.
    Spline will be interpolated if n_control_points = n_points and
    approximated if not.

    Parameters
    ----------
    points: (m + 1 x dim) array
        points to be interpolated/approximated

    degrees: [int int]
        degree of spline in both directions

    n_control_points: int
        number of control points

    knot_vectors: [list list]
        desired knot vectors of spline in both directions

    target_spline: spline or [spline spline]
        spline used for interpolation/approximation
        can be given as one surface spline or two 1D splines

    associated_queries: [list list] or [array array]
        values where the splines will be evaluated (both directions!)
        will only be used if knot_vectors or target_splines is also given!

    centripetal: bool (default = True)
        if True -> centripetal parametrisation will be used

    interpolate_endpoints: bool
        if True -> endpoints are interpolated

    verbose_outputs:
        returns additional information as dict

    Returns
    -------
    fitted_spline: spline
        interpolated/approximated spline

    residual: float
        residual (coefficient_matrix @ control_points - points)
    """
    if target_splines is None:
        target_splines = [None, None]

    elif len(target_splines) == 2:
        if target_splines[0].para_dim != 1 or target_splines[1].para_dim != 1:
            raise ValueError(
                "Parametric dimension of each target spline must be 1!"
            )
        if target_splines[0].whatami != target_splines[1].whatami:
            raise TypeError(
                "Given splines do not match regarding type, dim or para_dim"
            )

    elif len(target_splines) == 1:
        # do something if one surface spline is given instead of 2 single splines
        if target_splines.para_dim != 2:
            raise ValueError(
                "If target spline is given as one spline, "
                "parametric dimension must be 2!"
            )
        # copy splines to avoid changing original spline
        original_target_splines = target_splines.copy()

        target_splines = [None, None]
        # extract spline in each direction (knot_vectors, weights etc.
        # of original spline), works for all type of splines!
        # extraction point is mid-point of parameter space
        target_splines[0] = original_target_splines.extract.spline(
            0,
            (
                original_target_splines.parametric_bounds[0, 0]
                + original_target_splines.parametric_bounds[1, 0] / 2
            ),
        )
        target_splines[1] = original_target_splines.extract.spline(
            1,
            (
                original_target_splines.parametric_bounds[0, 1]
                + original_target_splines.parametric_bounds[1, 1] / 2
            ),
        )

    target_points = points
    u_k = [None, None]

    if associated_queries is not None and (
        knot_vectors is not None or target_splines is not None
    ):
        u_k[0] = _np.array(associated_queries[0])
        u_k[1] = _np.array(associated_queries[1])

        # check dimensions of queries
        if u_k[0].shape[1] != 1 or u_k[1].shape[1] != 1:
            raise ValueError(
                "Parametric dimension of target_spline "
                "and associated_queries do not match!"
            )

    else:
        u_k = parametrise(
            points=target_points, size=size, centripetal=centripetal
        )

    if knot_vectors is None:
        knot_vectors = [None, None]

    # curve fit for every i in n_points_u
    calculated_control_points = _np.empty(
        (n_control_points[0] * size[1], target_points.shape[1])
    )
    residual_u = _np.empty(size[0])

    # first spline is fitted with given values
    fitted_spline_u, residual_u[0] = curve(
        points=target_points[0 : size[0]],
        degree=degrees[0],
        n_control_points=n_control_points[0],
        knot_vector=knot_vectors[0],
        target_spline=target_splines[0],
        associated_queries=u_k[0],
        centripetal=centripetal,
        interpolate_endpoints=interpolate_endpoints,
    )

    calculated_control_points[
        0 : n_control_points[0]
    ] = fitted_spline_u.control_points

    for i in range(1, size[1]):
        # previously fitted spline is used for the other fits
        fitted_spline_u, residual_u[i] = curve(
            points=target_points[i * size[0] : (i + 1) * size[0]],
            target_spline=fitted_spline_u,
            associated_queries=u_k[0],
            centripetal=centripetal,
            interpolate_endpoints=interpolate_endpoints,
        )
        # cps in u direction (later fitted in v direction)
        calculated_control_points[
            i * n_control_points[0] : (i + 1) * n_control_points[0]
        ] = fitted_spline_u.control_points

    # curve fit for every k in n_control_points_u
    residual_v = _np.empty(n_control_points[0])

    # first spline is fitted with given values
    fitted_spline_v, residual_v[0] = curve(
        points=calculated_control_points[
            range(0, n_control_points[0] * size[1], n_control_points[0])
        ],
        degree=degrees[1],
        n_control_points=n_control_points[1],
        knot_vector=knot_vectors[1],
        target_spline=target_splines[1],
        associated_queries=u_k[1],
        centripetal=centripetal,
        interpolate_endpoints=interpolate_endpoints,
    )

    calculated_control_points[
        range(
            0,
            n_control_points[0] * n_control_points[1],
            n_control_points[0],
        )
    ] = fitted_spline_v.control_points

    for k in range(1, n_control_points[0]):
        # previously fitted spline is used for the other fits
        fitted_spline_v, residual_v[k] = curve(
            points=calculated_control_points[
                range(k, n_control_points[0] * size[1], n_control_points[0])
            ],
            target_spline=fitted_spline_v,
            associated_queries=u_k[1],
            centripetal=centripetal,
            interpolate_endpoints=interpolate_endpoints,
        )

        calculated_control_points[
            range(
                k,
                n_control_points[0] * n_control_points[1],
                n_control_points[0],
            )
        ] = fitted_spline_v.control_points

    fitted_spline = _settings.NAME_TO_TYPE["BSpline"](
        degrees=degrees,
        knot_vectors=[
            fitted_spline_u.knot_vectors[0],
            fitted_spline_v.knot_vectors[0],
        ],
        control_points=calculated_control_points[: _np.prod(n_control_points)],
    )

    residual = _np.linalg.norm(
        (_np.linalg.norm(residual_u), _np.linalg.norm(residual_v))
    )

    if verbose_output:
        return {
            "knot_vectors": fitted_spline.knot_vector,
            "degrees": fitted_spline.degrees,
            "control_points": fitted_spline.control_points,
            "points": points,
            "associated_queries": u_k,
            "residual": residual,
        }
    return fitted_spline, residual
