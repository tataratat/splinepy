import numpy as _np
from scipy.sparse.linalg import lsqr as _lsqr

from splinepy import settings as _settings
from splinepy import utils as _utils
from splinepy.utils import log as _log
from splinepy.utils.data import make_matrix as _make_matrix


def parametrize_curve(points, n_points, centripetal):
    """
    Parametrizes the given points to be used in interpolation/approximation.

    Parameters
    ----------
    points: (m + 1 x dim) array
        points to be interpolated/approximated

    n_points: int
        number of points

    centripetal: bool
        if True -> centripetal parametrization will be used

    Returns
    -------
    u_k: (n x 1) array
        array containing parametrization values
    """

    chord_lengths = _np.empty((n_points - 1, 1))
    u_k = _np.empty((n_points, 1))

    total_chord_length = 0.0

    for i in range(1, n_points):
        ith_chord_length = _np.linalg.norm(
            _np.array(points[i, :] - points[(i - 1), :])
        )

        if centripetal:
            ith_chord_length = _np.sqrt(ith_chord_length)

        total_chord_length += ith_chord_length
        chord_lengths[i - 1] = ith_chord_length

    chord_lengths = chord_lengths / total_chord_length
    u_k[0] = 0.0
    u_k[-1] = 1.0

    for i in range(1, n_points - 1):
        u_k[i] = u_k[i - 1] + chord_lengths[i - 1]

    return u_k.reshape(-1, 1)


def compute_knot_vector(degree, n_control_points, u_k, n_points):
    """
    Computes the knot_vector for a spline interpolation/approximation.

    Parameters
    ----------
    degree: int
        degree of spline

    n_control_points: int
        number of control points

    u_k: (n x 1) array
        array containing parametrization values

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
        for i in range(1, n_control_points - degree):
            knot_vector[i + degree] = _np.sum(u_k[i : (i + degree)]) / degree

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
    points, target_spline, queries, interpolate_endpoints
):
    """
    Builds and solves the system to calculate control points.

    Parameters
    ----------
    points: (m + 1 x dim) array
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
        residual of least squares approximation
    """

    # build matrix
    coefficient_matrix = _make_matrix(
        *target_spline.basis_and_support(queries),
        target_spline.control_points.shape[0],
    )

    if (
        interpolate_endpoints
        and target_spline.control_points.shape != points.shape
    ):
        # reduction of matrix necessary for interpolation of endpoints
        target_spline.control_points[0] = points[0]
        # [-1] index doesn't work -> therefore exact index is used
        target_spline.control_points[
            target_spline.control_points.shape[0] - 1
        ] = points[-1]

        if target_spline.control_points.shape[0] == 2:
            # control points = endpoints
            residual = 0
        elif _utils.data.has_scipy:
            # move known values to RHS
            points[1:-1, :] += (
                -coefficient_matrix[1:-1, [0]] * points[0, :]
                - coefficient_matrix[1:-1, [-1]] * points[-1, :]
            ).todense()
            residual = _np.empty(points.shape[1])
            # solve system
            for i in range(0, points.shape[1]):
                (
                    target_spline.control_points[1:-1, i],
                    _,
                    _,
                    residual[i],
                ) = _lsqr(coefficient_matrix[1:-1, 1:-1], points[1:-1, i])[:4]

        else:
            # move known values to RHS
            points[1:-1, :] += (
                -coefficient_matrix[1:-1, [0]] * points[0, :]
                - coefficient_matrix[1:-1, [-1]] * points[-1, :]
            )
            # solve system
            (
                target_spline.control_points[1:-1],
                residual,
            ) = _np.linalg.lstsq(
                coefficient_matrix[1:-1, 1:-1], points[1:-1]
            )[:2]
            residual = _np.sqrt(residual)

    # interpolation or approximation with original system
    elif _utils.data.has_scipy:
        residual = _np.empty(points.shape[1])
        for i in range(0, points.shape[1]):
            target_spline.control_points[:, i], _, _, residual[i] = _lsqr(
                coefficient_matrix, points[:, i]
            )[:4]

    else:
        target_spline.control_points, residual = _np.linalg.lstsq(
            coefficient_matrix, points
        )[:2]
        residual = _np.sqrt(residual)

    residual = _np.linalg.norm(residual)

    return residual


def fit_curve(
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

    degree:

    n_control_points:

    knot_vector:

    target_spline: spline
        spline used for interpolation/approximation

    associated_queries: (n x 1) array
        values where the spline will be evaluated
        will only be used if a knot_vector is also given!

    centripetal:

    interpolate_endpoints: bool
        if True -> endpoints are interpolated

    verbose_outputs:

    Returns
    -------
    fitted_spline: spline
        interpolated/approximated spline

    residual: float
        norm of residuals of least squares approximation
    """

    # calculate n_points due to multiple usage
    n_points = points.shape[0]

    if associated_queries is not None and (
        knot_vector is not None or target_spline is not None
    ):
        # if associated queries are used, knot_vector must be given too!
        u_k = associated_queries

    else:
        u_k = parametrize_curve(
            points=points, n_points=n_points, centripetal=centripetal
        )

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
        if degree and n_control_points and knot_vector is not None:
            _log.warning(
                "Problem is over specified: n_control_points will be "
                "calculated with degree and knot_vector"
            )
            n_control_points = len(knot_vector) - degree - 1

        elif degree is None:
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

        elif n_control_points is None:
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

        elif knot_vector is None:
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
        points=points,
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


class Fitter:
    """
    Helper class to create interpolated or approximated splines
    from points.
    """

    def __init__(self, spl):
        self.spline = spl

    def fit_curve(self, *args, **kwargs):
        return fit_curve(self.spline, *args, **kwargs)

    def parametrize_curve(self, *args, **kwargs):
        return parametrize_curve(self.spline, *args, **kwargs)

    def compute_knot_vector(self, *args, **kwargs):
        return compute_knot_vector(self.spline, *args, **kwargs)

    def solve_for_control_points(self, *args, **kwargs):
        return solve_for_control_points(self.spline, *args, **kwargs)


# Use function docstrings in Extractor functions
Fitter.fit_curve.__doc__ = fit_curve.__doc__
Fitter.parametrize_curve.__doc__ = parametrize_curve.__doc__
Fitter.compute_knot_vector.__doc__ = compute_knot_vector.__doc__
Fitter.solve_for_control_points.__doc__ = solve_for_control_points.__doc__
