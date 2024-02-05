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
    fitting_points:(m + 1 x dim) array
        points to be interpolated/approximated

    size:array-like
        number of points per parametric dimension

    centripetal:bool
        if True -> centripetal parameterization will be used

    Returns
    -------
    parametric_coordinates:array-like
        array containing parameterization values per parametric dimension
    """

    def parameterize_line(fitting_points, n_fitting_points, centripetal):
        u_k = _np.empty((n_fitting_points, 1))

        chord_lengths = _np.linalg.norm(
            _np.diff(fitting_points, axis=0), axis=1
        )

        if centripetal:
            chord_lengths = _np.sqrt(chord_lengths)

        total_chord_length = _np.sum(chord_lengths)

        chord_lengths = chord_lengths / total_chord_length
        u_k[0] = 0.0
        u_k[-1] = 1.0

        u_k[1:-1] = _np.cumsum(chord_lengths[:-1]).reshape(-1, 1)
        return u_k.reshape(-1, 1)

    mi = MultiIndex(size)

    # Loop over all dimensions and append each para_coords to list
    n_dimensions = len(size)
    if n_dimensions == 1:
        u_k = parameterize_line(
            fitting_points=fitting_points,
            n_fitting_points=size[0],
            centripetal=centripetal,
        )
        parametric_coordinates = [u_k]
    else:
        parametric_coordinates = []
        for k in range(n_dimensions):
            u_k = _np.zeros((size[k], 1))
            entry_indices = [slice(None) for _ in range(n_dimensions)]
            for i in range(size[1 - k]):
                entry_indices[1 - k] = i
                u_k += parameterize_line(
                    fitting_points=fitting_points[mi[tuple(entry_indices)]],
                    n_fitting_points=size[k],
                    centripetal=centripetal,
                )
            parametric_coordinates.append(u_k / size[1 - k])

    return parametric_coordinates


def compute_knot_vector(degree, n_control_points, u_k, n_fitting_points):
    """
    Computes the knot_vector for a spline.

    Parameters
    ----------
    degree:int
        degree of spline

    n_control_points:int
        number of control points

    u_k:(n x 1) array-like
        array containing parameterization values

    n_fitting_points:int
        number of fitting_points

    Returns
    -------
    knot_vector:((degree + n_control_points + 1) x 1) array
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
        # NURBS Book eq. (9.69)
        for j in range(1, n_control_points - degree):
            i = int(j * d)
            alpha = (j * d) - i
            knot_vector[j + degree] = (1 - alpha) * u_k[i - 1] + alpha * u_k[i]
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
        # [-1] index doesn't work -> therefore actual index is used
        fitting_spline.control_points[
            fitting_spline.control_points.shape[0] - 1
        ] = rhs[-1]

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
    fitting_points:(m + 1 x dim) array
        points to be interpolated/approximated

    degree:int
        degree of spline

    n_control_points:int
        number of control points

    knot_vector:list
        desired knot vector of spline

    fitting_spline:spline
        spline used for interpolation/approximation

    associated_queries:(n x 1) array
        values where the spline will be evaluated

    centripetal:bool (default = True)
        if True -> centripetal parameterization will be used

    interpolate_endpoints:bool
        if True -> endpoints are interpolated

    verbose_outputs:
        returns additional information as dict

    Returns
    -------
    fitted_spline:spline
        interpolated/approximated spline

    residual:float
        residual (coefficient_matrix @ control_points - fitting_points)
    """
    # calculate n_fitting_points due to multiple usage
    n_fitting_points = fitting_points.shape[0]

    if associated_queries is not None:
        u_k = associated_queries

    else:
        u_k = parameterize(
            fitting_points=fitting_points,
            size=[n_fitting_points],
            centripetal=centripetal,
        )[0]

    # check dimension of associated queries
    if associated_queries is not None and associated_queries.shape[1] != 1:
        raise ValueError("Dimension of associated_queries must be 1")

    if fitting_spline is not None:
        # Sanity checks for fitting_spline mit raises
        # check compatibility of dimensions
        if fitting_points.shape[1] != fitting_spline.dim:
            raise ValueError(
                "Physical dimension of fitting_spline "
                "and fitting_points do not match!"
            )

        if fitting_spline.para_dim != 1:
            raise ValueError(
                "parametric dimension of fitting_spline "
                "for curve-fit must be 1!"
            )

        if degree is not None:
            _log.warning("Ignore/Overwrite degrees (fitting_spline given).")
        degree = fitting_spline.degrees[0]

        if n_control_points is not None:
            _log.warning(
                "Ignore/Overwrite n_control_points (fitting_spline given)."
            )
        n_control_points = fitting_spline.control_points.shape[0]

        if knot_vector is not None:
            _log.warning(
                "Ignore/Overwrite knot_vector (fitting_spline given)."
            )
        knot_vector = fitting_spline.knot_vectors[0]

    else:
        # sanity checks for given values
        if (
            degree is not None
            and n_control_points is not None
            and knot_vector is not None
        ):
            _log.warning(
                "Problem is over specified:n_control_points will be "
                "calculated with degree and knot_vector"
            )
            n_control_points = len(knot_vector) - degree - 1

        elif degree is None:
            if knot_vector is not None and n_control_points is not None:
                _log.info(
                    "Neither degree nor fitting_vector was given. Degree was "
                    "calculated with knot_vector and n_control_points."
                )
                degree = len(knot_vector) - n_control_points - 1

            else:
                raise ValueError(
                    "Neither degree nor fitting_vector was given "
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
                    "Neither n_control_points nor fitting_vector was given "
                    "and degree or knot_vector is None "
                    "-> set n_control_points = n_fitting_points ."
                )
                n_control_points = n_fitting_points

        if knot_vector is None:
            if degree >= n_control_points:
                raise ValueError(
                    "Given degree must be lower than n_control_points!"
                )

            knot_vector = compute_knot_vector(
                degree=degree,
                n_control_points=n_control_points,
                u_k=u_k,
                n_fitting_points=n_fitting_points,
            )

        # build fitting_spline
        control_points = _np.empty((n_control_points, fitting_points.shape[1]))

        fitting_spline = _settings.NAME_TO_TYPE["BSpline"](
            degrees=[degree],
            knot_vectors=[knot_vector],
            control_points=control_points,
        )

    # check number of fitting_points and control points
    if n_control_points > n_fitting_points:
        raise ValueError(
            "(n_fitting_points >=  n_control_points) not satisfied!"
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

    centripetal: bool (default = True)
        if True -> centripetal parameterization will be used

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
    # initialize list for fitting splines
    # structure is needed due to function calls later
    fitting_splines = [None, None]
    if fitting_spline is not None:
        if fitting_spline.para_dim != 2:
            raise ValueError(
                "If fitting spline is given as one spline, "
                "parametric dimension must be 2!"
            )
        # extract spline in each direction (knot_vectors, weights etc.
        # of original spline), works for all type of splines!
        # extraction point is mid-point of parametric bounds
        para_bounds = fitting_spline.parametric_bounds
        fitting_splines[0] = fitting_spline.extract.spline(
            0,
            (para_bounds[0, 0] + para_bounds[1, 0]) / 2,
        )
        fitting_splines[1] = fitting_spline.extract.spline(
            1,
            (para_bounds[0, 1] + para_bounds[1, 1]) / 2,
        )

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
        u_k[0] = _np.array(associated_queries[0])
        u_k[1] = _np.array(associated_queries[1])

    else:
        u_k = parameterize(
            fitting_points=fitting_points, size=size, centripetal=centripetal
        )

    if knot_vectors is None:
        knot_vectors = [None, None]

    # curve fit for every j in n_points_v
    interim_control_points = _np.empty(
        (n_control_points[0] * size[1], fitting_points.shape[1])
    )
    residual_u = _np.empty(size[1])
    mi_pts = MultiIndex(size)
    mi_interim_cps = MultiIndex((n_control_points[0], size[1]))

    # first spline is fitted with given values
    fitted_spline_u, residual_u[0] = curve(
        fitting_points=fitting_points[mi_pts[:, 0]],
        degree=degrees[0],
        n_control_points=n_control_points[0],
        knot_vector=knot_vectors[0],
        fitting_spline=fitting_splines[0],
        associated_queries=u_k[0],
        centripetal=centripetal,
        interpolate_endpoints=interpolate_endpoints,
    )

    interim_control_points[mi_interim_cps[:, 0]] = (
        fitted_spline_u.control_points
    )

    for v in range(1, size[1]):
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

    # curve fit for every k in n_control_points_u
    residual_v = _np.empty(n_control_points[0])

    # first spline is fitted with given values
    fitted_spline_v, residual_v[0] = curve(
        fitting_points=interim_control_points[mi_interim_cps[0, :]],
        degree=degrees[1],
        n_control_points=n_control_points[1],
        knot_vector=knot_vectors[1],
        fitting_spline=fitting_splines[1],
        associated_queries=u_k[1],
        centripetal=centripetal,
        interpolate_endpoints=interpolate_endpoints,
    )

    interim_control_points[mi_interim_cps[0, : n_control_points[1]]] = (
        fitted_spline_v.control_points
    )

    for u in range(1, n_control_points[0]):
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

    fitted_spline = _settings.NAME_TO_TYPE["BSpline"](
        degrees=degrees,
        knot_vectors=[
            fitted_spline_u.knot_vectors[0],
            fitted_spline_v.knot_vectors[0],
        ],
        control_points=interim_control_points[
            mi_interim_cps[: n_control_points[0], : n_control_points[1]]
        ],
    )

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
