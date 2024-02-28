from functools import wraps as _wraps

import numpy as _np

from splinepy import settings as _settings


def valid_queries(spline, queries):
    """Check queries

    Check if queries are valid for a specific spline, so not to exceed
    parametric bounds for requests

    Parameters
    ----------
    spline : Spline-Type
      spline that sets the bounds
    queries : np.array
      Requested queries

    Returns
    -------
    valid : true
    """
    bounds = spline.parametric_bounds
    # Check dimensions
    if bounds.shape[1] != queries.shape[1]:
        raise ValueError(
            f"Dimension mismatch between parametric dimension of "
            f"spline ({spline.para_dim}), and query-request "
            f"({queries.shape[1]})."
        )

    # Check minimum value
    min_query = _np.min(queries, axis=0)
    if _np.any(bounds[0, :] > min_query):
        error_dim = _np.where(bounds[1, :] > min_query)[0][0]
        error_query = _np.argmin(queries, axis=0)[error_dim]
        raise ValueError(
            f"Query request out of bounds in parametric dimension "
            f"{error_dim}. Detected query {queries[error_query,:]} at "
            f"positions {error_query}, which is out of bounds with "
            f"minimum values {bounds[1,:]}."
        )

    # Check maximum value
    max_query = _np.max(queries, axis=0)
    if _np.any(bounds[1, :] < max_query):
        error_dim = _np.where(bounds[1, :] < max_query)[0][0]
        error_query = _np.argmax(queries, axis=0)[error_dim]
        raise ValueError(
            f"Query request out of bounds in parametric dimension "
            f"{error_dim}. Detected query {queries[error_query,:]} at "
            f"positions {error_query}, which is out of bounds with "
            f"maximum values {bounds[1,:]}."
        )
    return True


def clamped_knot_vectors(spline, warning=True):
    """
    Checks if knot vector is clamped. This also referred as open knot vector.
    If spline doesn't have enough information, this will return None.

    Parameters
    ---------
    spline: Spline

    Returns
    -------
    is_clamped: bool
    """
    if not spline.has_knot_vectors:
        return True

    degrees = spline.degrees
    knot_vectors = spline.knot_vectors
    if degrees is None or knot_vectors is None:
        return None

    for d, kv in zip(degrees, knot_vectors):
        kv_arr = _np.asanyarray(kv)

        front = all(abs(kv_arr[: (d + 1)] - kv_arr[0]) < _settings.TOLERANCE)
        end = all(abs(kv_arr[-(d + 1) :] - kv_arr[-1]) < _settings.TOLERANCE)

        if not front or not end:
            if warning:
                spline._logw(
                    "Spline has an unclamped/closed knot vector. "
                    "This is an experimental feature and may cause "
                    "unexpected behavior."
                )
            return False

    return True


class Checker:
    """Helper class to allow direct extraction from spline obj (BSpline or
    NURBS). Internal use only.

    Examples
    ---------

    .. code-block :: python

      spline.check.valid_queries(queries)

    Parameters
    ----------
    spline : Spline
      Parent spline, that is to be checked
    """

    __slots__ = ("_helpee",)

    def __init__(self, spl):
        self._helpee = spl

    @_wraps(valid_queries)
    def valid_queries(self, queries):
        return valid_queries(self._helpee, queries)

    @_wraps(clamped_knot_vectors)
    def clamped_knot_vectors(self, warning=True):
        return clamped_knot_vectors(self._helpee, warning)
