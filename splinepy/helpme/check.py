import numpy as np


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
    min_query = np.min(queries, axis=0)
    if np.any(bounds[0, :] > min_query):
        error_dim = np.where(bounds[1, :] > min_query)[0][0]
        error_query = np.argmin(queries, axis=0)[error_dim]
        raise ValueError(
            f"Query request out of bounds in parametric dimension "
            f"{error_dim}. Detected query {queries[error_query,:]} at "
            f"positions {error_query}, which is out of bounds with "
            f"minimum values {bounds[1,:]}."
        )

    # Check maximum value
    max_query = np.max(queries, axis=0)
    if np.any(bounds[1, :] < max_query):
        error_dim = np.where(bounds[1, :] < max_query)[0][0]
        error_query = np.argmax(queries, axis=0)[error_dim]
        raise ValueError(
            f"Query request out of bounds in parametric dimension "
            f"{error_dim}. Detected query {queries[error_query,:]} at "
            f"positions {error_query}, which is out of bounds with "
            f"maximum values {bounds[1,:]}."
        )
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

    def __init__(self, spl):
        self._spline = spl

    def valid_queries(self, queries):
        return valid_queries(self._spline, queries)
