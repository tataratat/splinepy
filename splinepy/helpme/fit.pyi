from _typeshed import Incomplete

from splinepy.helpme.multi_index import MultiIndex as MultiIndex

def parameterize(fitting_points, size, centripetal): ...
def compute_knot_vector(degree, n_control_points, u_k, n_fitting_points): ...
def solve_for_control_points(
    fitting_points, fitting_spline, queries, interpolate_endpoints
): ...
def curve(
    fitting_points,
    degree: Incomplete | None = None,
    n_control_points: Incomplete | None = None,
    knot_vector: Incomplete | None = None,
    fitting_spline: Incomplete | None = None,
    associated_queries: Incomplete | None = None,
    centripetal: bool = True,
    interpolate_endpoints: bool = True,
    verbose_output: bool = False,
): ...
def surface(
    fitting_points,
    size,
    degrees: Incomplete | None = None,
    n_control_points: Incomplete | None = None,
    knot_vectors: Incomplete | None = None,
    fitting_spline: Incomplete | None = None,
    associated_queries: Incomplete | None = None,
    centripetal: bool = True,
    interpolate_endpoints: bool = True,
    verbose_output: bool = False,
): ...
