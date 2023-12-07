import numpy as _np
from scipy.sparse.linalg import spsolve as _spsolve
from splinepy import settings as _settings
from splinepy.utils import log as _log
from splinepy.utils.data import make_matrix as _make_matrix
from splinepy.utils.data import has_scipy as _has_scipy

# add utils functionalities


def parametrize_curve(points, centripetal, n_points): 
    """
    Parametrizes the given points to be used in interpolation/approximation.

    Parameters
    ----------

    Returns
    -------
    """
    
    chord_lengths = _np.empty((n_points, 1))
    chord_lengths[0] = 0.0
    chord_lengths[-1] = 1.0

    total_chord_length = 0.0

    for i in range(1, n_points - 1):
        chord_lengths[i] = (
            chord_lengths[i - 1] 
            + _np.linalg.norm(_np.array(points[i, :] 
            - points[(i - 1), :])))
        
        if centripetal:
            chord_lengths[i] = _np.sqrt(chord_lengths[i])
        
        total_chord_length += chord_lengths[i]
    
    u_k = chord_lengths / total_chord_length

    return u_k.reshape(-1,1)

def parametrize_surface(points, size_u, size_v, centripetal):

    # Compute u_k
    temp_u_k = []
    # v - direction
    for v in range(size_v):
        points_u = points[(size_u * v) : (size_u * (v + 1)),:] 
        temp_temp_u_k = parametrize_curve(points_u, size_u, centripetal) 
        temp_u_k = _np.vstack((temp_u_k,temp_temp_u_k))
    
    # Average u - direction
    u_k = _np.sum(temp_u_k,axis=0) / size_v

    # Compute v_l
    temp_v_l = []
    # u - direction
    for u in range(size_u):
        points_v = points[(size_v * u) : (size_v * (u + 1)),:] 
        temp_temp_v_l = parametrize_curve(points_v, size_v, centripetal) 
        temp_v_l = _np.vstack((temp_v_l,temp_temp_v_l))
    
    # Average v - direction
    v_l = _np.sum(temp_v_l,axis=0) / size_u

    return u_k, v_l


def compute_knot_vector(degree, n_control_points, u_k, n_points):

    knot_vector = _np.empty(degree + n_control_points + 1)

    knot_vector[:(degree + 1)] = 0.0
    knot_vector[-(degree + 1):] = 1.0 # max(u_k), ist aber immer = 1?

    # interpolation or approximation
    if n_points == n_control_points:
        # interpolation
        # NURBS Book eq. 9.8
        for i in range(1, n_control_points - degree):
            knot_vector[i + degree] = _np.sum(u_k[i:(i + degree -1)]) / degree

    else:
        # approximation
        # NURBS Book eq. (9.68) (n_control_points = n + 1)
        d = (n_points) / (n_control_points - degree)
        # NURBS Book eq. (9.69)
        for j in range(1, n_control_points - degree):
            i = int(j * d)
            alpha = (j * d) - i
            knot_vector[j + degree] = (
                (1 - alpha) * u_k[i - 1] + alpha * u_k[i]
                )

    return knot_vector

# fitting 

def solve_for_control_points(points, target_spline, queries):

    # choose available solver
    if _has_scipy:
        solving_function = _spsolve
    else:
        solving_function = _np.linalg.solve

    coefficient_matrix = _make_matrix(
        *target_spline.basis_and_support(queries), # ?????? 
         target_spline.control_points.shape[0]
        )

    target_spline.control_points[:, 0] = solving_function(
        coefficient_matrix, points.T
    )

def fit_curve(points, 
              degree = None, 
              n_control_points = None,
              knot_vector = None,
              target_spline = None,
              associated_queries = None,
              centripetal = True,
              interpolate_endpoints = True, 
              verbose_output = False):

    n_points = points.shape[0]
    u_k = parametrize_curve(points, centripetal, n_points)

    # check compatibility of dimensions
    if points.shape[1] != target_spline.dim:
        raise ValueError("Physical dimension of target_spline"
                         "and points do not match!")
    
    if associated_queries.shape[1] != target_spline.para_dim:
        raise ValueError("Parametric dimension of target_spline"
                         "and associated_queries do not match!")  
    
    if target_spline is not None:
        # Sanity checks for target_spline mit raises
        if degree is not None:
            _log._warning("Overwrite degrees (target_spline given).")
        degree = target_spline.degrees

        if n_control_points is not None:
            _log._warning("Ignore n_control_points (target_spline given).")
        n_control_points = target_spline.control_points.shape[0]

        if knot_vector is not None:
            _log._warning("Overwrite knot_vector (target_spline given).")
        knot_vector = target_spline.knot_vector

    else:
        # Sanity checks for given values
        if degree and n_control_points and knot_vector is not None:
            _log.warning("Problem is over specified: n_cps will be" 
                         "calculated with degree and knot_vector")
            n_control_points = len(knot_vector) - degree - 1

        if degree is None:
            if knot_vector and n_control_points is not None:
                degree = len(knot_vector) - n_control_points - 1
                _log.info(
                    "Neither degree nor target_vector was given."
                    "Degree was calculated with knot_vector and n_cps.")
                
            else:
                raise ValueError(
                    "Neither degree nor target_vector was given and n_cps or "
                    "knot_vector is None -> degree couldn't be calculated.")

        if n_control_points is None:
            if knot_vector and degree is None:
                _log.info(
                    "n_control_points was not given and therefore calculated"
                    "with knot_vector and degree")
                n_control_points = len(knot_vector) - degree - 1
            else:
                _log.info(
                    "Neither n_cps nor target_vector was given and degree or "
                    "knot_vector is None -> n_cps is set to n_points .")
                n_control_points = n_points

        if knot_vector is None:
            knot_vector = compute_knot_vector(
                degree, n_control_points, u_k, n_points
                )

        # Check knot_vector dimensions and update if required
        control_points = _np.empty((n_control_points, points.shape[1]))

        # initialize spline 
        target_spline = _settings.NAME_TO_TYPE["BSpline"](
            degrees = [degree],
            knot_vectors = [knot_vector],
            control_points = control_points)
        
    if associated_queries is None:
        # use u_k as queries
        associated_queries = _np.reshape(u_k, (len(u_k),1))
        
    if (n_control_points == n_points):    
        # calculate control_points
        solve_for_control_points(points, target_spline, u_k, associated_queries)
    else:
        approximate_curve(points, target_spline, degree, n_control_points)

def approximate_curve(points, target_spline, degree, n_control_points):
    # Number of variable num_control_points (non-predetermined by P0 and Pm)
    n_variable_control_points = n_control_points - 2
    # Calculate R_k eq. (9.63)
    n_points = points.shape[0]


def fit_surface(points, 
                degree_u = 3, 
                degree_v = 3, 
                n_control_points_u = None,
                n_control_points_v = None,
                knot_vectors = None,
                target_spline = None,
                associated_queries = None,
                centripetal = True,
                interpolate_endpoints = True, 
                verbose_output = False):
    # Parametrize Surface
    u_k, v_l = parametrize_surface(points, n_control_points_u, n_control_points_v, centripetal)
    n_points = points.shape[0]

    knot_vector_u = compute_knot_vector(degree_u, n_points,n_control_points_u,u_k)
    knot_vector_v = compute_knot_vector(degree_v, n_points,n_control_points_v,v_l)

    dimension = points.shape[1]

    control_points_u = _np.empty((n_control_points_u, dimension))
    control_points_v = _np.empty((n_control_points_v, dimension))

    temp_control_points = []
    # u-direction global interpolation
    for v in range(0,n_control_points_v):
        points_u = points[(n_control_points_u * v) : (n_control_points_u * (v + 1)),:]
        spline = _settings.NAME_TO_TYPE["Bezier"](degrees=degree_u, knot_vectors=knot_vector_u, control_points=control_points_u)
        solve_for_control_points(spline,points_u)
        temp_control_points = _np.vstack(temp_control_points,spline.control_points)
    # v-direction global interpolation
    for u in range(0,n_control_points_u):
        points_v = points[(n_control_points_v * u) : (n_control_points_v * (u + 1)),:]
        spline = _settings.NAME_TO_TYPE["Bezier"](degrees=degree_v, knot_vectors=knot_vector_v, control_points=control_points_v)
        solve_for_control_points(spline,points_v)
        temp_control_points = _np.vstack(temp_control_points,spline.control_points)

    return None

# outer layer

def interpolate_curve (points, degree, centripetal, knot_vector):
    pass
def approximate_curve (points, degree, n_control_points,centripetal, knot_vector):
    pass
def interpolate_surface (points, size_u, size_v, degree_u, degree_v, centripetal):
    pass
def approximate_surface (points, n_points_u, n_points_v, size_u, size_v, degree_u, degree_v, centripetal):
    pass