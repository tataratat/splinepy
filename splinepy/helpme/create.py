from functools import wraps as _wraps

import numpy as _np
from gustaf.utils import arr as _arr

from splinepy import settings as _settings
from splinepy.utils import log as _log
from splinepy.utils.data import make_matrix as _make_matrix


def embedded(spline, new_dimension):
    """Embeds Spline in given dimension.
    If dim > current_dim, it will add zeros columns to the control points.
    If dim < current_dim, it takes first n-columns of the control points.
    If dim == current_dim, it returns a copy.

    Parameters
    ----------
    spline: Spline
      (`self`-argument if called via extract member of a spline)
    new_dimension : integer
      new physical dimension into which the spline control points will be set

    Returns
    -------
    embedded_spline : Spline
    """
    from splinepy.spline import Spline as _Spline

    # Check input type
    if not isinstance(spline, _Spline):
        raise NotImplementedError("Embedded only works for splines")

    if new_dimension < 1:
        raise ValueError("Cannot reduce dimension below 1")

    if spline.dim > new_dimension:
        spline_dict = spline.todict()
        spline_dict["control_points"] = spline_dict["control_points"][
            :, :new_dimension
        ]
        return type(spline)(**spline_dict)
    elif spline.dim < new_dimension:
        spline_dict = spline.todict()
        spline_dict["control_points"] = _np.hstack(
            (
                spline_dict["control_points"][:, :new_dimension],
                _np.zeros((spline.cps.shape[0], new_dimension - spline.dim)),
            )
        )
        return type(spline)(**spline_dict)
    else:
        # use todict() to avoid unnecessary data copies from `copy()`
        return type(spline)(**spline.todict())


def extruded(spline, extrusion_vector=None):
    """Extrudes Splines.

    Linear Extrusion of a given spline along an extrusion vector. This will
    extend the parametric dimensionality by one. If the dimension of the
    extrusion vector is higher than the spline's dimension, the spline will be
    extended into the new dimension.

    Parameters
    ----------
    spline: Spline
      (`self`-argument if called via extract member of a spline)
    extrusion_vector: np.ndarray

    Returns
    -------
    extruded_spline : Spline
    """
    from splinepy.spline import Spline as _Spline

    # Check input type
    if not isinstance(spline, _Spline):
        raise NotImplementedError("Extrude only works for splines")

    # Check extrusion_vector
    if extrusion_vector is not None:
        # make flat extrusion_vector
        extrusion_vector = _np.asarray(extrusion_vector).ravel()
    else:
        raise ValueError("No extrusion extrusion_vector given")

    # Check extrusion_vector dimension
    # formulate correct cps
    if spline.dim == extrusion_vector.shape[0]:
        cps = spline.control_points
    elif spline.dim < extrusion_vector.shape[0]:
        expansion_dimension = extrusion_vector.shape[0] - spline.dim
        # one smaller dim is allowed
        # warn that we assume new dim is all zero
        _log.debug(
            f"Given extrusion vector is {expansion_dimension} dimension "
            "bigger than spline's dim. Assuming 0.0 entries for "
            "new dimension.",
        )
        cps = _np.hstack(
            (
                spline.control_points,
                _np.zeros((len(spline.control_points), expansion_dimension)),
            )
        )
    else:
        raise ValueError(
            "Dimension Mismatch between extrusion extrusion vector "
            "and spline."
        )

    # Start Extrusion
    spline_dict = {}

    spline_dict["degrees"] = _np.concatenate((spline.degrees, [1]))
    spline_dict["control_points"] = _np.vstack((cps, cps + extrusion_vector))
    if spline.has_knot_vectors:
        spline_dict["knot_vectors"] = spline.knot_vectors + [[0, 0, 1, 1]]
    if spline.is_rational:
        spline_dict["weights"] = _np.concatenate(
            (spline.weights, spline.weights)
        )

    return type(spline)(**spline_dict)


def revolved(
    spline, axis=None, center=None, angle=None, n_knot_spans=None, degree=True
):
    """Revolve spline around an axis and extend its parametric dimension.

    Parameters
    ----------
    spline: Spline
      (`self`-argument if called via extract member of a spline)
    axis : np.ndarray
      Axis of revolution
    center : np.ndarray
      Center of revolution
    angle : float
      angle of the revolution.
    n_knot_spans : int
      number of non-zero knot-elements for result-spline (if applicable)
    degree : bool
      use degrees instead of radiant

    Returns
    -------
    spline : Spline
    """
    from splinepy.spline import Spline as _Spline

    # Check input type
    if not isinstance(spline, _Spline):
        raise NotImplementedError("Revolutions only works for splines")

    # Check axis
    if axis is not None:
        # Transform into numpy array
        axis = _np.asarray(axis).ravel()
        # Check Axis dimension
        if spline.control_points.shape[1] > axis.shape[0]:
            raise ValueError(
                "Dimension Mismatch between extrusion axis and spline."
            )
        elif spline.control_points.shape[1] < axis.shape[0]:
            _log.debug(
                "Control Point dimension is smaller than axis dimension,"
                " filling with zeros"
            )
            expansion_dimension = axis.shape[0] - spline.dim
            cps = _np.hstack(
                (
                    spline.control_points,
                    _np.zeros(
                        (len(spline.control_points), expansion_dimension)
                    ),
                )
            )
        else:
            cps = _np.copy(spline.control_points)

        # Make sure axis is normalized

        axis_norm = _np.linalg.norm(axis)
        if not _np.isclose(axis_norm, 0, atol=_settings.TOLERANCE):
            axis = axis / axis_norm
        else:
            raise ValueError("Axis-norm is too close to zero.")
    else:
        cps = _np.copy(spline.control_points)
        if spline.control_points.shape[1] == 3:
            raise ValueError("No rotation axis given")

    # Set Problem dimension
    problem_dimension = cps.shape[1]

    # Make sure axis is ignored for 2D
    if problem_dimension == 2:
        axis = None

    # Update angle
    if angle is None:
        spline._logd("No angle given for the revolution. Using 360 degrees.")
        angle = 360 if degree else _np.pi * 2

    if degree:
        angle = _np.radians(angle)

    # Init center
    if center is not None:
        center = _np.asarray(center).ravel()
        # Check Axis dimension
        if not (problem_dimension == center.shape[0]):
            raise ValueError(
                "Dimension Mismatch between axis and center of rotation."
            )
        cps -= center

    # The parametric dimension is independent of the revolution but the
    # rotation-matrix is only implemented for 2D and 3D problems
    if cps.shape[1] not in {2, 3}:
        raise NotImplementedError(
            "Sorry," "revolutions only implemented for 2D and 3D splines"
        )

    # Angle must be (0, pi) non including
    # Rotation is always performed in half steps
    PI = _np.pi
    minimum_n_knot_spans = int(
        _np.ceil((abs(angle) + abs(_settings.TOLERANCE)) / PI)
    )
    if n_knot_spans is None or (n_knot_spans < minimum_n_knot_spans):
        n_knot_spans = minimum_n_knot_spans

    if "Bezier" in spline.name and n_knot_spans > 1:
        raise ValueError(
            "Revolutions are only supported for angles up to 180 "
            "degrees for Bezier type splines as they consist of only "
            "one knot span"
        )

    # Determine auxiliary values
    rot_a = angle / (2 * n_knot_spans)
    half_counter_angle = PI / 2 - rot_a
    weight = _np.sin(half_counter_angle)
    factor = 1 / weight

    # Determine rotation matrix
    rotation_matrix = _arr.rotation_matrix_around_axis(
        axis=axis, rotation=rot_a, degree=False
    ).T

    # Start Extrusion
    spline_dict = {}

    spline_dict["degrees"] = _np.concatenate((spline.degrees, [2]))

    spline_dict["control_points"] = cps
    end_points = cps
    for _i_segment in range(n_knot_spans):
        # Rotate around axis
        mid_points = _np.matmul(end_points, rotation_matrix)
        end_points = _np.matmul(mid_points, rotation_matrix)
        # Move away from axis using dot-product tricks
        if problem_dimension == 3:
            mp_scale = axis * _np.dot(mid_points, axis).reshape(-1, 1)
            mid_points = (mid_points - mp_scale) * factor + mp_scale
        else:
            mid_points *= factor
        spline_dict["control_points"] = _np.concatenate(
            (spline_dict["control_points"], mid_points, end_points)
        )

    if spline.has_knot_vectors:
        kv = [0, 0, 0]
        [kv.extend([i + 1, i + 1]) for i in range(n_knot_spans - 1)]
        spline_dict["knot_vectors"] = spline.knot_vectors + [
            kv + [n_knot_spans] * 3
        ]
    if spline.is_rational:
        mid_weights = spline.weights * weight
        spline_dict["weights"] = spline.weights
        for _i_segment in range(n_knot_spans):
            spline_dict["weights"] = _np.concatenate(
                (spline_dict["weights"], mid_weights, spline.weights)
            )
    else:
        _log.debug(
            "True revolutions are only possible for rational spline types.",
            "Creating Approximation.",
        )

    if center is not None:
        spline_dict["control_points"] += center

    return type(spline)(**spline_dict)


def swept(
    cross_section,
    trajectory,
    cross_section_normal=None,
    set_on_trajectory=False,
    rotation_adaption=None,
):
    """
    Sweeps a cross-section along a trajectory. The cross-section
    receives rotation into the direction of the trajectory tangent
    vector and is then placed either at the evaluation points of the
    trajectory's knots or at the trajectory's control points. This
    depends on the value of the set_on_trajectory parameter. It can
    create both a surface or a solid, depending on the dimension of
    the cross-section.

    The sweeping process has some limitations, since the cross-section
    cannot be preserved exactly along the whole trajectory.

    The maths behind can be found in the NURBS Book, Piegl & Tiller,
    2nd edition, chapter 10.4 Swept Surfaces.

    Parameters
    ----------
    cross_section : Spline
      Cross-section to be swept
    trajectory : Spline
      Trajectory along which the cross-section is swept
    cross_section_normal : np.array
      Normal vector of the cross-section
      Default is [0, 0, 1]
    set_on_trajectory : bool
      If True, the cross-section will be placed at the evaluation
      points of the trajectory's knots. If False, the cross-section
      will be placed at the control points of the trajectory.
      Default is False.
    rotation_adaption : float
      Angle in radians by which the cross-section is rotated around
      the trajectory tangent vector. This is an additional rotation
      if the user wants to adapt the cross-section rotation.
      Example with rectangular cross-section:
      .         x
      .       x   x                                      x x x x x
      .     x       x                                    x       x
      .   x           x    --> rotation around pi/4 -->  x       x
      .     x       x                                    x       x
      .       x   x                                      x x x x x
      .         x

    Returns
    -------
    swept_spline : Spline
      Spline resulting from the sweep
    """

    from splinepy import NURBS as _NURBS
    from splinepy import BSpline as _BSpline

    ### INPUT CHECKS ###

    if not isinstance(cross_section, (_BSpline, _NURBS)):
        raise TypeError(
            "cross_section must be an instance of BSpline or NURBS"
        )
    if not isinstance(trajectory, (_BSpline, _NURBS)):
        raise TypeError("trajectory must be an instance of BSpline or NURBS")
    if not trajectory.para_dim == 1:
        raise TypeError("trajectory must be of parametric dimension 1")
    if not isinstance(set_on_trajectory, bool):
        raise TypeError("set_on_trajectory must be a boolean")
    if not isinstance(rotation_adaption, (float, int, type(None))):
        raise TypeError("rotation_adaption must be a float or int")

    if cross_section_normal is not None and not len(cross_section_normal) == 3:
        raise ValueError("cross_section_normal must be a 3D vector")
    # setting default value for cross_section_normal
    if cross_section_normal is None:
        cross_section_normal = _np.array([0, 0, 1])
        # add debug message
        _log.debug("No cross_section_normal given. Defaulting to [0, 0, 1].")

    ### STARTING CALCULATIONS ###

    # make copies so we can work on it inplace
    trajectory = trajectory.create.embedded(3)
    cross_section = cross_section.create.embedded(3)

    # initialize parameter values
    par_value = trajectory.greville_abscissae()
    par_value = par_value.reshape(-1, 1)

    ### TRANSFORMATION MATRICES ###

    # tangent vector 'e1' of trajectory at parameter value 0
    e1 = trajectory.derivative([par_value[0]], [1])
    e1 = (e1 / _np.linalg.norm(e1)).ravel()

    # evaluating a vector normal to e1
    vec = [-e1[1], e1[0], -e1[2]]
    B = []
    # avoid dividing by zero
    if _np.linalg.norm(_np.cross(e1, vec)) > _settings.TOLERANCE:
        B.append(_np.cross(e1, vec) / _np.linalg.norm(_np.cross(e1, vec)))
    else:
        vec = [e1[2], -e1[1], e1[0]]
        B.append(_np.cross(e1, vec) / _np.linalg.norm(_np.cross(e1, vec)))
        # add debug message
        _log.debug(
            "Division by zero occurred. Using alternative vector for B."
        )

    # initialize transformation matrices and tangent-vector-collection
    T = []
    A = []
    tang_collection = []

    # evaluating transformation matrices for each trajectory point
    for i in range(len(par_value)):
        # calculation according to NURBS Book, eq. 10.27
        # tangent vector e1 on trajectory at parameter value i
        e1 = trajectory.derivative([par_value[i]], [1])
        e1 = (e1 / _np.linalg.norm(e1)).ravel()
        # collecting tangent vectors for later use
        tang_collection.append(e1)

        # projecting B_(i) onto the plane normal to e1
        B.append(B[i] - _np.dot(B[i], e1) * e1)
        B[i + 1] = B[i + 1] / _np.linalg.norm(B[i + 1])

        # defining e2 and e3 vectors
        e3 = B[i + 1]
        e2 = _np.cross(e3, e1)

        # array of transformation matrices from global to local coordinates
        T.append(_np.vstack((e1, e2, e3)))

        # array of transformation matrices from local to global coordinates
        A.append(_np.linalg.inv(T[i]))

    # separate procedure, if trajectory is closed and B[0] != B[-1]
    # recalculate B-vector and middle the values between B and B_rec
    # according to NURBS Book, Piegl & Tiller, 2nd edition, p. 483
    is_trajectory_closed = _np.array_equal(
        trajectory.evaluate([[0]]), trajectory.evaluate([par_value[-1]])
    )
    is_B_start_equal_B_end = _np.array_equal(B[0], B[-1])

    if is_trajectory_closed and not is_B_start_equal_B_end:
        # reset transformation matrices
        T = []
        A = []
        # preallocate B_rec
        B_rec = [None] * len(B)
        # make sure start of B_rec is equal to the end of B
        B_rec[0] = B[-1]
        # redo the calculation of B using tang_collection from before
        # in order to avoid recalculating the tangent vectors;
        # calculation according to NURBS Book, eq. 10.27
        for i in range(len(par_value)):
            B_rec[i + 1] = (
                B_rec[i]
                - _np.dot(B_rec[i], tang_collection[i]) * tang_collection[i]
            )
            B_rec[i + 1] = B_rec[i + 1] / _np.linalg.norm(B_rec[i + 1])
            # middle point between B and B_rec
            B_rec[i + 1] = (B[i + 1] + B_rec[i + 1]) * 0.5
            # normalizing B_rec
            B_rec[i + 1] = B_rec[i + 1] / _np.linalg.norm(B_rec[i + 1])
            # defining e2 and e3 axis-vectors
            e3 = B_rec[i + 1]
            e2 = _np.cross(e3, tang_collection[i])

            # array of transformation matrices from global to local coordinates
            T.append(_np.vstack((tang_collection[i], e2, e3)))

            # array of transformation matrices from local to global coordinates
            A.append(_np.linalg.inv(T[i]))

        # check if the beginning and the end of the B-vector are the same
        if not _np.allclose(B_rec[0], B_rec[-1], rtol=1e-3):
            _log.warning(
                "Vector calculation is not exact due to the "
                "trajectory being closed in an uncommon way."
            )

    ### ROTATION MATRIX ###
    # rotates cross-section normal vector into global e1 direction

    # skip calculation if cross-section normal vector is default
    if _np.array_equal(cross_section_normal, _np.array([0, 0, 1])):
        R = _np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]])
    # else, calculate rotation matrix for given cross-section normal vector
    else:
        # calculate angle of cross-section normal vector around e2-axis
        angle_of_cs_normal_1 = _np.arctan2(
            cross_section_normal[2], cross_section_normal[0]
        )

        # calculate angle of cross-section normal vector around e3-axis
        angle_of_cs_normal_2 = _np.arctan2(
            cross_section_normal[1], cross_section_normal[2]
        )

        # calculate rotation matrix for cross-section normal vector
        R1 = _arr.rotation_matrix_around_axis(
            axis=[0, 1, 0], rotation=angle_of_cs_normal_1, degree=False
        )
        R2 = _arr.rotation_matrix_around_axis(
            axis=[0, 0, 1], rotation=angle_of_cs_normal_2, degree=False
        )
        R = _np.matmul(R2, R1)

    # rotate cross-section around trajectory tangent vector (e1) if wanted
    if rotation_adaption is not None:
        R = _np.matmul(
            _arr.rotation_matrix_around_axis(
                axis=[1, 0, 0], rotation=rotation_adaption, degree=False
            ),
            R,
        )

    # remove numerical noise
    R = _np.where(_np.abs(R) < _settings.TOLERANCE, 0, R)

    ### SWEEPING PROCESS ###

    # evaluate center of cross-section and translate to origin
    cross_para_center = _np.mean(cross_section.parametric_bounds, axis=0)
    cs_center = cross_section.evaluate(
        cross_para_center.reshape(-1, cross_section.para_dim)
    ).ravel()
    cross_section.control_points -= cs_center

    # set cross-section control points along trajectory
    swept_spline_cps = []
    for i, par_val in enumerate(par_value):
        temp_csp = []
        # evaluate trajectory if user wants to set cross-section on trajectory.
        if set_on_trajectory:
            evals = trajectory.evaluate([par_val]).ravel()
        # place every control point of cross-section separately
        for cscp in cross_section.control_points:
            # rotate cross-section in trajectory direction
            normal_cscp = _np.matmul(R, cscp)
            # transform cross-section to global coordinates
            normal_cscp = _np.matmul(A[i], normal_cscp)
            # check if user wants to place cross-section at
            # evaluation points or control points of trajectory
            if set_on_trajectory:
                # translate cross-section to trajectory evaluation point
                normal_cscp += evals
            else:
                # translate cross-section to trajectory control point
                normal_cscp += trajectory.control_points[i]
            # append control point to list
            temp_csp.append(normal_cscp)

        # collect all control points
        swept_spline_cps.append(_np.array(temp_csp))

    # create spline dictionary
    dict_swept_spline = {
        "degrees": [*cross_section.degrees, *trajectory.degrees],
        "knot_vectors": [
            *cross_section.knot_vectors,
            *trajectory.knot_vectors,
        ],
        "control_points": _np.asarray(swept_spline_cps).reshape(
            -1, cross_section.dim
        ),
    }

    # add weights properly if spline is rational
    if cross_section.is_rational or trajectory.is_rational:

        def weights(spline):
            if spline.is_rational:
                return spline.weights
            return _np.ones(spline.control_points.shape[0])

        trajectory_weights = weights(trajectory)
        cross_section_weights = weights(cross_section)
        dict_swept_spline["weights"] = _np.outer(
            trajectory_weights, cross_section_weights
        ).reshape(-1, 1)
        spline_type = _NURBS
    else:
        spline_type = _BSpline

    # create swept spline
    swept_spline = spline_type(**dict_swept_spline)

    return swept_spline


def from_bounds(parametric_bounds, physical_bounds):
    """Creates a minimal spline with given parametric bounds, physical bounds.
    Physical bounds can have less or equal number of
    dimension as parametric bounds. (Greater is not supported)

    Parameters
    -----------
    parametric_bounds: (2, n) array-like
    physical_bounds: (2, n) array-like

    Returns
    --------
    spline: BSpline
    """
    physical_bounds = _np.asanyarray(physical_bounds).reshape(2, -1)
    parametric_bounds = _np.asanyarray(parametric_bounds).reshape(2, -1)

    # get correctly sized bez box
    phys_size = physical_bounds[1] - physical_bounds[0]

    # minimal bspline
    bspline_box = box(*phys_size).bspline  # kvs are in [0, 1]
    bspline_box.cps += physical_bounds[0]  # apply offset

    # update parametric bounds
    new_kvs = []
    para_size = parametric_bounds[1] - parametric_bounds[0]
    for i, kv in enumerate(bspline_box.kvs):
        # apply scale and offset
        new_kv = (kv * para_size[i]) + parametric_bounds[0][i]
        new_kvs.append(new_kv)

    # update at once
    bspline_box.kvs = new_kvs

    return bspline_box


def determinant_spline(spline):
    """
    Creates a spline representing the projection of the jacobian determinant

    dmin > 0 --> spline is not tangled
    dmin <= 0 --> spline is COULD BE tangled

    Only works if (parameter dimension == physical dimension).
    A definitive statement about the entanglement can only be made with
    non-rational splines or rational splines with equal weights!
    Otherwise, the resulting spline, and therefore the entanglement check
    is just an approximation.

    references:
    Mantzaflaris, A., Jüttler, B., Khoromskij, B. N., & Langer, U. (2017)
    Limkilde, A., Evgrafov, A., Gravesen, J., & Mantzaflaris, A. (2021)

    Parameters
    ----------
    spline: Spline

    Returns
    -------
    determinant_projection: BSpline
      Spline which represents the jacobian determinant of the
      given spline object.
    """
    from splinepy.utils.data import has_scipy as _has_scipy

    # choose solver (numpy/scipy)
    if _has_scipy:
        from scipy.sparse.linalg import spsolve as _spsolve

        solving_function = _spsolve

    else:
        solving_function = _np.linalg.solve

    # checks if dimensions are equal.
    if spline.dim != spline.para_dim:
        raise ValueError(
            "Dimension of parameter space and physical space"
            "are not the same ({spline.whatami})!"
        )

    # checks if spline is rational and if weights differ
    # otherwise NURBS can be treated like regular BSpline

    if spline.is_rational and not _np.allclose(
        spline.weights, spline.weights.ravel()[0]
    ):
        _log.warning(
            f"Spline is rational ({spline.whatami}) with unequal weights!"
            "Result is only an approximation."
        )

    # calculate degrees of det spline
    degrees_determinant_spline = spline.degrees * spline.dim - 1
    # calculate necessary additional knot multiplicity due to degree elevation
    multiplicity_increase = degrees_determinant_spline - spline.degrees

    # knot_vectors
    if spline.has_knot_vectors:
        # initialize lists for knot vectors
        knot_vectors_determinant_spline = []

        for u_kv, kn_m, m in zip(
            spline.unique_knots,
            spline.knot_multiplicities,
            multiplicity_increase,
        ):
            # increase knot multiplicities:
            #   @ each inner knot -> mult_inc + 1
            #   @ begin & end -> mult_inc
            # since continuity at inner knots further decreases by 1
            kn_m[1:-1] += 1
            temp_knot_vector = _np.repeat(u_kv, kn_m + m)
            knot_vectors_determinant_spline.append(temp_knot_vector)

        # number of cpts
        n_control_points = _np.prod(
            [
                len(kvs_ds) - d_ds - 1
                for kvs_ds, d_ds in zip(
                    knot_vectors_determinant_spline, degrees_determinant_spline
                )
            ]
        )
        determinant_projection = _settings.NAME_TO_TYPE["BSpline"](
            degrees=degrees_determinant_spline,
            knot_vectors=knot_vectors_determinant_spline,
            control_points=_np.empty((n_control_points, 1)),
        )

    else:
        # if Bezier Spline
        n_control_points = _np.prod(degrees_determinant_spline + 1)
        determinant_projection = _settings.NAME_TO_TYPE["Bezier"](
            degrees=degrees_determinant_spline,
            control_points=_np.empty((n_control_points, 1)),
        )

    # get det(J) of spline at greville pts
    # to calculate cpts of determinant Spline
    sample_queries = determinant_projection.greville_abscissae(
        duplicate_tolerance=_settings.TOLERANCE
    )
    jacobian_determinants = _np.linalg.det(spline.jacobian(sample_queries))
    coefficient_matrix = _make_matrix(
        *determinant_projection.basis_and_support(sample_queries),
        n_control_points,
    )
    determinant_projection.control_points[:, 0] = solving_function(
        coefficient_matrix, jacobian_determinants
    )

    return determinant_projection


def parametric_view(spline, axes=True, conform=False):
    """Create parametric view of given spline. Previously called
    `naive_spline()`. Degrees are always 1 and knot multiplicity is not
    preserved. Returns BSpline, as BSpline and NURBS should look the same as
    parametric view.
    Will take shallow copy of underlying data of spline_data and show_options
    from original spline.
    However, if conforming basis is desired, set confrom=True.

    Parameters
    -----------
    spline: BSpline or NURBS
    axes: bool
      If True, will configure axes settings, it is supported.
    conform: bool
      Default is False


    Returns
    --------
    para_spline: BSpline
    """
    p_bounds = spline.parametric_bounds
    para_spline = from_bounds(
        parametric_bounds=p_bounds, physical_bounds=p_bounds
    )

    # process to create conforming para_view splines
    if conform:
        # conforming spline - rational?
        if spline.is_rational:
            para_spline = para_spline.nurbs

        # conforming spline - bezier?
        if not spline.has_knot_vectors:
            dict_spline = para_spline.todict()
            dict_spline.pop("knot_vectors")
            para_spline = type(spline)(**dict_spline)

        # conforming degrees
        d_diff = spline.degrees - para_spline.degrees
        if any(d_diff):
            # turn difference in degrees into degree query
            d_query = []
            for i, dq in enumerate(d_diff):
                d_query += [i] * dq
            para_spline.elevate_degrees(d_query)

        # process knots to insert
        if spline.has_knot_vectors:
            for i, (kv, d) in enumerate(
                zip(spline.knot_vectors, spline.degrees)
            ):
                n_repeating = int(d + 1)
                query = kv[n_repeating:-n_repeating]
                if len(query) > 0:
                    para_spline.insert_knots(i, query)

        # conforming weights
        if spline.is_rational:
            para_spline.ws[:] = spline.ws

    elif spline.has_knot_vectors:
        # in non-conform case, we can just add missing unique knots
        for i, uk in enumerate(spline.unique_knots):
            query = uk[1:-1]
            if len(query) > 0:
                para_spline.insert_knots(i, query)

    # take shallow copy
    para_spline.spline_data._saved = spline.spline_data._saved.copy()
    para_spline.show_options._options = spline.show_options._options.copy()

    if axes and "axes" in spline.show_options.valid_keys():
        # configure axes
        bs = p_bounds
        bs_diff_001 = (bs[1] - bs[0]) * 0.001
        lower_b = bs[0] - bs_diff_001
        upper_b = bs[1] + bs_diff_001
        axes_config = {
            "xtitle": "u",
            "xrange": [lower_b[0], upper_b[0]],
            "tip_size": 0,
            "xminor_ticks": 3,
            "xygrid": False,
        }
        if spline.para_dim == 2:
            axes_config.update(
                ytitle="v",
                yrange=[lower_b[1], upper_b[1]],
                yminor_ticks=3,
                yzgrid=False,
            )

        if spline.para_dim == 3:
            axes_config.update(
                ztitle="w",
                zrange=[lower_b[2], upper_b[2]],
                zminor_ticks=3,
                zxgrid=False,
            )

        para_spline.show_options["axes"] = axes_config
        # it is a view, so cps won't be realistic
        para_spline.show_options["control_points"] = False
        para_spline.show_options["lighting"] = "off"

    return para_spline


def line(points):
    """Create a spline with the provided points as control points.

    Parameters
    ----------
    points: (n, d) numpy.ndarray
      npoints x ndims array of control points

    Returns
    -------
    line: BSpline
      Spline degree [1].
    """
    # lines have degree 1
    degree = 1

    cps = _np.array(points)
    nknots = cps.shape[0] + degree + 1

    knots = _np.concatenate(
        (
            _np.full(degree, 0.0),
            _np.linspace(0.0, 1.0, nknots - 2 * degree),
            _np.full(degree, 1.0),
        )
    )

    spline = _settings.NAME_TO_TYPE["BSpline"](
        control_points=cps, knot_vectors=[knots], degrees=[degree]
    )

    return spline


def arc(
    radius=1.0,
    angle=90.0,
    n_knot_spans=-1,
    start_angle=0.0,
    degree=True,
):
    """Creates a 1-D arc as Rational Bezier or NURBS with given radius and
    angle. The arc lies in the x-y plane and rotates around the z-axis.

    Parameters
    ----------
    radius : float, optional
      radius of the arc, defaults to 1
    angle : float, optional
      angle of the section of the arc, defaults to 90 degrees
    n_knot_spans : int
      Number of knot spans, by default minimum number for angle is used.
    start_angle : float, optional
      starting point of the angle, by default 0.
    degree: bool, optional
      degrees for angle used, by default True

    Returns
    -------
    arc: NURBS or RationalBezier
    """
    # Define point spline of degree 0 at starting point of the arc
    if degree:
        start_angle = _np.radians(start_angle)
        angle = _np.radians(angle)
    start_point = [
        radius * _np.cos(start_angle),
        radius * _np.sin(start_angle),
    ]
    point_spline = _settings.NAME_TO_TYPE["RationalBezier"](
        degrees=[0], control_points=[start_point], weights=[1.0]
    )
    # Bezier splines only support angles lower than 180 degrees
    if abs(angle) >= _np.pi or n_knot_spans > 1:
        point_spline = point_spline.nurbs

    # Revolve - set degree to False, since all the angles are converted to rad
    arc_attrib = point_spline.create.revolved(
        angle=angle, n_knot_spans=n_knot_spans, degree=False
    ).todict()
    # Remove the first parametric dimensions, which is only a point and
    # only used for the revolution
    arc_attrib["degrees"] = list(arc_attrib["degrees"])[1:]
    if point_spline.has_knot_vectors:
        arc_attrib["knot_vectors"] = list(arc_attrib["knot_vectors"])[1:]

    return type(point_spline)(**arc_attrib)


def circle(radius=1.0, n_knot_spans=3):
    """Circle (parametric dim = 1) with radius r in the x-y plane around the
    origin. The spline has an open knot vector and degree 2.

    Parameters
    ----------
    radius : float, optional
      radius, defaults to one
    n_knots_spans : int, optional
      number of knot spans, defaults to 3

    Returns
    -------
    circle: NURBS
    """
    return arc(radius=radius, angle=360, n_knot_spans=n_knot_spans)


def box(*lengths):
    """N-D box (hyper rectangle).

    Parameters
    ----------
    *lengths: list(float)

    Returns
    -------
    nd_box: Bezier
    """
    # may dim check here?
    # starting point
    nd_box = _settings.NAME_TO_TYPE["Bezier"](
        degrees=[1], control_points=[[0], [lengths[0]]]
    )
    # use extrude
    for i, len_ in enumerate(lengths[1:]):
        nd_box = nd_box.create.extruded([0] * int(i + 1) + [len_])

    return nd_box


def plate(radius=1.0):
    """Creates a biquadratic 2-D spline in the shape of a plate with given
    radius.

    Parameters
    ----------
    radius : float, optional
      Radius of the plate, defaults to one

    Returns
    -------
    plate: RationalBezier
    """
    degrees = [2, 2]
    control_points = (
        _np.array(
            [
                [-0.5, -0.5],
                [0.0, -1.0],
                [0.5, -0.5],
                [-1.0, 0.0],
                [0.0, 0.0],
                [1.0, 0.0],
                [-0.5, 0.5],
                [0.0, 1.0],
                [0.5, 0.5],
            ]
        )
        * radius
    )
    weights = _np.tile([1.0, 1 / _np.sqrt(2)], 5)[:-1]

    return _settings.NAME_TO_TYPE["RationalBezier"](
        degrees=degrees, control_points=control_points, weights=weights
    )


def disk(
    outer_radius,
    inner_radius=None,
    angle=360.0,
    n_knot_spans=4,
    degree=True,
):
    """Surface spline describing a potentially hollow disk with quadratic
    degree along curved dimension and linear along thickness. The angle
    describes the returned part of the disk.

    Parameters
    ----------
    outer_radius : float
      Outer radius of the disk
    inner_radius : float, optional
      Inner radius of the disk, in case of hollow disk, by default 0.
    angle : float, optional
      Rotational angle, by default 360. describing a complete revolution
    n_knot_spans : int, optional
      Number of knot spans, by default 4

    Returns
    -------
    disk: NURBS
      Surface NURBS of degrees (1,2)
    """
    if inner_radius is None:
        inner_radius = 0.0

    cps = _np.array([[inner_radius, 0.0], [outer_radius, 0.0]])
    weights = _np.ones([cps.shape[0]])
    knots = _np.repeat([0.0, 1.0], 2)

    return _settings.NAME_TO_TYPE["NURBS"](
        control_points=cps, knot_vectors=[knots], degrees=[1], weights=weights
    ).create.revolved(angle=angle, n_knot_spans=n_knot_spans, degree=degree)


def torus(
    torus_radius,
    section_outer_radius,
    section_inner_radius=None,
    torus_angle=None,
    section_angle=None,
    section_n_knot_spans=4,
    torus_n_knot_spans=4,
    degree=True,
):
    """Creates a volumetric NURBS spline describing a torus revolved around the
    x-axis. Possible cross-sections are plate, disk (yielding a tube) and
    section of a disk.

    Parameters
    ----------
    torus_radius : float
      Radius of the torus
    section_outer_radius : float
      Radius of the section of the torus
    section_inner_radius : float, optional
      Inner radius in case of hollow torus, by default 0.
    torus_angle : float, optional
      Rotational angle of the torus, by default None, giving a complete
      revolution
    section_angle : float, optional
      Rotational angle, by default None, yielding a complete revolution
    section_n_knot_spans : float, optional
      Number of knot spans along the cross-section, by default 4
    torus_n_knot_spans : float, optional
      Number of knot spans along the torus, by default 4

    Returns
    -------
    torus: NURBS
      Volumetric spline in the shape of a torus with degrees (1,2,2)
    """

    if torus_angle is None:
        torus_angle = 2 * _np.pi
        degree = False

    if section_angle is None:
        section_angle = 2 * _np.pi
        section_angle_flag = False
        degree = False
    else:
        section_angle_flag = True

    if section_inner_radius is None:
        section_inner_radius = 0
        section_inner_radius_flag = False
    else:
        section_inner_radius_flag = True

    # Create the cross-section
    if not section_angle_flag and not section_inner_radius_flag:
        cross_section = plate(section_outer_radius)
        # For more than 180 degree only NURBS can be used
        if abs(torus_angle) >= _np.pi:
            cross_section = cross_section.nurbs
    else:
        cross_section = disk(
            outer_radius=section_outer_radius,
            inner_radius=section_inner_radius,
            n_knot_spans=section_n_knot_spans,
            angle=section_angle,
            degree=degree,
        )

    # Create a surface spline representing a disk and move it from the origin
    cross_section.control_points[:, 1] += torus_radius

    return cross_section.create.revolved(
        axis=[1.0, 0, 0],
        center=_np.zeros(3),
        angle=torus_angle,
        n_knot_spans=torus_n_knot_spans,
        degree=degree,
    )


def sphere(
    outer_radius,
    inner_radius=None,
    angle=360.0,
    n_knot_spans=-1,
    degree=True,
):
    """Creates a volumetric spline describing a sphere with radius R.

    Parameters
    ----------
    outer_radius : float
      Outer radius of the sphere
    inner_radius : float, optional
      Inner radius of the potentially hollow sphere.
    angle : float
      Rotational angle around x-axis, by default each 360
      (describing a complete revolution)
    n_knot_spans : int
      Number of knot spans

    Returns
    -------
    sphere: NURBS
      Volumetric NURBS with degrees (1,2,2)
    """
    if inner_radius is not None:
        inner_radius = float(inner_radius)

    sphere = disk(
        outer_radius, inner_radius, angle=180, degree=True
    ).nurbs.create.revolved(
        axis=[1, 0, 0],
        angle=angle,
        n_knot_spans=n_knot_spans,
        degree=degree,
    )

    return sphere


def surface_circle(outer_radius):
    """
    Create a circle consisting of a single surface patch

    Parameters
    ----------
    outer_rarius : float
      radius

    Returns
    -------
    patch : spline
      Circle spline
    """
    aux_0_w = 2**-0.5
    aux_0 = outer_radius * aux_0_w

    return _settings.NAME_TO_TYPE["RationalBezier"](
        degrees=[2, 2],
        control_points=[
            [-aux_0, -aux_0],
            [0, -2 * aux_0],
            [aux_0, -aux_0],
            [-2 * aux_0, 0],
            [0, 0],
            [2 * aux_0, 0],
            [-aux_0, aux_0],
            [0, 2 * aux_0],
            [aux_0, aux_0],
        ],
        weights=[1, aux_0_w, 1, aux_0_w, 1, aux_0_w, 1, aux_0_w, 1],
    )


def cone(
    outer_radius,
    height,
    inner_radius=None,
    volumetric=True,
    angle=360.0,
    degree=True,
):
    """Creates a cone with circular base.

    Parameters
    ----------
    radius : float
      Radius of the base
    height : float
      Height of the cone
    volumetric : bool, optional
      Parameter whether surface or volume spline, by default True
    angle : float
      Rotation angle in degrees, only used for solid model

    Returns
    -------
    cone: NURBS
      Volumetric or surface NURBS describing a cone
    """

    if volumetric:
        ground = disk(
            outer_radius, inner_radius=inner_radius, angle=angle, degree=degree
        )
    else:
        ground = circle(outer_radius)

    # Extrude in z
    cone = ground.create.extruded([0, 0, height])
    # Move all upper control points to one
    cone.control_points[_np.isclose(cone.control_points[:, -1], height)] = [
        0,
        0,
        height,
    ]

    return cone


def pyramid(width, length, height):
    """Creates a volumetric spline in the shape of a pyramid with linear degree
    in every direction.

    Parameters
    ----------
    width : float
      Dimension of base in x-axis
    length : float
      Dimension of base in y-axis
    height : float
      Height in z-direction

    Returns
    -------
    pyramid: BSpline
      Volumetric linear spline in the shape of a pyramid
    """

    # Create box
    p = box(width, length, height)

    # Collapse all upper points on one control point
    p.control_points[4:, :] = [width / 2, length / 2, height]

    return p


class Creator:
    """Helper class to build new splines from existing geometries.

    Examples
    ---------
    >>> spline_faces = my_spline.create.extrude(vector=[3, 1, 3])
    """

    __slots__ = ("_helpee",)

    def __init__(self, spl):
        self._helpee = spl

    @_wraps(extruded)
    def extruded(self, *args, **kwargs):
        return extruded(self._helpee, *args, **kwargs)

    @_wraps(embedded)
    def embedded(self, *args, **kwargs):
        return embedded(self._helpee, *args, **kwargs)

    @_wraps(revolved)
    def revolved(self, *args, **kwargs):
        return revolved(self._helpee, *args, **kwargs)

    @_wraps(parametric_view)
    def parametric_view(self, *args, **kwargs):
        return parametric_view(self._helpee, *args, **kwargs)

    @_wraps(determinant_spline)
    def determinant_spline(self, *args, **kwargs):
        return determinant_spline(self._helpee, *args, **kwargs)
