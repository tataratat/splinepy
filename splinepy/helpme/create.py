import numpy as np
from gustaf.utils import arr

from splinepy import settings
from splinepy.utils import log


def extruded(spline, extrusion_vector=None):
    """Extrudes Splines.

    Parameters
    ----------
    spline: Spline
    extrusion_vector: np.ndarray
    """
    from splinepy.spline import Spline

    # Check input type
    if not isinstance(spline, Spline):
        raise NotImplementedError("Extrude only works for splines")

    # Check extrusion_vector
    if extrusion_vector is not None:
        # make flat extrusion_vector
        extrusion_vector = np.asarray(extrusion_vector).ravel()
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
        log.debug(
            f"Given extrusion vector is {expansion_dimension} dimension "
            "bigger than spline's dim. Assuming 0.0 entries for "
            "new dimension.",
        )
        cps = np.hstack(
            (
                spline.control_points,
                np.zeros((len(spline.control_points), expansion_dimension)),
            )
        )
    else:
        raise ValueError(
            "Dimension Mismatch between extrusion extrusion vector "
            "and spline."
        )

    # Start Extrusion
    spline_dict = dict()

    spline_dict["degrees"] = np.concatenate((spline.degrees, [1]))
    spline_dict["control_points"] = np.vstack((cps, cps + extrusion_vector))
    if spline.has_knot_vectors:
        spline_dict["knot_vectors"] = spline.knot_vectors + [[0, 0, 1, 1]]
    if spline.is_rational:
        spline_dict["weights"] = np.concatenate(
            (spline.weights, spline.weights)
        )

    return type(spline)(**spline_dict)


def revolved(
    spline, axis=None, center=None, angle=None, n_knot_spans=None, degree=True
):
    """Revolve spline around an axis and extend its parametric dimension.

    Parameters
    ----------
    spline : Spline
      Basis-Spline to be revolved
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
    from splinepy.spline import Spline

    # Check input type
    if not isinstance(spline, Spline):
        raise NotImplementedError("Revolutions only works for splines")

    # Check axis
    if axis is not None:
        # Transform into numpy array
        axis = np.asarray(axis).ravel()
        # Check Axis dimension
        if spline.control_points.shape[1] > axis.shape[0]:
            raise ValueError(
                "Dimension Mismatch between extrusion axis and spline."
            )
        elif spline.control_points.shape[1] < axis.shape[0]:
            log.debug(
                "Control Point dimension is smaller than axis dimension,"
                " filling with zeros"
            )
            expansion_dimension = axis.shape[0] - spline.dim
            cps = np.hstack(
                (
                    spline.control_points,
                    np.zeros(
                        (len(spline.control_points), expansion_dimension)
                    ),
                )
            )
        else:
            cps = np.copy(spline.control_points)

        # Make sure axis is normalized

        axis_norm = np.linalg.norm(axis)
        if not np.isclose(axis_norm, 0, atol=settings.TOLERANCE):
            axis = axis / axis_norm
        else:
            raise ValueError("Axis-norm is too close to zero.")
    else:
        cps = np.copy(spline.control_points)
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
        angle = 360

    if degree:
        angle = np.radians(angle)

    # Init center
    if center is not None:
        center = np.asarray(center).ravel()
        # Check Axis dimension
        if not (problem_dimension == center.shape[0]):
            raise ValueError(
                "Dimension Mismatch between axis and center of rotation."
            )
        cps -= center

    # The parametric dimension is independent of the revolution but the
    # rotation-matrix is only implemented for 2D and 3D problems
    if not (cps.shape[1] == 2 or cps.shape[1] == 3):
        raise NotImplementedError(
            "Sorry," "revolutions only implemented for 2D and 3D splines"
        )

    # Angle must be (0, pi) non including
    # Rotation is always performed in half steps
    PI = np.pi
    minimum_n_knot_spans = int(
        np.ceil(np.abs((angle + settings.TOLERANCE) / PI))
    )
    if (n_knot_spans) is None or (n_knot_spans < minimum_n_knot_spans):
        n_knot_spans = minimum_n_knot_spans

    if "Bezier" in spline.name:
        if n_knot_spans > 1:
            raise ValueError(
                "Revolutions are only supported for angles up to 180 "
                "degrees for Bezier type splines as they consist of only "
                "one knot span"
            )

    # Determine auxiliary values
    rot_a = angle / (2 * n_knot_spans)
    half_counter_angle = PI / 2 - rot_a
    weight = np.sin(half_counter_angle)
    factor = 1 / weight

    # Determine rotation matrix
    rotation_matrix = arr.rotation_matrix_around_axis(
        axis=axis, rotation=rot_a, degree=False
    ).T

    # Start Extrusion
    spline_dict = dict()

    spline_dict["degrees"] = np.concatenate((spline.degrees, [2]))

    spline_dict["control_points"] = cps
    end_points = cps
    for i_segment in range(n_knot_spans):
        # Rotate around axis
        mid_points = np.matmul(end_points, rotation_matrix)
        end_points = np.matmul(mid_points, rotation_matrix)
        # Move away from axis using dot-product tricks
        if problem_dimension == 3:
            mp_scale = axis * np.dot(mid_points, axis).reshape(-1, 1)
            mid_points = (mid_points - mp_scale) * factor + mp_scale
        else:
            mid_points *= factor
        spline_dict["control_points"] = np.concatenate(
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
        for i_segment in range(n_knot_spans):
            spline_dict["weights"] = np.concatenate(
                (spline_dict["weights"], mid_weights, spline.weights)
            )
    else:
        log.debug(
            "True revolutions are only possible for rational spline types.",
            "Creating Approximation.",
        )

    if center is not None:
        spline_dict["control_points"] += center

    return type(spline)(**spline_dict)


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
    physical_bounds = np.asanyarray(physical_bounds).reshape(2, -1)
    parametric_bounds = np.asanyarray(parametric_bounds).reshape(2, -1)

    # get correctly sized bez box
    phys_size = physical_bounds[1] - physical_bounds[0]

    # minimal bspline
    bspline_box = box(*phys_size).bspline  # kvs are in [0, 1]
    bspline_box.cps += physical_bounds[0]  # apply offset

    # update parametric bounds
    for i, kv in enumerate(bspline_box.kvs):
        # apply scale and offset
        new_kv = (kv * parametric_bounds[1][i]) + parametric_bounds[0][i]
        bspline_box.kvs[i] = new_kv

    return bspline_box


def parametric_view(spline, axes=True):
    """Create parametric view of given spline. Previously called
    `naive_spline()`. Degrees are always 1 and knot multiplicity is not
    preserved. Returns BSpline, as BSpline and NURBS should look the same as
    parametric view.
    Will take shallow copy of underlying data of spline_data and show_options
    from original spline.

    Parameters
    -----------
    spline: BSpline or NURBS
    axes: bool
     If True, will configure axes settings, it is supported.

    Returns
    --------
    para_spline: BSpline
    """
    p_bounds = spline.parametric_bounds
    para_spline = from_bounds(
        parametric_bounds=p_bounds, physical_bounds=p_bounds
    )

    # loop through knot vectors and insert missing knots
    for i, kv in enumerate(spline.unique_knots):
        if len(kv) > 2:
            para_spline.insert_knots(i, kv[1:-1])

    # take shallow copy
    para_spline.spline_data._saved = spline.spline_data._saved.copy()
    para_spline.show_options._options[
        para_spline._show_options._backend
    ] = spline.show_options._options[spline._show_options._backend].copy()

    if axes and "axes" in spline.show_options.valid_keys():
        # configure axes
        bs = p_bounds
        bs_diff_001 = (bs[1] - bs[0]) * 0.001
        lower_b = bs[0] - bs_diff_001
        upper_b = bs[1] + bs_diff_001
        axes_config = dict(
            xtitle="u",
            ytitle="v",
            xrange=[lower_b[0], upper_b[0]],
            yrange=[lower_b[1], upper_b[1]],
            tip_size=0,
            xminor_ticks=3,
            yminor_ticks=3,
            xygrid=False,
            yzgrid=False,
        )
        if spline.para_dim == 3:
            axes_config.update(ztitle="w")
            axes_config.update(zrange=[lower_b[2], upper_b[2]])
            axes_config.update(zminor_ticks=3)
            axes_config.update(zxgrid=False)

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

    cps = np.array(points)
    nknots = cps.shape[0] + degree + 1

    knots = np.concatenate(
        (
            np.full(degree, 0.0),
            np.linspace(0.0, 1.0, nknots - 2 * degree),
            np.full(degree, 1.0),
        )
    )

    spline = settings.NAME_TO_TYPE["BSpline"](
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
        start_angle = np.radians(start_angle)
        angle = np.radians(angle)
    start_point = [radius * np.cos(start_angle), radius * np.sin(start_angle)]
    point_spline = settings.NAME_TO_TYPE["RationalBezier"](
        degrees=[0], control_points=[start_point], weights=[1.0]
    )
    # Bezier splines only support angles lower than 180 degrees
    if abs(angle) >= np.pi or n_knot_spans > 1:
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
    """ND box (hyper rectangle).

    Parameters
    ----------
    *lengths: list(float)

    Returns
    -------
    nd_box: Bezier
    """
    # may dim check here?
    # starting point
    nd_box = settings.NAME_TO_TYPE["Bezier"](
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
        np.array(
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
    weights = np.tile([1.0, 1 / np.sqrt(2)], 5)[:-1]

    return settings.NAME_TO_TYPE["RationalBezier"](
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

    cps = np.array([[inner_radius, 0.0], [outer_radius, 0.0]])
    weights = np.ones([cps.shape[0]])
    knots = np.repeat([0.0, 1.0], 2)

    return settings.NAME_TO_TYPE["NURBS"](
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
        torus_angle = 2 * np.pi
        degree = False

    if section_angle is None:
        section_angle = 2 * np.pi
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
        if abs(torus_angle) >= np.pi:
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
        center=np.zeros(3),
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

    return settings.NAME_TO_TYPE["RationalBezier"](
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
    cone.control_points[np.isclose(cone.control_points[:, -1], height)] = [
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
    >>> my_spline = <your-spline>
    >>> spline_faces = my_spline.create.extrude(vector=[3,1,3])
    """

    def __init__(self, spl):
        self.spline = spl

    def extruded(self, *args, **kwargs):
        return extruded(self.spline, *args, **kwargs)

    def revolved(self, *args, **kwargs):
        return revolved(self.spline, *args, **kwargs)

    def parametric_view(self, *args, **kwargs):
        return parametric_view(self.spline, *args, **kwargs)
