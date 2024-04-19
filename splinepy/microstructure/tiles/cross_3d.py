import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class Cross3D(_TileBase):
    """Simple crosstile with linear-quadratic branches and a trilinear
    center spline.

    .. raw:: html

        <p><a href="../_static/Cross3D.html">Fullscreen</a>.</p>
        <embed type="text/html" width="100%" height="400" src="../_static/Cross3D.html" />

    """  # noqa: E501

    _dim = 3
    _para_dim = 3
    _evaluation_points = _np.array(
        [
            [0.0, 0.5, 0.5],
            [1.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 1.0, 0.5],
            [0.5, 0.5, 0.0],
            [0.5, 0.5, 1.0],
        ]
    )
    _n_info_per_eval_point = 1

    def _closing_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        closure=None,
        boundary_width=0.1,
        filling_height=0.5,
        **kwargs,  # noqa ARG002
    ):
        """Create a closing tile to match with closed surface.

        Parameters
        ----------
        parameters : np.ndarray(6, 1)
          Six evaluation points with one parameter is used. This parameter
          describes the radius of the cylinder at the evaluation point.
          The parameters must be a two-dimensional np.array, where the
          value must be between 0.01 and 0.49
        parameter_sensitivities: np.ndarray(6, 1, para_dim)
          Describes the parameter sensitivities with respect to some design
          variable. In case the design variables directly apply to the
          parameter itself, they evaluate as delta_ij.
        boundary_width : float
          with of the boundary surrounding branch
        filling_height : float
          portion of the height that is filled in parametric domain
        closure : str
          parametric dimension that needs to be closed.
          Must be {"z_min", "z_max"}

        Returns
        -------
        list_of_splines : list
        derivative_list : list / None
        """
        # Check parameters
        if closure is None:
            raise ValueError("No closing direction given")

        if parameters is None:
            self._logd("Tile request is not parametrized, setting default 0.2")
            parameters = _np.array(
                _np.ones(
                    (len(self._evaluation_points), self._n_info_per_eval_point)
                )
                * 0.2
            )

        if not (_np.all(parameters > 0) and _np.all(parameters < 0.5)):
            raise ValueError("Thickness out of range (0, .5)")

        # Check if user requests derivative splines
        if parameter_sensitivities is not None:
            self.check_param_derivatives(parameter_sensitivities)
            n_derivatives = parameter_sensitivities.shape[2]
            derivatives = []
        else:
            n_derivatives = 0
            derivatives = None

        if not (0.0 < float(boundary_width) < 0.5):
            raise ValueError("Boundary Width is out of range")

        if not (0.0 < float(filling_height) < 1.0):
            raise ValueError("Filling must  be in (0,1)")

        splines = []
        for i_derivative in range(n_derivatives + 1):
            if i_derivative == 0:
                bound_width = boundary_width
                fill_height = filling_height
                inv_boundary_width = 1.0 - boundary_width
                inv_filling_height = 1.0 - filling_height
                center_width = 1.0 - (2 * boundary_width)
                ctps_mid_height_top = (1.0 + filling_height) * 0.5
                ctps_mid_height_bottom = 1.0 - ctps_mid_height_top
                r_center = center_width * 0.5
                v_zero = 0.0
                v_one_half = 0.5
                v_one = 1.0
            else:
                bound_width = 0.0
                fill_height = 0.0
                inv_boundary_width = 0.0
                inv_filling_height = 0.0
                center_width = 0.0
                ctps_mid_height_top = 0.0
                ctps_mid_height_bottom = 0.0
                r_center = center_width * 0.0
                v_zero = 0.0
                v_one_half = 0.0
                v_one = 0.0

            spline_list = []
            if closure == "z_min":
                branch_thickness = (
                    parameters[5, 0]
                    if i_derivative == 0
                    else parameter_sensitivities[5, 0, i_derivative - 1]
                )

                ctps_corner = _np.array(
                    [
                        [v_zero, v_zero, v_zero],
                        [bound_width, v_zero, v_zero],
                        [v_zero, bound_width, v_zero],
                        [bound_width, bound_width, v_zero],
                        [v_zero, v_zero, fill_height],
                        [bound_width, v_zero, fill_height],
                        [v_zero, bound_width, fill_height],
                        [bound_width, bound_width, fill_height],
                    ]
                )
                spline_list.append(
                    _Bezier(degrees=[1, 1, 1], control_points=ctps_corner)
                )

                spline_list.append(
                    _Bezier(
                        degrees=[1, 1, 1],
                        control_points=(
                            ctps_corner
                            + _np.array([v_zero, inv_boundary_width, v_zero])
                        ),
                    )
                )

                spline_list.append(
                    _Bezier(
                        degrees=[1, 1, 1],
                        control_points=(
                            ctps_corner
                            + _np.array([inv_boundary_width, v_zero, v_zero])
                        ),
                    )
                )

                spline_list.append(
                    _Bezier(
                        degrees=[1, 1, 1],
                        control_points=(
                            ctps_corner
                            + _np.array(
                                [
                                    inv_boundary_width,
                                    inv_boundary_width,
                                    v_zero,
                                ]
                            )
                        ),
                    )
                )

                center_ctps = _np.array(
                    [
                        [bound_width, bound_width, v_zero],
                        [inv_boundary_width, bound_width, v_zero],
                        [bound_width, inv_boundary_width, v_zero],
                        [inv_boundary_width, inv_boundary_width, v_zero],
                        [bound_width, bound_width, fill_height],
                        [inv_boundary_width, bound_width, fill_height],
                        [bound_width, inv_boundary_width, fill_height],
                        [inv_boundary_width, inv_boundary_width, fill_height],
                    ]
                )
                spline_list.append(
                    _Bezier(degrees=[1, 1, 1], control_points=center_ctps)
                )

                spline_list.append(
                    _Bezier(
                        degrees=[1, 1, 1],
                        control_points=(
                            _np.maximum(
                                center_ctps
                                - _np.array([center_width, v_zero, v_zero]),
                                v_zero,
                            )
                            if i_derivative == 0
                            else _np.zeros((8, 3))
                        ),
                    )
                )

                spline_list.append(
                    _Bezier(
                        degrees=[1, 1, 1],
                        control_points=(
                            _np.maximum(
                                center_ctps
                                - _np.array([v_zero, center_width, v_zero]),
                                v_zero,
                            )
                            if i_derivative == 0
                            else _np.zeros((8, 3))
                        ),
                    )
                )

                spline_list.append(
                    _Bezier(
                        degrees=[1, 1, 1],
                        control_points=(
                            _np.minimum(
                                center_ctps
                                + _np.array([center_width, v_zero, v_zero]),
                                v_one,
                            )
                            if i_derivative == 0
                            else _np.zeros((8, 3))
                        ),
                    )
                )

                spline_list.append(
                    _Bezier(
                        degrees=[1, 1, 1],
                        control_points=(
                            _np.minimum(
                                center_ctps
                                + _np.array([v_zero, center_width, v_zero]),
                                v_one,
                            )
                            if i_derivative == 0
                            else _np.zeros((8, 3))
                        ),
                    )
                )

                branch_ctps = _np.array(
                    [
                        [-r_center, -r_center, fill_height],
                        [r_center, -r_center, fill_height],
                        [-r_center, r_center, fill_height],
                        [r_center, r_center, fill_height],
                        [
                            -branch_thickness,
                            -branch_thickness,
                            ctps_mid_height_top,
                        ],
                        [
                            branch_thickness,
                            -branch_thickness,
                            ctps_mid_height_top,
                        ],
                        [
                            -branch_thickness,
                            branch_thickness,
                            ctps_mid_height_top,
                        ],
                        [
                            branch_thickness,
                            branch_thickness,
                            ctps_mid_height_top,
                        ],
                        [-branch_thickness, -branch_thickness, v_one],
                        [branch_thickness, -branch_thickness, v_one],
                        [-branch_thickness, branch_thickness, v_one],
                        [branch_thickness, branch_thickness, v_one],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])
                spline_list.append(
                    _Bezier(degrees=[1, 1, 2], control_points=branch_ctps)
                )
            elif closure == "z_max":
                branch_thickness = (
                    parameters[4, 0]
                    if i_derivative == 0
                    else parameter_sensitivities[4, 0, i_derivative - 1]
                )
                ctps_corner = _np.array(
                    [
                        [v_zero, v_zero, inv_filling_height],
                        [bound_width, v_zero, inv_filling_height],
                        [v_zero, bound_width, inv_filling_height],
                        [bound_width, bound_width, inv_filling_height],
                        [v_zero, v_zero, v_one],
                        [bound_width, v_zero, v_one],
                        [v_zero, bound_width, v_one],
                        [bound_width, bound_width, v_one],
                    ]
                )

                spline_list.append(
                    _Bezier(degrees=[1, 1, 1], control_points=ctps_corner)
                )

                spline_list.append(
                    _Bezier(
                        degrees=[1, 1, 1],
                        control_points=(
                            ctps_corner
                            + _np.array([v_zero, inv_boundary_width, v_zero])
                        ),
                    )
                )

                spline_list.append(
                    _Bezier(
                        degrees=[1, 1, 1],
                        control_points=(
                            ctps_corner
                            + _np.array([inv_boundary_width, v_zero, v_zero])
                        ),
                    )
                )

                spline_list.append(
                    _Bezier(
                        degrees=[1, 1, 1],
                        control_points=(
                            ctps_corner
                            + _np.array(
                                [
                                    inv_boundary_width,
                                    inv_boundary_width,
                                    v_zero,
                                ]
                            )
                        ),
                    )
                )

                center_ctps = _np.array(
                    [
                        [bound_width, bound_width, inv_filling_height],
                        [inv_boundary_width, bound_width, inv_filling_height],
                        [bound_width, inv_boundary_width, inv_filling_height],
                        [
                            inv_boundary_width,
                            inv_boundary_width,
                            inv_filling_height,
                        ],
                        [bound_width, bound_width, v_one],
                        [inv_boundary_width, bound_width, v_one],
                        [bound_width, inv_boundary_width, v_one],
                        [inv_boundary_width, inv_boundary_width, v_one],
                    ]
                )
                spline_list.append(
                    _Bezier(degrees=[1, 1, 1], control_points=center_ctps)
                )

                spline_list.append(
                    _Bezier(
                        degrees=[1, 1, 1],
                        control_points=(
                            _np.maximum(
                                center_ctps
                                - _np.array([center_width, v_zero, v_zero]),
                                v_zero,
                            )
                            if i_derivative == 0
                            else _np.zeros((8, 3))
                        ),
                    )
                )

                spline_list.append(
                    _Bezier(
                        degrees=[1, 1, 1],
                        control_points=(
                            _np.maximum(
                                center_ctps
                                - _np.array([v_zero, center_width, v_zero]),
                                v_zero,
                            )
                            if i_derivative == 0
                            else _np.zeros((8, 3))
                        ),
                    )
                )

                spline_list.append(
                    _Bezier(
                        degrees=[1, 1, 1],
                        control_points=(
                            _np.minimum(
                                center_ctps
                                + _np.array([center_width, v_zero, v_zero]),
                                v_one,
                            )
                            if i_derivative == 0
                            else _np.zeros((8, 3))
                        ),
                    )
                )

                spline_list.append(
                    _Bezier(
                        degrees=[1, 1, 1],
                        control_points=(
                            _np.minimum(
                                center_ctps
                                + _np.array([v_zero, center_width, v_zero]),
                                v_one,
                            )
                            if i_derivative == 0
                            else _np.zeros((8, 3))
                        ),
                    )
                )

                branch_ctps = _np.array(
                    [
                        [-branch_thickness, -branch_thickness, v_zero],
                        [branch_thickness, -branch_thickness, v_zero],
                        [-branch_thickness, branch_thickness, v_zero],
                        [branch_thickness, branch_thickness, v_zero],
                        [
                            -branch_thickness,
                            -branch_thickness,
                            ctps_mid_height_bottom,
                        ],
                        [
                            branch_thickness,
                            -branch_thickness,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -branch_thickness,
                            branch_thickness,
                            ctps_mid_height_bottom,
                        ],
                        [
                            branch_thickness,
                            branch_thickness,
                            ctps_mid_height_bottom,
                        ],
                        [-r_center, -r_center, inv_filling_height],
                        [r_center, -r_center, inv_filling_height],
                        [-r_center, r_center, inv_filling_height],
                        [r_center, r_center, inv_filling_height],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])
                spline_list.append(
                    _Bezier(degrees=[1, 1, 2], control_points=branch_ctps)
                )
            else:
                raise NotImplementedError(
                    "Requested closing dimension is not supported"
                )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        return (splines, derivatives)

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        center_expansion=1.0,
        closure=None,
        **kwargs,  # noqa ARG002
    ):
        """Create a microtile based on the parameters that describe the branch
        thicknesses.

        Thickness parameters are used to describe the inner radius of the
        outward facing branches

        Parameters
        ----------
        parameters : np.array(6, 1)
          Six evaluation points with one parameter is used. This parameter
          describes the radius of the cylinder at the evaluation point.
          The parameters must be a two-dimensional np.array, where the
          value must be between 0.01 and 0.49
        parameter_sensitivities: np.ndarray(6, 1, para_dim)
          Describes the parameter sensitivities with respect to some design
          variable. In case the design variables directly apply to the
          parameter itself, they evaluate as delta_ij.
        center_expansion : float
          thickness of center is expanded by a factor
        closure : str
          parametric dimension that needs to be closed.
          Must be {"z_min", "z_max"}
        **kwargs
          Will be passed to _closing_tile if given

        Returns
        -------
        microtile_list : list(splines)
        derivative_list : list / None
        """

        if not isinstance(center_expansion, float):
            raise ValueError("Invalid Type")

        if not ((center_expansion > 0.5) and (center_expansion < 1.5)):
            raise ValueError("Center Expansion must be in (.5,1.5)")

        # Max radius, so there is no tanglement in the crosstile
        max_radius = min(0.5, (0.5 / center_expansion))

        # set to default if nothing is given
        if parameters is None:
            self._logd("Setting branch thickness to default 0.2")
            parameters = (
                _np.ones(
                    (
                        self._evaluation_points.shape[0],
                        self._n_info_per_eval_point,
                    )
                )
                * 0.2
            )

        # Check for type and consistency
        self.check_params(parameters)
        if _np.any(parameters <= 0) or _np.any(parameters > max_radius):
            raise ValueError(
                f"Radii must be in (0,{max_radius}) for "
                f"center_expansion {center_expansion}"
            )

        if parameter_sensitivities is not None:
            self.check_param_derivatives(parameter_sensitivities)
            n_derivatives = parameter_sensitivities.shape[2]
            derivatives = []
        else:
            n_derivatives = 0
            derivatives = None

        if closure is not None:
            return self._closing_tile(
                parameters=parameters,
                parameter_sensitivities=parameter_sensitivities,
                closure=closure,
                **kwargs,
            )

        splines = []
        for i_derivative in range(n_derivatives + 1):
            # Constant auxiliary values
            if i_derivative == 0:
                [x_min_r, x_max_r, y_min_r, y_max_r, z_min_r, z_max_r] = (
                    parameters[:, 0]
                )
                v_one_half = 0.5
                center_r = center_expansion * _np.mean(parameters[:, 0])
                hd_center = 0.5 * (v_one_half + center_r)
            else:

                [x_min_r, x_max_r, y_min_r, y_max_r, z_min_r, z_max_r] = (
                    parameter_sensitivities[:, :, i_derivative - 1].flatten()
                )
                v_one_half = 0.0
                center_r = center_expansion * _np.mean(
                    parameter_sensitivities[:, :, i_derivative - 1].flatten()
                )
                hd_center = 0.5 * center_r

            spline_list, center = [], [v_one_half] * 3

            # Create the center-tile
            center_points = _np.array(
                [
                    [-center_r, -center_r, -center_r],
                    [center_r, -center_r, -center_r],
                    [-center_r, center_r, -center_r],
                    [center_r, center_r, -center_r],
                    [-center_r, -center_r, center_r],
                    [center_r, -center_r, center_r],
                    [-center_r, center_r, center_r],
                    [center_r, center_r, center_r],
                ]
            )
            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1], control_points=center_points + center
                )
            )

            # X-Axis branches
            # X-Min-Branch
            x_min_ctps = _np.array(
                [
                    [-v_one_half, -x_min_r, -x_min_r],
                    [-hd_center, -x_min_r, -x_min_r],
                    center_points[0, :],
                    [-v_one_half, x_min_r, -x_min_r],
                    [-hd_center, x_min_r, -x_min_r],
                    center_points[2, :],
                    [-v_one_half, -x_min_r, x_min_r],
                    [-hd_center, -x_min_r, x_min_r],
                    center_points[4, :],
                    [-v_one_half, x_min_r, x_min_r],
                    [-hd_center, x_min_r, x_min_r],
                    center_points[6, :],
                ]
            )
            spline_list.append(
                _Bezier(degrees=[2, 1, 1], control_points=x_min_ctps + center)
            )

            # X-Max-Branch
            x_max_ctps = _np.array(
                [
                    center_points[1, :],
                    [hd_center, -x_max_r, -x_max_r],
                    [v_one_half, -x_max_r, -x_max_r],
                    center_points[3, :],
                    [hd_center, x_max_r, -x_max_r],
                    [v_one_half, x_max_r, -x_max_r],
                    center_points[5, :],
                    [hd_center, -x_max_r, x_max_r],
                    [v_one_half, -x_max_r, x_max_r],
                    center_points[7, :],
                    [hd_center, x_max_r, x_max_r],
                    [v_one_half, x_max_r, x_max_r],
                ]
            )
            spline_list.append(
                _Bezier(degrees=[2, 1, 1], control_points=x_max_ctps + center)
            )

            # Y-Axis branches
            # Y-Min-Branch
            y_min_ctps = _np.array(
                [
                    [-y_min_r, -v_one_half, -y_min_r],
                    [y_min_r, -v_one_half, -y_min_r],
                    [-y_min_r, -hd_center, -y_min_r],
                    [y_min_r, -hd_center, -y_min_r],
                    center_points[0, :],
                    center_points[1, :],
                    [-y_min_r, -v_one_half, y_min_r],
                    [y_min_r, -v_one_half, y_min_r],
                    [-y_min_r, -hd_center, y_min_r],
                    [y_min_r, -hd_center, y_min_r],
                    center_points[4, :],
                    center_points[5, :],
                ]
            )
            spline_list.append(
                _Bezier(degrees=[1, 2, 1], control_points=y_min_ctps + center)
            )

            # Y-Max-Branch
            y_max_ctps = _np.array(
                [
                    center_points[2, :],
                    center_points[3, :],
                    [-y_max_r, hd_center, -y_max_r],
                    [y_max_r, hd_center, -y_max_r],
                    [-y_max_r, v_one_half, -y_max_r],
                    [y_max_r, v_one_half, -y_max_r],
                    center_points[6, :],
                    center_points[7, :],
                    [-y_max_r, hd_center, y_max_r],
                    [y_max_r, hd_center, y_max_r],
                    [-y_max_r, v_one_half, y_max_r],
                    [y_max_r, v_one_half, y_max_r],
                ]
            )
            spline_list.append(
                _Bezier(degrees=[1, 2, 1], control_points=y_max_ctps + center)
            )

            # Z-Axis branches
            # Z-Min-Branch
            z_min_ctps = _np.array(
                [
                    [-z_min_r, -z_min_r, -v_one_half],
                    [z_min_r, -z_min_r, -v_one_half],
                    [-z_min_r, z_min_r, -v_one_half],
                    [z_min_r, z_min_r, -v_one_half],
                    [-z_min_r, -z_min_r, -hd_center],
                    [z_min_r, -z_min_r, -hd_center],
                    [-z_min_r, z_min_r, -hd_center],
                    [z_min_r, z_min_r, -hd_center],
                    center_points[0, :],
                    center_points[1, :],
                    center_points[2, :],
                    center_points[3, :],
                ]
            )
            spline_list.append(
                _Bezier(degrees=[1, 1, 2], control_points=z_min_ctps + center)
            )

            # Z-Max-Branch
            z_max_ctps = _np.array(
                [
                    center_points[4, :],
                    center_points[5, :],
                    center_points[6, :],
                    center_points[7, :],
                    [-z_max_r, -z_max_r, hd_center],
                    [z_max_r, -z_max_r, hd_center],
                    [-z_max_r, z_max_r, hd_center],
                    [z_max_r, z_max_r, hd_center],
                    [-z_max_r, -z_max_r, v_one_half],
                    [z_max_r, -z_max_r, v_one_half],
                    [-z_max_r, z_max_r, v_one_half],
                    [z_max_r, z_max_r, v_one_half],
                ]
            )
            spline_list.append(
                _Bezier(degrees=[1, 1, 2], control_points=z_max_ctps + center)
            )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        return (splines, derivatives)
