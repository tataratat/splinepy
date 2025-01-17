import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class Cross2D(_TileBase):
    """Simple crosstile with linear-quadratic branches and a bilinear
    center spline.

    .. raw:: html

        <p><a href="../_static/Cross2D.html">Fullscreen</a>.</p>
        <embed type="text/html" width="100%" height="400" src="../_static/Cross2D.html" />

    """  # noqa: E501

    _dim = 2
    _para_dim = 2
    _evaluation_points = _np.array(
        [
            [0.0, 0.5],
            [1.0, 0.5],
            [0.5, 0.0],
            [0.5, 1.0],
        ]
    )
    _n_info_per_eval_point = 1
    _sensitivities_implemented = True
    _closure_directions = ["x_min", "x_max", "y_min", "y_max"]
    _parameter_bounds = [
        [0.0, 0.5]
    ] * 4  # valid for default center_expansion=1.0
    _parameters_shape = (4, 1)
    _default_parameter_value = 0.2

    _BOUNDARY_WIDTH_BOUNDS = [0.0, 0.5]
    _FILLING_HEIGHT_BOUNDS = [0.0, 1.0]
    _CENTER_EXPANSION_BOUNDS = [0.5, 1.5]

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
        parameters : np.ndarray(4, 1)
          Four evaluation poinst with one parameter is used. This parameter
          describes the radius of the cylinder at the evaluation point.
          The parameters must be a two-dimensional np.array, where the
          value must be between 0.01 and 0.49
        parameter_sensitivities: np.ndarray
          Describes the parameter sensitivities with respect to some design
          variable. In case the design variables directly apply to the
          parameter itself, they evaluate as delta_ij
        boundary_width : float
          with of the boundary surrounding branch
        filling_height : float
          portion of the height that is filled in parametric domain
        closure : str
          parametric dimension that needs to be closed, given in the form
          "x_min", "x_max", etc.

        Returns
        -------
        list_of_splines : list
        derivative_list : list / None
        """
        # Check parameters
        if closure is None:
            raise ValueError("No closing direction given")

        parameters, n_derivatives, derivatives = self._process_input(
            parameters=parameters,
            parameter_sensitivities=parameter_sensitivities,
        )

        self._check_custom_parameter(
            boundary_width, "boundary width", self._BOUNDARY_WIDTH_BOUNDS
        )
        self._check_custom_parameter(
            filling_height, "filling height", self._FILLING_HEIGHT_BOUNDS
        )

        splines = []
        for i_derivative in range(n_derivatives + 1):
            # Constant auxiliary values
            if i_derivative == 0:
                fill_height = filling_height
                bound_width = boundary_width
                inv_bound_width = 1.0 - bound_width
                inv_fill_height = 1.0 - fill_height
                ctps_mid_height_top = (1 + fill_height) * 0.5
                ctps_mid_height_bottom = 1.0 - ctps_mid_height_top
                v_one_half = 0.5
                v_one = 1.0
                v_zero = 0.0
                parameters = parameters[:, 0]
            else:
                # Set constant values to zero for derivatives
                fill_height = 0.0
                bound_width = 0.0
                inv_bound_width = 0.0
                inv_fill_height = 0.0
                ctps_mid_height_top = 0.0
                ctps_mid_height_bottom = 0.0
                v_one_half = 0.0
                v_one = 0.0
                v_zero = 0.0
                parameters = parameter_sensitivities[:, 0, i_derivative - 1]

            spline_list = []
            if closure == "x_min":
                # Minimum x position
                branch_thickness = parameters[1]

                block0_ctps = _np.array(
                    [
                        [v_zero, v_zero],
                        [fill_height, v_zero],
                        [v_zero, bound_width],
                        [fill_height, bound_width],
                    ]
                )

                block1_ctps = _np.array(
                    [
                        [v_zero, bound_width],
                        [fill_height, bound_width],
                        [v_zero, inv_bound_width],
                        [fill_height, inv_bound_width],
                    ]
                )

                block2_ctps = _np.array(
                    [
                        [v_zero, inv_bound_width],
                        [fill_height, inv_bound_width],
                        [v_zero, v_one],
                        [fill_height, v_one],
                    ]
                )

                branch_ctps = _np.array(
                    [
                        [fill_height, bound_width],
                        [ctps_mid_height_top, v_one_half - branch_thickness],
                        [v_one, v_one_half - branch_thickness],
                        [fill_height, inv_bound_width],
                        [ctps_mid_height_top, v_one_half + branch_thickness],
                        [v_one, v_one_half + branch_thickness],
                    ]
                )

                spline_list.append(
                    _Bezier(degrees=[1, 1], control_points=block0_ctps)
                )

                spline_list.append(
                    _Bezier(degrees=[1, 1], control_points=block1_ctps)
                )

                spline_list.append(
                    _Bezier(degrees=[1, 1], control_points=block2_ctps)
                )

                spline_list.append(
                    _Bezier(degrees=[2, 1], control_points=branch_ctps)
                )
            elif closure == "x_max":
                # Maximum x position
                branch_thickness = parameters[0]

                block0_ctps = _np.array(
                    [
                        [inv_fill_height, v_zero],
                        [v_one, v_zero],
                        [inv_fill_height, bound_width],
                        [v_one, bound_width],
                    ]
                )

                block1_ctps = _np.array(
                    [
                        [inv_fill_height, bound_width],
                        [v_one, bound_width],
                        [inv_fill_height, inv_bound_width],
                        [v_one, inv_bound_width],
                    ]
                )

                block2_ctps = _np.array(
                    [
                        [inv_fill_height, inv_bound_width],
                        [v_one, inv_bound_width],
                        [inv_fill_height, v_one],
                        [v_one, v_one],
                    ]
                )

                branch_ctps = _np.array(
                    [
                        [0, v_one_half - branch_thickness],
                        [
                            ctps_mid_height_bottom,
                            v_one_half - branch_thickness,
                        ],
                        [inv_fill_height, bound_width],
                        [v_zero, v_one_half + branch_thickness],
                        [
                            ctps_mid_height_bottom,
                            v_one_half + branch_thickness,
                        ],
                        [inv_fill_height, inv_bound_width],
                    ]
                )

                spline_list.append(
                    _Bezier(degrees=[1, 1], control_points=block0_ctps)
                )

                spline_list.append(
                    _Bezier(degrees=[1, 1], control_points=block1_ctps)
                )

                spline_list.append(
                    _Bezier(degrees=[1, 1], control_points=block2_ctps)
                )

                spline_list.append(
                    _Bezier(degrees=[2, 1], control_points=branch_ctps)
                )
            elif closure == "y_min":
                # Minimum y position
                branch_thickness = parameters[3]

                block0_ctps = _np.array(
                    [
                        [v_zero, v_zero],
                        [bound_width, v_zero],
                        [v_zero, fill_height],
                        [bound_width, fill_height],
                    ]
                )

                block1_ctps = _np.array(
                    [
                        [bound_width, v_zero],
                        [inv_bound_width, v_zero],
                        [bound_width, fill_height],
                        [inv_bound_width, fill_height],
                    ]
                )

                block2_ctps = _np.array(
                    [
                        [inv_bound_width, v_zero],
                        [v_one, v_zero],
                        [inv_bound_width, fill_height],
                        [v_one, fill_height],
                    ]
                )

                branch_ctps = _np.array(
                    [
                        [bound_width, fill_height],
                        [inv_bound_width, fill_height],
                        [v_one_half - branch_thickness, ctps_mid_height_top],
                        [v_one_half + branch_thickness, ctps_mid_height_top],
                        [v_one_half - branch_thickness, v_one],
                        [v_one_half + branch_thickness, v_one],
                    ]
                )

                spline_list.append(
                    _Bezier(degrees=[1, 1], control_points=block0_ctps)
                )

                spline_list.append(
                    _Bezier(degrees=[1, 1], control_points=block1_ctps)
                )

                spline_list.append(
                    _Bezier(degrees=[1, 1], control_points=block2_ctps)
                )

                spline_list.append(
                    _Bezier(degrees=[1, 2], control_points=branch_ctps)
                )
            elif closure == "y_max":
                # Maximum y position
                branch_thickness = parameters[2]

                block0_ctps = _np.array(
                    [
                        [v_zero, inv_fill_height],
                        [bound_width, inv_fill_height],
                        [v_zero, v_one],
                        [bound_width, v_one],
                    ]
                )

                block1_ctps = _np.array(
                    [
                        [bound_width, inv_fill_height],
                        [inv_bound_width, inv_fill_height],
                        [bound_width, v_one],
                        [inv_bound_width, v_one],
                    ]
                )

                block2_ctps = _np.array(
                    [
                        [inv_bound_width, inv_fill_height],
                        [v_one, inv_fill_height],
                        [inv_bound_width, v_one],
                        [v_one, v_one],
                    ]
                )

                branch_ctps = _np.array(
                    [
                        [v_one_half - branch_thickness, v_zero],
                        [v_one_half + branch_thickness, v_zero],
                        [
                            v_one_half - branch_thickness,
                            ctps_mid_height_bottom,
                        ],
                        [
                            v_one_half + branch_thickness,
                            ctps_mid_height_bottom,
                        ],
                        [bound_width, inv_fill_height],
                        [inv_bound_width, inv_fill_height],
                    ]
                )

                spline_list.append(
                    _Bezier(degrees=[1, 1], control_points=block0_ctps)
                )

                spline_list.append(
                    _Bezier(degrees=[1, 1], control_points=block1_ctps)
                )

                spline_list.append(
                    _Bezier(degrees=[1, 1], control_points=block2_ctps)
                )

                spline_list.append(
                    _Bezier(degrees=[1, 2], control_points=branch_ctps)
                )
            else:
                raise NotImplementedError(
                    "Requested closing dimension is not supported"
                )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)
        # Return results
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
        parameters : np.array(4, 1)
          Four evaluation point with one parameter is used. This parameter
          describes the radius of the cylinder at the evaluation point.
          The parameters must be a two-dimensional np.array, where the
          value must be between 0.01 and 0.49
        parameter_sensitivities: list(np.ndarray)
          Describes the parameter sensitivities with respect to some design
          variable. In case the design variables directly apply to the
          parameter itself, they evaluate as delta_ij
        center_expansion : float
            thickness of center is expanded by a factor
        closure : str
          parametric dimension that needs to be closed, given in the form
          "x_min", "x_max", etc.

        Returns
        -------
        microtile_list : list(splines)
        derivative_list : list / None
        """

        self._check_custom_parameter(
            center_expansion, "center expansion", self._CENTER_EXPANSION_BOUNDS
        )

        parameters, n_derivatives, derivatives = self._process_input(
            parameters=parameters,
            parameter_sensitivities=parameter_sensitivities,
        )

        max_radius = min(0.5, (0.5 / center_expansion))
        if not (_np.all(parameters > 0) and _np.all(parameters < max_radius)):
            raise ValueError(f"Thickness out of range (0, {max_radius})")

        # Closure requested, pass to function
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
                [x_min_r, x_max_r, y_min_r, y_max_r] = parameters[:, 0]
                v_one_half = 0.5
                # center radius
                center_r = (
                    (x_min_r + x_max_r + y_min_r + y_max_r)
                    / 4.0
                    * center_expansion
                ).item()
                hd_center = 0.5 * (0.5 + center_r)
            else:
                [x_min_r, x_max_r, y_min_r, y_max_r] = (
                    parameter_sensitivities[:, :, i_derivative - 1]
                ).flatten()
                v_one_half = 0.0
                # center radius
                center_r = (
                    (x_min_r + x_max_r + y_min_r + y_max_r)
                    / 4.0
                    * center_expansion
                )
                hd_center = 0.5 * center_r

            # Init return value
            spline_list = []

            # Create the center-tile
            center_points = _np.array(
                [
                    [-center_r, -center_r],
                    [center_r, -center_r],
                    [-center_r, center_r],
                    [center_r, center_r],
                ]
            ) + _np.array([v_one_half, v_one_half])

            y_min_ctps = _np.array(
                [
                    [-y_min_r, -v_one_half],
                    [y_min_r, -v_one_half],
                    [-y_min_r, -hd_center],
                    [y_min_r, -hd_center],
                    [-center_r, -center_r],
                    [center_r, -center_r],
                ]
            ) + _np.array([v_one_half, v_one_half])

            y_max_ctps = _np.array(
                [
                    [-center_r, center_r],
                    [center_r, center_r],
                    [-y_max_r, hd_center],
                    [y_max_r, hd_center],
                    [-y_max_r, v_one_half],
                    [y_max_r, v_one_half],
                ]
            ) + _np.array([v_one_half, v_one_half])

            x_min_ctps = _np.array(
                [
                    [-v_one_half, -x_min_r],
                    [-hd_center, -x_min_r],
                    [-center_r, -center_r],
                    [-v_one_half, x_min_r],
                    [-hd_center, x_min_r],
                    [-center_r, center_r],
                ]
            ) + _np.array([v_one_half, v_one_half])

            x_max_ctps = _np.array(
                [
                    [center_r, -center_r],
                    [hd_center, -x_max_r],
                    [v_one_half, -x_max_r],
                    [center_r, center_r],
                    [hd_center, x_max_r],
                    [v_one_half, x_max_r],
                ]
            ) + _np.array([v_one_half, v_one_half])

            spline_list.append(
                _Bezier(degrees=[1, 1], control_points=center_points)
            )

            spline_list.append(
                _Bezier(degrees=[2, 1], control_points=x_min_ctps)
            )

            spline_list.append(
                _Bezier(degrees=[2, 1], control_points=x_max_ctps)
            )

            spline_list.append(
                _Bezier(degrees=[1, 2], control_points=y_min_ctps)
            )

            spline_list.append(
                _Bezier(degrees=[1, 2], control_points=y_max_ctps)
            )

            # Pass to output
            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)
        # Return results
        return (splines, derivatives)
