import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class InverseCross3D(_TileBase):
    """Simple inverse crosstile to tile with linear-quadratic branches and
    a trilinear center spline.

    Class that provides necessary functions to create inverse microtile,
    that can be used to describe the domain surrounding cross-tile
    microstructure.

    .. raw:: html

        <p><a href="../_static/InverseCross3D.html">Fullscreen</a>.</p>
        <embed type="text/html" width="100%" height="400" src="../_static/InverseCross3D.html" />

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
    # TODO: implemented sensitivities are not correct
    _sensitivities_implemented = True
    _closure_directions = ["z_min", "z_max"]
    _parameter_bounds = [[0.2, 0.3]] * 6  # For default values
    _parameters_shape = (6, 1)
    _default_parameter_value = 0.21

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
        separator_distance=0.3,
        **kwargs,  # noqa ARG002
    ):
        """Create a closing tile to match with closed surface.

        Parameters
        ----------
        parameters : tuple(np.ndarray)
          radii of fitting cylinder at evaluation points
        boundary_width : float
          with of the boundary surrounding branch
        filling_height : float
          portion of the height that is filled in parametric domain
        separator_distance : float
          Describes the position of the separator layer in the control point
          domain
        closure : str
          parametric dimension that needs to be closed.
          Must be one of {"z_min", "z_max"}

        Returns
        -------
        list_of_splines : list
        derivatives: list<list<splines>> / None
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
            if i_derivative == 0:
                # Auxiliary values
                fill_height_aux = filling_height
                sep_distance_aux = separator_distance
                inv_filling_height = 1.0 - filling_height
                ctps_mid_height_top = (1 + filling_height) * 0.5
                ctps_mid_height_bottom = 1.0 - ctps_mid_height_top
                center_width = 1.0 - 2 * boundary_width
                r_center = center_width * 0.5
                half_r_center = (r_center + 0.5) * 0.5
                aux_column_width = 0.5 - 2 * (0.5 - separator_distance)
                v_zero = 0.0
                v_one_half = 0.5
                v_one = 1.0
                if closure == "z_min":
                    branch_thickness = parameters.flatten()[5]
                elif closure == "z_max":
                    branch_thickness = parameters.flatten()[4]
            else:
                fill_height_aux = 0.0
                sep_distance_aux = 0.0
                inv_filling_height = 0.0
                ctps_mid_height_top = 0.0
                ctps_mid_height_bottom = 0.0
                center_width = 0.0
                r_center = 0.0
                half_r_center = 0.0
                aux_column_width = 0.0
                v_zero, v_one_half, v_one = [0.0] * 3
                if closure == "z_min":
                    branch_thickness = parameter_sensitivities.flatten()[5]
                elif closure == "z_max":
                    branch_thickness = parameter_sensitivities.flatten()[4]
            spline_list = []
            if closure == "z_min":
                branch_neighbor_x_min_ctps = _np.array(
                    [
                        [-v_one_half, -r_center, fill_height_aux],
                        [-half_r_center, -r_center, fill_height_aux],
                        [-r_center, -r_center, fill_height_aux],
                        [-v_one_half, r_center, fill_height_aux],
                        [-half_r_center, r_center, fill_height_aux],
                        [-r_center, r_center, fill_height_aux],
                        [-v_one_half, -aux_column_width, ctps_mid_height_top],
                        [
                            -sep_distance_aux,
                            -aux_column_width,
                            ctps_mid_height_top,
                        ],
                        [
                            -branch_thickness,
                            -branch_thickness,
                            ctps_mid_height_top,
                        ],
                        [-v_one_half, aux_column_width, ctps_mid_height_top],
                        [
                            -sep_distance_aux,
                            aux_column_width,
                            ctps_mid_height_top,
                        ],
                        [
                            -branch_thickness,
                            branch_thickness,
                            ctps_mid_height_top,
                        ],
                        [-v_one_half, -aux_column_width, v_one],
                        [-sep_distance_aux, -aux_column_width, v_one],
                        [-branch_thickness, -branch_thickness, v_one],
                        [-v_one_half, aux_column_width, v_one],
                        [-sep_distance_aux, aux_column_width, v_one],
                        [-branch_thickness, branch_thickness, v_one],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[2, 1, 2],
                        control_points=branch_neighbor_x_min_ctps,
                    )
                )

                branch_neighbor_x_max_ctps = _np.array(
                    [
                        [r_center, -r_center, fill_height_aux],
                        [half_r_center, -r_center, fill_height_aux],
                        [v_one_half, -r_center, fill_height_aux],
                        [r_center, r_center, fill_height_aux],
                        [half_r_center, r_center, fill_height_aux],
                        [v_one_half, r_center, fill_height_aux],
                        [
                            branch_thickness,
                            -branch_thickness,
                            ctps_mid_height_top,
                        ],
                        [
                            sep_distance_aux,
                            -aux_column_width,
                            ctps_mid_height_top,
                        ],
                        [v_one_half, -aux_column_width, ctps_mid_height_top],
                        [
                            branch_thickness,
                            branch_thickness,
                            ctps_mid_height_top,
                        ],
                        [
                            sep_distance_aux,
                            aux_column_width,
                            ctps_mid_height_top,
                        ],
                        [v_one_half, aux_column_width, ctps_mid_height_top],
                        [branch_thickness, -branch_thickness, v_one],
                        [sep_distance_aux, -aux_column_width, v_one],
                        [v_one_half, -aux_column_width, v_one],
                        [branch_thickness, branch_thickness, v_one],
                        [sep_distance_aux, aux_column_width, v_one],
                        [v_one_half, aux_column_width, v_one],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[2, 1, 2],
                        control_points=branch_neighbor_x_max_ctps,
                    )
                )

                branch_neighbor_y_min_ctps = _np.array(
                    [
                        [-r_center, -v_one_half, fill_height_aux],
                        [r_center, -v_one_half, fill_height_aux],
                        [-r_center, -half_r_center, fill_height_aux],
                        [r_center, -half_r_center, fill_height_aux],
                        [-r_center, -r_center, fill_height_aux],
                        [r_center, -r_center, fill_height_aux],
                        [-aux_column_width, -v_one_half, ctps_mid_height_top],
                        [aux_column_width, -v_one_half, ctps_mid_height_top],
                        [
                            -aux_column_width,
                            -sep_distance_aux,
                            ctps_mid_height_top,
                        ],
                        [
                            aux_column_width,
                            -sep_distance_aux,
                            ctps_mid_height_top,
                        ],
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
                        [-aux_column_width, -v_one_half, v_one],
                        [aux_column_width, -v_one_half, v_one],
                        [-aux_column_width, -sep_distance_aux, v_one],
                        [aux_column_width, -sep_distance_aux, v_one],
                        [-branch_thickness, -branch_thickness, v_one],
                        [branch_thickness, -branch_thickness, v_one],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[1, 2, 2],
                        control_points=branch_neighbor_y_min_ctps,
                    )
                )

                branch_neighbor_y_max_ctps = _np.array(
                    [
                        [-r_center, r_center, fill_height_aux],
                        [r_center, r_center, fill_height_aux],
                        [-r_center, half_r_center, fill_height_aux],
                        [r_center, half_r_center, fill_height_aux],
                        [-r_center, v_one_half, fill_height_aux],
                        [r_center, v_one_half, fill_height_aux],
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
                        [
                            -aux_column_width,
                            sep_distance_aux,
                            ctps_mid_height_top,
                        ],
                        [
                            aux_column_width,
                            sep_distance_aux,
                            ctps_mid_height_top,
                        ],
                        [-aux_column_width, v_one_half, ctps_mid_height_top],
                        [aux_column_width, v_one_half, ctps_mid_height_top],
                        [-branch_thickness, branch_thickness, v_one],
                        [branch_thickness, branch_thickness, v_one],
                        [-aux_column_width, sep_distance_aux, v_one],
                        [aux_column_width, sep_distance_aux, v_one],
                        [-aux_column_width, v_one_half, v_one],
                        [aux_column_width, v_one_half, v_one],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[1, 2, 2],
                        control_points=branch_neighbor_y_max_ctps,
                    )
                )

                branch_x_min_y_min_ctps = _np.array(
                    [
                        [-v_one_half, -v_one_half, fill_height_aux],
                        [-half_r_center, -v_one_half, fill_height_aux],
                        [-r_center, -v_one_half, fill_height_aux],
                        [-v_one_half, -half_r_center, fill_height_aux],
                        [-half_r_center, -half_r_center, fill_height_aux],
                        [-r_center, -half_r_center, fill_height_aux],
                        [-v_one_half, -r_center, fill_height_aux],
                        [-half_r_center, -r_center, fill_height_aux],
                        [-r_center, -r_center, fill_height_aux],
                        [-v_one_half, -v_one_half, ctps_mid_height_top],
                        [-sep_distance_aux, -v_one_half, ctps_mid_height_top],
                        [-aux_column_width, -v_one_half, ctps_mid_height_top],
                        [-v_one_half, -sep_distance_aux, ctps_mid_height_top],
                        [
                            -sep_distance_aux,
                            -sep_distance_aux,
                            ctps_mid_height_top,
                        ],
                        [
                            -aux_column_width,
                            -sep_distance_aux,
                            ctps_mid_height_top,
                        ],
                        [-v_one_half, -aux_column_width, ctps_mid_height_top],
                        [
                            -sep_distance_aux,
                            -aux_column_width,
                            ctps_mid_height_top,
                        ],
                        [
                            -branch_thickness,
                            -branch_thickness,
                            ctps_mid_height_top,
                        ],
                        [-v_one_half, -v_one_half, v_one],
                        [-sep_distance_aux, -v_one_half, v_one],
                        [-aux_column_width, -v_one_half, v_one],
                        [-v_one_half, -sep_distance_aux, v_one],
                        [-sep_distance_aux, -sep_distance_aux, v_one],
                        [-aux_column_width, -sep_distance_aux, v_one],
                        [-v_one_half, -aux_column_width, v_one],
                        [-sep_distance_aux, -aux_column_width, v_one],
                        [-branch_thickness, -branch_thickness, v_one],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[2, 2, 2],
                        control_points=branch_x_min_y_min_ctps,
                    )
                )

                branch_x_min_y_max_ctps = _np.array(
                    [
                        [-v_one_half, r_center, fill_height_aux],
                        [-half_r_center, r_center, fill_height_aux],
                        [-r_center, r_center, fill_height_aux],
                        [-v_one_half, half_r_center, fill_height_aux],
                        [-half_r_center, half_r_center, fill_height_aux],
                        [-r_center, half_r_center, fill_height_aux],
                        [-v_one_half, v_one_half, fill_height_aux],
                        [-half_r_center, v_one_half, fill_height_aux],
                        [-r_center, v_one_half, fill_height_aux],
                        [-v_one_half, aux_column_width, ctps_mid_height_top],
                        [
                            -sep_distance_aux,
                            aux_column_width,
                            ctps_mid_height_top,
                        ],
                        [
                            -branch_thickness,
                            branch_thickness,
                            ctps_mid_height_top,
                        ],
                        [-v_one_half, sep_distance_aux, ctps_mid_height_top],
                        [
                            -sep_distance_aux,
                            sep_distance_aux,
                            ctps_mid_height_top,
                        ],
                        [
                            -aux_column_width,
                            sep_distance_aux,
                            ctps_mid_height_top,
                        ],
                        [-v_one_half, v_one_half, ctps_mid_height_top],
                        [-sep_distance_aux, v_one_half, ctps_mid_height_top],
                        [-aux_column_width, v_one_half, ctps_mid_height_top],
                        [-v_one_half, aux_column_width, v_one],
                        [-sep_distance_aux, aux_column_width, v_one],
                        [-branch_thickness, branch_thickness, v_one],
                        [-v_one_half, sep_distance_aux, v_one],
                        [-sep_distance_aux, sep_distance_aux, v_one],
                        [-aux_column_width, sep_distance_aux, v_one],
                        [-v_one_half, v_one_half, v_one],
                        [-sep_distance_aux, v_one_half, v_one],
                        [-aux_column_width, v_one_half, v_one],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[2, 2, 2],
                        control_points=branch_x_min_y_max_ctps,
                    )
                )

                branch_x_max_y_min_ctps = _np.array(
                    [
                        [r_center, -v_one_half, fill_height_aux],
                        [half_r_center, -v_one_half, fill_height_aux],
                        [v_one_half, -v_one_half, fill_height_aux],
                        [r_center, -half_r_center, fill_height_aux],
                        [half_r_center, -half_r_center, fill_height_aux],
                        [v_one_half, -half_r_center, fill_height_aux],
                        [r_center, -r_center, fill_height_aux],
                        [half_r_center, -r_center, fill_height_aux],
                        [v_one_half, -r_center, fill_height_aux],
                        [aux_column_width, -v_one_half, ctps_mid_height_top],
                        [sep_distance_aux, -v_one_half, ctps_mid_height_top],
                        [v_one_half, -v_one_half, ctps_mid_height_top],
                        [
                            aux_column_width,
                            -sep_distance_aux,
                            ctps_mid_height_top,
                        ],
                        [
                            sep_distance_aux,
                            -sep_distance_aux,
                            ctps_mid_height_top,
                        ],
                        [v_one_half, -sep_distance_aux, ctps_mid_height_top],
                        [
                            branch_thickness,
                            -branch_thickness,
                            ctps_mid_height_top,
                        ],
                        [
                            sep_distance_aux,
                            -aux_column_width,
                            ctps_mid_height_top,
                        ],
                        [v_one_half, -aux_column_width, ctps_mid_height_top],
                        [aux_column_width, -v_one_half, v_one],
                        [sep_distance_aux, -v_one_half, v_one],
                        [v_one_half, -v_one_half, v_one],
                        [aux_column_width, -sep_distance_aux, v_one],
                        [sep_distance_aux, -sep_distance_aux, v_one],
                        [v_one_half, -sep_distance_aux, v_one],
                        [branch_thickness, -branch_thickness, v_one],
                        [sep_distance_aux, -aux_column_width, v_one],
                        [v_one_half, -aux_column_width, v_one],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[2, 2, 2],
                        control_points=branch_x_max_y_min_ctps,
                    )
                )

                branch_x_max_y_max_ctps = _np.array(
                    [
                        [r_center, r_center, fill_height_aux],
                        [half_r_center, r_center, fill_height_aux],
                        [v_one_half, r_center, fill_height_aux],
                        [r_center, half_r_center, fill_height_aux],
                        [half_r_center, half_r_center, fill_height_aux],
                        [v_one_half, half_r_center, fill_height_aux],
                        [r_center, v_one_half, fill_height_aux],
                        [half_r_center, v_one_half, fill_height_aux],
                        [v_one_half, v_one_half, fill_height_aux],
                        [
                            branch_thickness,
                            branch_thickness,
                            ctps_mid_height_top,
                        ],
                        [
                            sep_distance_aux,
                            aux_column_width,
                            ctps_mid_height_top,
                        ],
                        [v_one_half, aux_column_width, ctps_mid_height_top],
                        [
                            aux_column_width,
                            sep_distance_aux,
                            ctps_mid_height_top,
                        ],
                        [
                            sep_distance_aux,
                            sep_distance_aux,
                            ctps_mid_height_top,
                        ],
                        [v_one_half, sep_distance_aux, ctps_mid_height_top],
                        [aux_column_width, v_one_half, ctps_mid_height_top],
                        [sep_distance_aux, v_one_half, ctps_mid_height_top],
                        [v_one_half, v_one_half, ctps_mid_height_top],
                        [branch_thickness, branch_thickness, v_one],
                        [sep_distance_aux, aux_column_width, v_one],
                        [v_one_half, aux_column_width, v_one],
                        [aux_column_width, sep_distance_aux, v_one],
                        [sep_distance_aux, sep_distance_aux, v_one],
                        [v_one_half, sep_distance_aux, v_one],
                        [aux_column_width, v_one_half, v_one],
                        [sep_distance_aux, v_one_half, v_one],
                        [v_one_half, v_one_half, v_one],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[2, 2, 2],
                        control_points=branch_x_max_y_max_ctps,
                    )
                )

            elif closure == "z_max":
                branch_neighbor_x_min_ctps = _np.array(
                    [
                        [-v_one_half, -aux_column_width, v_zero],
                        [-sep_distance_aux, -aux_column_width, v_zero],
                        [-branch_thickness, -branch_thickness, v_zero],
                        [-v_one_half, aux_column_width, v_zero],
                        [-sep_distance_aux, aux_column_width, v_zero],
                        [-branch_thickness, branch_thickness, v_zero],
                        [
                            -v_one_half,
                            -aux_column_width,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -sep_distance_aux,
                            -aux_column_width,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -branch_thickness,
                            -branch_thickness,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -v_one_half,
                            aux_column_width,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -sep_distance_aux,
                            aux_column_width,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -branch_thickness,
                            branch_thickness,
                            ctps_mid_height_bottom,
                        ],
                        [-v_one_half, -r_center, inv_filling_height],
                        [-half_r_center, -r_center, inv_filling_height],
                        [-r_center, -r_center, inv_filling_height],
                        [-v_one_half, r_center, inv_filling_height],
                        [-half_r_center, r_center, inv_filling_height],
                        [-r_center, r_center, inv_filling_height],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[2, 1, 2],
                        control_points=branch_neighbor_x_min_ctps,
                    )
                )

                branch_neighbor_x_max_ctps = _np.array(
                    [
                        [branch_thickness, -branch_thickness, v_zero],
                        [sep_distance_aux, -aux_column_width, v_zero],
                        [v_one_half, -aux_column_width, v_zero],
                        [branch_thickness, branch_thickness, v_zero],
                        [sep_distance_aux, aux_column_width, v_zero],
                        [v_one_half, aux_column_width, v_zero],
                        [
                            branch_thickness,
                            -branch_thickness,
                            ctps_mid_height_bottom,
                        ],
                        [
                            sep_distance_aux,
                            -aux_column_width,
                            ctps_mid_height_bottom,
                        ],
                        [
                            v_one_half,
                            -aux_column_width,
                            ctps_mid_height_bottom,
                        ],
                        [
                            branch_thickness,
                            branch_thickness,
                            ctps_mid_height_bottom,
                        ],
                        [
                            sep_distance_aux,
                            aux_column_width,
                            ctps_mid_height_bottom,
                        ],
                        [v_one_half, aux_column_width, ctps_mid_height_bottom],
                        [r_center, -r_center, inv_filling_height],
                        [half_r_center, -r_center, inv_filling_height],
                        [v_one_half, -r_center, inv_filling_height],
                        [r_center, r_center, inv_filling_height],
                        [half_r_center, r_center, inv_filling_height],
                        [v_one_half, r_center, inv_filling_height],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[2, 1, 2],
                        control_points=branch_neighbor_x_max_ctps,
                    )
                )

                branch_neighbor_y_min_ctps = _np.array(
                    [
                        [-aux_column_width, -v_one_half, v_zero],
                        [aux_column_width, -v_one_half, v_zero],
                        [-aux_column_width, -sep_distance_aux, v_zero],
                        [aux_column_width, -sep_distance_aux, v_zero],
                        [-branch_thickness, -branch_thickness, v_zero],
                        [branch_thickness, -branch_thickness, v_zero],
                        [
                            -aux_column_width,
                            -v_one_half,
                            ctps_mid_height_bottom,
                        ],
                        [
                            aux_column_width,
                            -v_one_half,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -aux_column_width,
                            -sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
                        [
                            aux_column_width,
                            -sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
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
                        [-r_center, -v_one_half, inv_filling_height],
                        [r_center, -v_one_half, inv_filling_height],
                        [-r_center, -half_r_center, inv_filling_height],
                        [r_center, -half_r_center, inv_filling_height],
                        [-r_center, -r_center, inv_filling_height],
                        [r_center, -r_center, inv_filling_height],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[1, 2, 2],
                        control_points=branch_neighbor_y_min_ctps,
                    )
                )

                branch_neighbor_y_max_ctps = _np.array(
                    [
                        [-branch_thickness, branch_thickness, v_zero],
                        [branch_thickness, branch_thickness, v_zero],
                        [-aux_column_width, sep_distance_aux, v_zero],
                        [aux_column_width, sep_distance_aux, v_zero],
                        [-aux_column_width, v_one_half, v_zero],
                        [aux_column_width, v_one_half, v_zero],
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
                        [
                            -aux_column_width,
                            sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
                        [
                            aux_column_width,
                            sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -aux_column_width,
                            v_one_half,
                            ctps_mid_height_bottom,
                        ],
                        [aux_column_width, v_one_half, ctps_mid_height_bottom],
                        [-r_center, r_center, inv_filling_height],
                        [r_center, r_center, inv_filling_height],
                        [-r_center, half_r_center, inv_filling_height],
                        [r_center, half_r_center, inv_filling_height],
                        [-r_center, v_one_half, inv_filling_height],
                        [r_center, v_one_half, inv_filling_height],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[1, 2, 2],
                        control_points=branch_neighbor_y_max_ctps,
                    )
                )

                branch_x_min_y_min_ctps = _np.array(
                    [
                        [-v_one_half, -v_one_half, v_zero],
                        [-sep_distance_aux, -v_one_half, v_zero],
                        [-aux_column_width, -v_one_half, v_zero],
                        [-v_one_half, -sep_distance_aux, v_zero],
                        [-sep_distance_aux, -sep_distance_aux, v_zero],
                        [-aux_column_width, -sep_distance_aux, v_zero],
                        [-v_one_half, -aux_column_width, v_zero],
                        [-sep_distance_aux, -aux_column_width, v_zero],
                        [-branch_thickness, -branch_thickness, v_zero],
                        [-v_one_half, -v_one_half, ctps_mid_height_bottom],
                        [
                            -sep_distance_aux,
                            -v_one_half,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -aux_column_width,
                            -v_one_half,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -v_one_half,
                            -sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -sep_distance_aux,
                            -sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -aux_column_width,
                            -sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -v_one_half,
                            -aux_column_width,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -sep_distance_aux,
                            -aux_column_width,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -branch_thickness,
                            -branch_thickness,
                            ctps_mid_height_bottom,
                        ],
                        [-v_one_half, -v_one_half, inv_filling_height],
                        [-half_r_center, -v_one_half, inv_filling_height],
                        [-r_center, -v_one_half, inv_filling_height],
                        [-v_one_half, -half_r_center, inv_filling_height],
                        [-half_r_center, -half_r_center, inv_filling_height],
                        [-r_center, -half_r_center, inv_filling_height],
                        [-v_one_half, -r_center, inv_filling_height],
                        [-half_r_center, -r_center, inv_filling_height],
                        [-r_center, -r_center, inv_filling_height],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[2, 2, 2],
                        control_points=branch_x_min_y_min_ctps,
                    )
                )

                branch_x_max_y_max_ctps = _np.array(
                    [
                        [branch_thickness, branch_thickness, v_zero],
                        [sep_distance_aux, aux_column_width, v_zero],
                        [v_one_half, aux_column_width, v_zero],
                        [aux_column_width, sep_distance_aux, v_zero],
                        [sep_distance_aux, sep_distance_aux, v_zero],
                        [v_one_half, sep_distance_aux, v_zero],
                        [aux_column_width, v_one_half, v_zero],
                        [sep_distance_aux, v_one_half, v_zero],
                        [v_one_half, v_one_half, v_zero],
                        [
                            branch_thickness,
                            branch_thickness,
                            ctps_mid_height_bottom,
                        ],
                        [
                            sep_distance_aux,
                            aux_column_width,
                            ctps_mid_height_bottom,
                        ],
                        [v_one_half, aux_column_width, ctps_mid_height_bottom],
                        [
                            aux_column_width,
                            sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
                        [
                            sep_distance_aux,
                            sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
                        [v_one_half, sep_distance_aux, ctps_mid_height_bottom],
                        [aux_column_width, v_one_half, ctps_mid_height_bottom],
                        [sep_distance_aux, v_one_half, ctps_mid_height_bottom],
                        [v_one_half, v_one_half, ctps_mid_height_bottom],
                        [r_center, r_center, inv_filling_height],
                        [half_r_center, r_center, inv_filling_height],
                        [v_one_half, r_center, inv_filling_height],
                        [r_center, half_r_center, inv_filling_height],
                        [half_r_center, half_r_center, inv_filling_height],
                        [v_one_half, half_r_center, inv_filling_height],
                        [r_center, v_one_half, inv_filling_height],
                        [half_r_center, v_one_half, inv_filling_height],
                        [v_one_half, v_one_half, inv_filling_height],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[2, 2, 2],
                        control_points=branch_x_max_y_max_ctps,
                    )
                )

                branch_x_max_y_min_ctps = _np.array(
                    [
                        [aux_column_width, -v_one_half, v_zero],
                        [sep_distance_aux, -v_one_half, v_zero],
                        [v_one_half, -v_one_half, v_zero],
                        [aux_column_width, -sep_distance_aux, v_zero],
                        [sep_distance_aux, -sep_distance_aux, v_zero],
                        [v_one_half, -sep_distance_aux, v_zero],
                        [branch_thickness, -branch_thickness, v_zero],
                        [sep_distance_aux, -aux_column_width, v_zero],
                        [v_one_half, -aux_column_width, v_zero],
                        [
                            aux_column_width,
                            -v_one_half,
                            ctps_mid_height_bottom,
                        ],
                        [
                            sep_distance_aux,
                            -v_one_half,
                            ctps_mid_height_bottom,
                        ],
                        [v_one_half, -v_one_half, ctps_mid_height_bottom],
                        [
                            aux_column_width,
                            -sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
                        [
                            sep_distance_aux,
                            -sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
                        [
                            v_one_half,
                            -sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
                        [
                            branch_thickness,
                            -branch_thickness,
                            ctps_mid_height_bottom,
                        ],
                        [
                            sep_distance_aux,
                            -aux_column_width,
                            ctps_mid_height_bottom,
                        ],
                        [
                            v_one_half,
                            -aux_column_width,
                            ctps_mid_height_bottom,
                        ],
                        [r_center, -v_one_half, inv_filling_height],
                        [half_r_center, -v_one_half, inv_filling_height],
                        [v_one_half, -v_one_half, inv_filling_height],
                        [r_center, -half_r_center, inv_filling_height],
                        [half_r_center, -half_r_center, inv_filling_height],
                        [v_one_half, -half_r_center, inv_filling_height],
                        [r_center, -r_center, inv_filling_height],
                        [half_r_center, -r_center, inv_filling_height],
                        [v_one_half, -r_center, inv_filling_height],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[2, 2, 2],
                        control_points=branch_x_max_y_min_ctps,
                    )
                )

                branch_x_min_y_max_ctps = _np.array(
                    [
                        [-v_one_half, aux_column_width, v_zero],
                        [-sep_distance_aux, aux_column_width, v_zero],
                        [-branch_thickness, branch_thickness, v_zero],
                        [-v_one_half, sep_distance_aux, v_zero],
                        [-sep_distance_aux, sep_distance_aux, v_zero],
                        [-aux_column_width, sep_distance_aux, v_zero],
                        [-v_one_half, v_one_half, v_zero],
                        [-sep_distance_aux, v_one_half, v_zero],
                        [-aux_column_width, v_one_half, v_zero],
                        [
                            -v_one_half,
                            aux_column_width,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -sep_distance_aux,
                            aux_column_width,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -branch_thickness,
                            branch_thickness,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -v_one_half,
                            sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -sep_distance_aux,
                            sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -aux_column_width,
                            sep_distance_aux,
                            ctps_mid_height_bottom,
                        ],
                        [-v_one_half, v_one_half, ctps_mid_height_bottom],
                        [
                            -sep_distance_aux,
                            v_one_half,
                            ctps_mid_height_bottom,
                        ],
                        [
                            -aux_column_width,
                            v_one_half,
                            ctps_mid_height_bottom,
                        ],
                        [-v_one_half, r_center, inv_filling_height],
                        [-half_r_center, r_center, inv_filling_height],
                        [-r_center, r_center, inv_filling_height],
                        [-v_one_half, half_r_center, inv_filling_height],
                        [-half_r_center, half_r_center, inv_filling_height],
                        [-r_center, half_r_center, inv_filling_height],
                        [-v_one_half, v_one_half, inv_filling_height],
                        [-half_r_center, v_one_half, inv_filling_height],
                        [-r_center, v_one_half, inv_filling_height],
                    ]
                ) + _np.array([v_one_half, v_one_half, v_zero])

                spline_list.append(
                    _Bezier(
                        degrees=[2, 2, 2],
                        control_points=branch_x_min_y_max_ctps,
                    )
                )
            else:
                raise ValueError("Corner Type not supported")

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        return (splines, derivatives)

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        separator_distance=0.3,
        center_expansion=1.0,
        closure=None,
        **kwargs,  # noqa ARG002
    ):
        """Create an inverse microtile based on the parameters that describe
        the branch thicknesses.

        Thickness parameters are used to describe the inner radius of the
        outward facing branches

        Parameters
        ----------
        parameters : np.array
          only first entry is used, defines the internal radii of the
          branches
        separator_distance : float
          Control point distance to separation layer of biquadratic degrees,
          determines the minimum branch thickness (defaults to 0.3)
        center_expansion : float
          thickness of center is expanded by a factor (default to 1.0), which
          determines the maximum branch thickness
        closure : str
          parametric dimension that needs to be closed.
          Must be one of {"z_min", "z_max"}

        Returns
        -------
        microtile_list : list(splines)
        """

        self._check_custom_parameter(
            center_expansion, "center expansion", self._CENTER_EXPANSION_BOUNDS
        )

        # Check if all radii are in allowed range
        max_radius = min(0.5, (0.5 / center_expansion))
        max_radius = min(max_radius, separator_distance)
        min_radius = max(0.5 - separator_distance, 0)

        parameters, n_derivatives, derivatives = self._process_input(
            parameters=parameters,
            parameter_sensitivities=parameter_sensitivities,
        )

        if _np.any(parameters < min_radius) or _np.any(
            parameters > max_radius
        ):
            raise ValueError(
                f"Radii must be in ({min_radius},{max_radius}) for "
                f"center_expansion {center_expansion}"
            )

        if closure is not None:
            return self._closing_tile(
                parameters=parameters,
                parameter_sensitivities=parameter_sensitivities,
                separator_distance=separator_distance,
                closure=closure,
                **kwargs,
            )

        splines = []
        for i_derivative in range(n_derivatives + 1):
            if i_derivative == 0:
                [
                    x_min_r,
                    x_max_r,
                    y_min_r,
                    y_max_r,
                    z_min_r,
                    z_max_r,
                ] = parameters.flatten()

                # center radius
                center_r = _np.sum(parameters) / 6.0 * center_expansion

                # Auxiliary values for smooothing (mid-branch thickness)
                [
                    aux_x_min,
                    aux_x_max,
                    aux_y_min,
                    aux_y_max,
                    aux_z_min,
                    aux_z_max,
                ] = _np.minimum(parameters.ravel(), center_r)
                # Branch midlength
                hd_center = 0.5 * (0.5 + center_r)
                aux_column_width = 0.5 - 2 * (0.5 - separator_distance)

                v_one_half = 0.5
                center_point = _np.array([0.5, 0.5, 0.5])

                sep_distance = separator_distance
            else:
                sensitivities_i = parameter_sensitivities[
                    :, 0, i_derivative - 1
                ]
                [
                    x_min_r,
                    x_max_r,
                    y_min_r,
                    y_max_r,
                    z_min_r,
                    z_max_r,
                ] = sensitivities_i

                center_r = _np.mean(sensitivities_i) * center_expansion
                [
                    aux_x_min,
                    aux_x_max,
                    aux_y_min,
                    aux_y_max,
                    aux_z_min,
                    aux_z_max,
                ] = _np.where(
                    parameters.ravel()
                    < _np.mean(parameters.ravel()) * center_expansion,
                    sensitivities_i,
                    center_r,
                )
                # Branch midlength
                hd_center = center_r / 2.0
                aux_column_width = 0.0
                v_one_half = 0.0
                center_point = _np.zeros(3)

                sep_distance = 0.0

            # Init return type
            spline_list = []

            # Start with branch interconnections
            x_min_y_min = (
                _np.array(
                    [
                        [-v_one_half, -v_one_half, -aux_column_width],
                        [-sep_distance, -v_one_half, -aux_column_width],
                        [-y_min_r, -v_one_half, -y_min_r],
                        [-v_one_half, -sep_distance, -aux_column_width],
                        [-hd_center, -hd_center, -aux_column_width],
                        [-aux_y_min, -hd_center, -aux_y_min],
                        [-v_one_half, -x_min_r, -x_min_r],
                        [-hd_center, -aux_x_min, -aux_x_min],
                        [-center_r, -center_r, -center_r],
                        [-v_one_half, -v_one_half, aux_column_width],
                        [-sep_distance, -v_one_half, aux_column_width],
                        [-y_min_r, -v_one_half, y_min_r],
                        [-v_one_half, -sep_distance, aux_column_width],
                        [-hd_center, -hd_center, aux_column_width],
                        [-aux_y_min, -hd_center, aux_y_min],
                        [-v_one_half, -x_min_r, x_min_r],
                        [-hd_center, -aux_x_min, aux_x_min],
                        [-center_r, -center_r, center_r],
                    ]
                )
                + center_point
            )

            x_max_y_min = (
                _np.array(
                    [
                        [y_min_r, -v_one_half, -y_min_r],
                        [sep_distance, -v_one_half, -aux_column_width],
                        [v_one_half, -v_one_half, -aux_column_width],
                        [aux_y_min, -hd_center, -aux_y_min],
                        [hd_center, -hd_center, -aux_column_width],
                        [v_one_half, -sep_distance, -aux_column_width],
                        [center_r, -center_r, -center_r],
                        [hd_center, -aux_x_max, -aux_x_max],
                        [v_one_half, -x_max_r, -x_max_r],
                        [y_min_r, -v_one_half, y_min_r],
                        [sep_distance, -v_one_half, aux_column_width],
                        [v_one_half, -v_one_half, aux_column_width],
                        [aux_y_min, -hd_center, aux_y_min],
                        [hd_center, -hd_center, aux_column_width],
                        [v_one_half, -sep_distance, aux_column_width],
                        [center_r, -center_r, center_r],
                        [hd_center, -aux_x_max, aux_x_max],
                        [v_one_half, -x_max_r, x_max_r],
                    ]
                )
                + center_point
            )

            x_min_y_max = (
                _np.array(
                    [
                        [-v_one_half, x_min_r, -x_min_r],
                        [-hd_center, aux_x_min, -aux_x_min],
                        [-center_r, center_r, -center_r],
                        [-v_one_half, sep_distance, -aux_column_width],
                        [-hd_center, hd_center, -aux_column_width],
                        [-aux_y_max, hd_center, -aux_y_max],
                        [-v_one_half, v_one_half, -aux_column_width],
                        [-sep_distance, v_one_half, -aux_column_width],
                        [-y_max_r, v_one_half, -y_max_r],
                        [-v_one_half, x_min_r, x_min_r],
                        [-hd_center, aux_x_min, aux_x_min],
                        [-center_r, center_r, center_r],
                        [-v_one_half, sep_distance, aux_column_width],
                        [-hd_center, hd_center, aux_column_width],
                        [-aux_y_max, hd_center, aux_y_max],
                        [-v_one_half, v_one_half, aux_column_width],
                        [-sep_distance, v_one_half, aux_column_width],
                        [-y_max_r, v_one_half, y_max_r],
                    ]
                )
                + center_point
            )

            x_max_y_max = (
                _np.array(
                    [
                        [center_r, center_r, -center_r],
                        [hd_center, aux_x_max, -aux_x_max],
                        [v_one_half, x_max_r, -x_max_r],
                        [aux_y_max, hd_center, -aux_y_max],
                        [hd_center, hd_center, -aux_column_width],
                        [v_one_half, sep_distance, -aux_column_width],
                        [y_max_r, v_one_half, -y_max_r],
                        [sep_distance, v_one_half, -aux_column_width],
                        [v_one_half, v_one_half, -aux_column_width],
                        [center_r, center_r, center_r],
                        [hd_center, aux_x_max, aux_x_max],
                        [v_one_half, x_max_r, x_max_r],
                        [aux_y_max, hd_center, aux_y_max],
                        [hd_center, hd_center, aux_column_width],
                        [v_one_half, sep_distance, aux_column_width],
                        [y_max_r, v_one_half, y_max_r],
                        [sep_distance, v_one_half, aux_column_width],
                        [v_one_half, v_one_half, aux_column_width],
                    ]
                )
                + center_point
            )

            x_min_z_min = (
                _np.array(
                    [
                        [-v_one_half, -aux_column_width, -v_one_half],
                        [-sep_distance, -aux_column_width, -v_one_half],
                        [-z_min_r, -z_min_r, -v_one_half],
                        [-v_one_half, aux_column_width, -v_one_half],
                        [-sep_distance, aux_column_width, -v_one_half],
                        [-z_min_r, z_min_r, -v_one_half],
                        [-v_one_half, -aux_column_width, -sep_distance],
                        [-hd_center, -aux_column_width, -hd_center],
                        [-aux_z_min, -aux_z_min, -hd_center],
                        [-v_one_half, aux_column_width, -sep_distance],
                        [-hd_center, aux_column_width, -hd_center],
                        [-aux_z_min, aux_z_min, -hd_center],
                        [-v_one_half, -x_min_r, -x_min_r],
                        [-hd_center, -aux_x_min, -aux_x_min],
                        [-center_r, -center_r, -center_r],
                        [-v_one_half, x_min_r, -x_min_r],
                        [-hd_center, aux_x_min, -aux_x_min],
                        [-center_r, center_r, -center_r],
                    ]
                )
                + center_point
            )

            x_max_z_min = (
                _np.array(
                    [
                        [z_min_r, -z_min_r, -v_one_half],
                        [sep_distance, -aux_column_width, -v_one_half],
                        [v_one_half, -aux_column_width, -v_one_half],
                        [z_min_r, z_min_r, -v_one_half],
                        [sep_distance, aux_column_width, -v_one_half],
                        [v_one_half, aux_column_width, -v_one_half],
                        [aux_z_min, -aux_z_min, -hd_center],
                        [hd_center, -aux_column_width, -hd_center],
                        [v_one_half, -aux_column_width, -sep_distance],
                        [aux_z_min, aux_z_min, -hd_center],
                        [hd_center, aux_column_width, -hd_center],
                        [v_one_half, aux_column_width, -sep_distance],
                        [center_r, -center_r, -center_r],
                        [hd_center, -aux_x_max, -aux_x_max],
                        [v_one_half, -x_max_r, -x_max_r],
                        [center_r, center_r, -center_r],
                        [hd_center, aux_x_max, -aux_x_max],
                        [v_one_half, x_max_r, -x_max_r],
                    ]
                )
                + center_point
            )

            x_min_z_max = (
                _np.array(
                    [
                        [-v_one_half, -x_min_r, x_min_r],
                        [-hd_center, -aux_x_min, aux_x_min],
                        [-center_r, -center_r, center_r],
                        [-v_one_half, x_min_r, x_min_r],
                        [-hd_center, aux_x_min, aux_x_min],
                        [-center_r, center_r, center_r],
                        [-v_one_half, -aux_column_width, sep_distance],
                        [-hd_center, -aux_column_width, hd_center],
                        [-aux_z_max, -aux_z_max, hd_center],
                        [-v_one_half, aux_column_width, sep_distance],
                        [-hd_center, aux_column_width, hd_center],
                        [-aux_z_max, aux_z_max, hd_center],
                        [-v_one_half, -aux_column_width, v_one_half],
                        [-sep_distance, -aux_column_width, v_one_half],
                        [-z_max_r, -z_max_r, v_one_half],
                        [-v_one_half, aux_column_width, v_one_half],
                        [-sep_distance, aux_column_width, v_one_half],
                        [-z_max_r, z_max_r, v_one_half],
                    ]
                )
                + center_point
            )

            x_max_z_max = (
                _np.array(
                    [
                        [center_r, -center_r, center_r],
                        [hd_center, -aux_x_max, aux_x_max],
                        [v_one_half, -x_max_r, x_max_r],
                        [center_r, center_r, center_r],
                        [hd_center, aux_x_max, aux_x_max],
                        [v_one_half, x_max_r, x_max_r],
                        [aux_z_max, -aux_z_max, hd_center],
                        [hd_center, -aux_column_width, hd_center],
                        [v_one_half, -aux_column_width, sep_distance],
                        [aux_z_max, aux_z_max, hd_center],
                        [hd_center, aux_column_width, hd_center],
                        [v_one_half, aux_column_width, sep_distance],
                        [z_max_r, -z_max_r, v_one_half],
                        [sep_distance, -aux_column_width, v_one_half],
                        [v_one_half, -aux_column_width, v_one_half],
                        [z_max_r, z_max_r, v_one_half],
                        [sep_distance, aux_column_width, v_one_half],
                        [v_one_half, aux_column_width, v_one_half],
                    ]
                )
                + center_point
            )

            y_min_z_min = (
                _np.array(
                    [
                        [-aux_column_width, -v_one_half, -v_one_half],
                        [aux_column_width, -v_one_half, -v_one_half],
                        [-aux_column_width, -sep_distance, -v_one_half],
                        [aux_column_width, -sep_distance, -v_one_half],
                        [-z_min_r, -z_min_r, -v_one_half],
                        [z_min_r, -z_min_r, -v_one_half],
                        [-aux_column_width, -v_one_half, -sep_distance],
                        [aux_column_width, -v_one_half, -sep_distance],
                        [-aux_column_width, -hd_center, -hd_center],
                        [aux_column_width, -hd_center, -hd_center],
                        [-aux_z_min, -aux_z_min, -hd_center],
                        [aux_z_min, -aux_z_min, -hd_center],
                        [-y_min_r, -v_one_half, -y_min_r],
                        [y_min_r, -v_one_half, -y_min_r],
                        [-aux_y_min, -hd_center, -aux_y_min],
                        [aux_y_min, -hd_center, -aux_y_min],
                        [-center_r, -center_r, -center_r],
                        [center_r, -center_r, -center_r],
                    ]
                )
                + center_point
            )

            y_max_z_min = (
                _np.array(
                    [
                        [-z_min_r, z_min_r, -v_one_half],
                        [z_min_r, z_min_r, -v_one_half],
                        [-aux_column_width, sep_distance, -v_one_half],
                        [aux_column_width, sep_distance, -v_one_half],
                        [-aux_column_width, v_one_half, -v_one_half],
                        [aux_column_width, v_one_half, -v_one_half],
                        [-aux_z_min, aux_z_min, -hd_center],
                        [aux_z_min, aux_z_min, -hd_center],
                        [-aux_column_width, hd_center, -hd_center],
                        [aux_column_width, hd_center, -hd_center],
                        [-aux_column_width, v_one_half, -sep_distance],
                        [aux_column_width, v_one_half, -sep_distance],
                        [-center_r, center_r, -center_r],
                        [center_r, center_r, -center_r],
                        [-aux_y_max, hd_center, -aux_y_max],
                        [aux_y_max, hd_center, -aux_y_max],
                        [-y_max_r, v_one_half, -y_max_r],
                        [y_max_r, v_one_half, -y_max_r],
                    ]
                )
                + center_point
            )

            y_min_z_max = (
                _np.array(
                    [
                        [-y_min_r, -v_one_half, y_min_r],
                        [y_min_r, -v_one_half, y_min_r],
                        [-aux_y_min, -hd_center, aux_y_min],
                        [aux_y_min, -hd_center, aux_y_min],
                        [-center_r, -center_r, center_r],
                        [center_r, -center_r, center_r],
                        [-aux_column_width, -v_one_half, sep_distance],
                        [aux_column_width, -v_one_half, sep_distance],
                        [-aux_column_width, -hd_center, hd_center],
                        [aux_column_width, -hd_center, hd_center],
                        [-aux_z_max, -aux_z_max, hd_center],
                        [aux_z_max, -aux_z_max, hd_center],
                        [-aux_column_width, -v_one_half, v_one_half],
                        [aux_column_width, -v_one_half, v_one_half],
                        [-aux_column_width, -sep_distance, v_one_half],
                        [aux_column_width, -sep_distance, v_one_half],
                        [-z_max_r, -z_max_r, v_one_half],
                        [z_max_r, -z_max_r, v_one_half],
                    ]
                )
                + center_point
            )

            y_max_z_max = (
                _np.array(
                    [
                        [-center_r, center_r, center_r],
                        [center_r, center_r, center_r],
                        [-aux_y_max, hd_center, aux_y_max],
                        [aux_y_max, hd_center, aux_y_max],
                        [-y_max_r, v_one_half, y_max_r],
                        [y_max_r, v_one_half, y_max_r],
                        [-aux_z_max, aux_z_max, hd_center],
                        [aux_z_max, aux_z_max, hd_center],
                        [-aux_column_width, hd_center, hd_center],
                        [aux_column_width, hd_center, hd_center],
                        [-aux_column_width, v_one_half, sep_distance],
                        [aux_column_width, v_one_half, sep_distance],
                        [-z_max_r, z_max_r, v_one_half],
                        [z_max_r, z_max_r, v_one_half],
                        [-aux_column_width, sep_distance, v_one_half],
                        [aux_column_width, sep_distance, v_one_half],
                        [-aux_column_width, v_one_half, v_one_half],
                        [aux_column_width, v_one_half, v_one_half],
                    ]
                )
                + center_point
            )

            x_min_y_min_z_min = (
                _np.array(
                    [
                        [-v_one_half, -v_one_half, -v_one_half],
                        [-sep_distance, -v_one_half, -v_one_half],
                        [-aux_column_width, -v_one_half, -v_one_half],
                        [-v_one_half, -sep_distance, -v_one_half],
                        [
                            -sep_distance,
                            -sep_distance,
                            -v_one_half,
                        ],
                        [-aux_column_width, -sep_distance, -v_one_half],
                        [-v_one_half, -aux_column_width, -v_one_half],
                        [-sep_distance, -aux_column_width, -v_one_half],
                        [-z_min_r, -z_min_r, -v_one_half],
                        [-v_one_half, -v_one_half, -sep_distance],
                        [
                            -sep_distance,
                            -v_one_half,
                            -sep_distance,
                        ],
                        [-aux_column_width, -v_one_half, -sep_distance],
                        [
                            -v_one_half,
                            -sep_distance,
                            -sep_distance,
                        ],
                        [-hd_center, -hd_center, -hd_center],
                        [-aux_column_width, -hd_center, -hd_center],
                        [-v_one_half, -aux_column_width, -sep_distance],
                        [-hd_center, -aux_column_width, -hd_center],
                        [-aux_z_min, -aux_z_min, -hd_center],
                        [-v_one_half, -v_one_half, -aux_column_width],
                        [-sep_distance, -v_one_half, -aux_column_width],
                        [-y_min_r, -v_one_half, -y_min_r],
                        [-v_one_half, -sep_distance, -aux_column_width],
                        [-hd_center, -hd_center, -aux_column_width],
                        [-aux_y_min, -hd_center, -aux_y_min],
                        [-v_one_half, -x_min_r, -x_min_r],
                        [-hd_center, -aux_x_min, -aux_x_min],
                        [-center_r, -center_r, -center_r],
                    ]
                )
                + center_point
            )

            x_max_y_min_z_min = (
                _np.array(
                    [
                        [aux_column_width, -v_one_half, -v_one_half],
                        [sep_distance, -v_one_half, -v_one_half],
                        [v_one_half, -v_one_half, -v_one_half],
                        [aux_column_width, -sep_distance, -v_one_half],
                        [sep_distance, -sep_distance, -v_one_half],
                        [v_one_half, -sep_distance, -v_one_half],
                        [z_min_r, -z_min_r, -v_one_half],
                        [sep_distance, -aux_column_width, -v_one_half],
                        [v_one_half, -aux_column_width, -v_one_half],
                        [aux_column_width, -v_one_half, -sep_distance],
                        [sep_distance, -v_one_half, -sep_distance],
                        [v_one_half, -v_one_half, -sep_distance],
                        [aux_column_width, -hd_center, -hd_center],
                        [hd_center, -hd_center, -hd_center],
                        [v_one_half, -sep_distance, -sep_distance],
                        [aux_z_min, -aux_z_min, -hd_center],
                        [hd_center, -aux_column_width, -hd_center],
                        [v_one_half, -aux_column_width, -sep_distance],
                        [y_min_r, -v_one_half, -y_min_r],
                        [sep_distance, -v_one_half, -aux_column_width],
                        [v_one_half, -v_one_half, -aux_column_width],
                        [aux_y_min, -hd_center, -aux_y_min],
                        [hd_center, -hd_center, -aux_column_width],
                        [v_one_half, -sep_distance, -aux_column_width],
                        [center_r, -center_r, -center_r],
                        [hd_center, -aux_x_max, -aux_x_max],
                        [v_one_half, -x_max_r, -x_max_r],
                    ]
                )
                + center_point
            )

            x_min_y_max_z_min = (
                _np.array(
                    [
                        [-v_one_half, aux_column_width, -v_one_half],
                        [-sep_distance, aux_column_width, -v_one_half],
                        [-z_min_r, z_min_r, -v_one_half],
                        [-v_one_half, sep_distance, -v_one_half],
                        [-sep_distance, sep_distance, -v_one_half],
                        [-aux_column_width, sep_distance, -v_one_half],
                        [-v_one_half, v_one_half, -v_one_half],
                        [-sep_distance, v_one_half, -v_one_half],
                        [-aux_column_width, v_one_half, -v_one_half],
                        [-v_one_half, aux_column_width, -sep_distance],
                        [-hd_center, aux_column_width, -hd_center],
                        [-aux_z_min, aux_z_min, -hd_center],
                        [-v_one_half, sep_distance, -sep_distance],
                        [-hd_center, hd_center, -hd_center],
                        [-aux_column_width, hd_center, -hd_center],
                        [-v_one_half, v_one_half, -sep_distance],
                        [-sep_distance, v_one_half, -sep_distance],
                        [-aux_column_width, v_one_half, -sep_distance],
                        [-v_one_half, x_min_r, -x_min_r],
                        [-hd_center, aux_x_min, -aux_x_min],
                        [-center_r, center_r, -center_r],
                        [-v_one_half, sep_distance, -aux_column_width],
                        [-hd_center, hd_center, -aux_column_width],
                        [-aux_y_max, hd_center, -aux_y_max],
                        [-v_one_half, v_one_half, -aux_column_width],
                        [-sep_distance, v_one_half, -aux_column_width],
                        [-y_max_r, v_one_half, -y_max_r],
                    ]
                )
                + center_point
            )

            x_max_y_max_z_min = (
                _np.array(
                    [
                        [z_min_r, z_min_r, -v_one_half],
                        [sep_distance, aux_column_width, -v_one_half],
                        [v_one_half, aux_column_width, -v_one_half],
                        [aux_column_width, sep_distance, -v_one_half],
                        [sep_distance, sep_distance, -v_one_half],
                        [v_one_half, sep_distance, -v_one_half],
                        [aux_column_width, v_one_half, -v_one_half],
                        [sep_distance, v_one_half, -v_one_half],
                        [v_one_half, v_one_half, -v_one_half],
                        [aux_z_min, aux_z_min, -hd_center],
                        [hd_center, aux_column_width, -hd_center],
                        [v_one_half, aux_column_width, -sep_distance],
                        [aux_column_width, hd_center, -hd_center],
                        [hd_center, hd_center, -hd_center],
                        [v_one_half, sep_distance, -sep_distance],
                        [aux_column_width, v_one_half, -sep_distance],
                        [sep_distance, v_one_half, -sep_distance],
                        [v_one_half, v_one_half, -sep_distance],
                        [center_r, center_r, -center_r],
                        [hd_center, aux_x_max, -aux_x_max],
                        [v_one_half, x_max_r, -x_max_r],
                        [aux_y_max, hd_center, -aux_y_max],
                        [hd_center, hd_center, -aux_column_width],
                        [v_one_half, sep_distance, -aux_column_width],
                        [y_max_r, v_one_half, -y_max_r],
                        [sep_distance, v_one_half, -aux_column_width],
                        [v_one_half, v_one_half, -aux_column_width],
                    ]
                )
                + center_point
            )

            x_min_y_min_z_max = (
                _np.array(
                    [
                        [-v_one_half, -v_one_half, aux_column_width],
                        [-sep_distance, -v_one_half, aux_column_width],
                        [-y_min_r, -v_one_half, y_min_r],
                        [-v_one_half, -sep_distance, aux_column_width],
                        [-hd_center, -hd_center, aux_column_width],
                        [-aux_y_min, -hd_center, aux_y_min],
                        [-v_one_half, -x_min_r, x_min_r],
                        [-hd_center, -aux_x_min, aux_x_min],
                        [-center_r, -center_r, center_r],
                        [-v_one_half, -v_one_half, sep_distance],
                        [-sep_distance, -v_one_half, sep_distance],
                        [-aux_column_width, -v_one_half, sep_distance],
                        [-v_one_half, -sep_distance, sep_distance],
                        [-hd_center, -hd_center, hd_center],
                        [-aux_column_width, -hd_center, hd_center],
                        [-v_one_half, -aux_column_width, sep_distance],
                        [-hd_center, -aux_column_width, hd_center],
                        [-aux_z_max, -aux_z_max, hd_center],
                        [-v_one_half, -v_one_half, v_one_half],
                        [-sep_distance, -v_one_half, v_one_half],
                        [-aux_column_width, -v_one_half, v_one_half],
                        [-v_one_half, -sep_distance, v_one_half],
                        [-sep_distance, -sep_distance, v_one_half],
                        [-aux_column_width, -sep_distance, v_one_half],
                        [-v_one_half, -aux_column_width, v_one_half],
                        [-sep_distance, -aux_column_width, v_one_half],
                        [-z_max_r, -z_max_r, v_one_half],
                    ]
                )
                + center_point
            )

            x_max_y_min_z_max = (
                _np.array(
                    [
                        [y_min_r, -v_one_half, y_min_r],
                        [sep_distance, -v_one_half, aux_column_width],
                        [v_one_half, -v_one_half, aux_column_width],
                        [aux_y_min, -hd_center, aux_y_min],
                        [hd_center, -hd_center, aux_column_width],
                        [v_one_half, -sep_distance, aux_column_width],
                        [center_r, -center_r, center_r],
                        [hd_center, -aux_x_max, aux_x_max],
                        [v_one_half, -x_max_r, x_max_r],
                        [aux_column_width, -v_one_half, sep_distance],
                        [sep_distance, -v_one_half, sep_distance],
                        [v_one_half, -v_one_half, sep_distance],
                        [aux_column_width, -hd_center, hd_center],
                        [hd_center, -hd_center, hd_center],
                        [v_one_half, -sep_distance, sep_distance],
                        [aux_z_max, -aux_z_max, hd_center],
                        [hd_center, -aux_column_width, hd_center],
                        [v_one_half, -aux_column_width, sep_distance],
                        [aux_column_width, -v_one_half, v_one_half],
                        [sep_distance, -v_one_half, v_one_half],
                        [v_one_half, -v_one_half, v_one_half],
                        [aux_column_width, -sep_distance, v_one_half],
                        [sep_distance, -sep_distance, v_one_half],
                        [v_one_half, -sep_distance, v_one_half],
                        [z_max_r, -z_max_r, v_one_half],
                        [sep_distance, -aux_column_width, v_one_half],
                        [v_one_half, -aux_column_width, v_one_half],
                    ]
                )
                + center_point
            )

            x_min_y_max_z_max = (
                _np.array(
                    [
                        [-v_one_half, x_min_r, x_min_r],
                        [-hd_center, aux_x_min, aux_x_min],
                        [-center_r, center_r, center_r],
                        [-v_one_half, sep_distance, aux_column_width],
                        [-hd_center, hd_center, aux_column_width],
                        [-aux_y_max, hd_center, aux_y_max],
                        [-v_one_half, v_one_half, aux_column_width],
                        [-sep_distance, v_one_half, aux_column_width],
                        [-y_max_r, v_one_half, y_max_r],
                        [-v_one_half, aux_column_width, sep_distance],
                        [-hd_center, aux_column_width, hd_center],
                        [-aux_z_max, aux_z_max, hd_center],
                        [-v_one_half, sep_distance, sep_distance],
                        [-hd_center, hd_center, hd_center],
                        [-aux_column_width, hd_center, hd_center],
                        [-v_one_half, v_one_half, sep_distance],
                        [-sep_distance, v_one_half, sep_distance],
                        [-aux_column_width, v_one_half, sep_distance],
                        [-v_one_half, aux_column_width, v_one_half],
                        [-sep_distance, aux_column_width, v_one_half],
                        [-z_max_r, z_max_r, v_one_half],
                        [-v_one_half, sep_distance, v_one_half],
                        [-sep_distance, sep_distance, v_one_half],
                        [-aux_column_width, sep_distance, v_one_half],
                        [-v_one_half, v_one_half, v_one_half],
                        [-sep_distance, v_one_half, v_one_half],
                        [-aux_column_width, v_one_half, v_one_half],
                    ]
                )
                + center_point
            )

            x_max_y_max_z_max = (
                _np.array(
                    [
                        [center_r, center_r, center_r],
                        [hd_center, aux_x_max, aux_x_max],
                        [v_one_half, x_max_r, x_max_r],
                        [aux_y_max, hd_center, aux_y_max],
                        [hd_center, hd_center, aux_column_width],
                        [v_one_half, sep_distance, aux_column_width],
                        [y_max_r, v_one_half, y_max_r],
                        [sep_distance, v_one_half, aux_column_width],
                        [v_one_half, v_one_half, aux_column_width],
                        [aux_z_max, aux_z_max, hd_center],
                        [hd_center, aux_column_width, hd_center],
                        [v_one_half, aux_column_width, sep_distance],
                        [aux_column_width, hd_center, hd_center],
                        [hd_center, hd_center, hd_center],
                        [v_one_half, sep_distance, sep_distance],
                        [aux_column_width, v_one_half, sep_distance],
                        [sep_distance, v_one_half, sep_distance],
                        [v_one_half, v_one_half, sep_distance],
                        [z_max_r, z_max_r, v_one_half],
                        [sep_distance, aux_column_width, v_one_half],
                        [v_one_half, aux_column_width, v_one_half],
                        [aux_column_width, sep_distance, v_one_half],
                        [sep_distance, sep_distance, v_one_half],
                        [v_one_half, sep_distance, v_one_half],
                        [aux_column_width, v_one_half, v_one_half],
                        [sep_distance, v_one_half, v_one_half],
                        [v_one_half, v_one_half, v_one_half],
                    ]
                )
                + center_point
            )

            # Append the control points to the spline list
            for control_points, degrees in [
                (x_min_y_min, [2, 2, 1]),
                (x_max_y_min, [2, 2, 1]),
                (x_min_y_max, [2, 2, 1]),
                (x_max_y_max, [2, 2, 1]),
                (x_min_z_min, [2, 1, 2]),
                (x_max_z_min, [2, 1, 2]),
                (x_min_z_max, [2, 1, 2]),
                (x_max_z_max, [2, 1, 2]),
                (y_min_z_min, [1, 2, 2]),
                (y_max_z_min, [1, 2, 2]),
                (y_min_z_max, [1, 2, 2]),
                (y_max_z_max, [1, 2, 2]),
                (x_min_y_min_z_min, [2, 2, 2]),
                (x_max_y_min_z_min, [2, 2, 2]),
                (x_min_y_max_z_min, [2, 2, 2]),
                (x_max_y_max_z_min, [2, 2, 2]),
                (x_min_y_min_z_max, [2, 2, 2]),
                (x_max_y_min_z_max, [2, 2, 2]),
                (x_min_y_max_z_max, [2, 2, 2]),
                (x_max_y_max_z_max, [2, 2, 2]),
            ]:
                spline_list.append(
                    _Bezier(degrees=degrees, control_points=control_points)
                )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)
        return (splines, derivatives)
