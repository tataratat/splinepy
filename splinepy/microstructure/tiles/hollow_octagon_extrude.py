import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class HollowOctagonExtrude(_TileBase):
    """Simple tile - looks like a nut, but in 3D.

    .. raw:: html

        <p><a href="../_static/HollowOctagonExtrude.html">Fullscreen</a>.</p>
        <embed type="text/html" width="100%" height="400" src="../_static/HollowOctagonExtrude.html" />

    """  # noqa: E501

    _dim = 3
    _para_dim = 3
    _evaluation_points = _np.array([[0.5, 0.5, 0.5]])
    _n_info_per_eval_point = 1
    _sensitivities_implemented = True
    _closure_directions = ["x_min", "x_max", "y_min", "y_max"]
    _parameter_bounds = [[0.0, 0.5]]
    _parameters_shape = (1, 1)
    _default_parameter_value = 0.2

    _CONTACT_LENGTH_BOUNDS = [0.0, 0.99]

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        contact_length=0.2,
        closure=None,
        **kwargs,  # noqa ARG002
    ):
        """Create a microtile based on the parameters that describe the wall
        thicknesses.

        Thickness parameters are used to describe the inner radius of the
        outward facing branches

        Parameters
        ----------
        parameters : np.array(1, 1)
          One evaluation point with one parameter is used. This parameter
          specifies the distance from the center to the inner edge, where
          the value must be between 0.01 and 0.49.
        parameter_sensitivities: np.ndarray
          Describes the parameter sensitivities with respect to some design
          variable. In case the design variables directly apply to the
          parameter itself, they evaluate as delta_ij
        contact_length : float
            the length of the wall that contacts the other microstructure
        Returns
        -------
        microtile_list : list(splines)
        derivatives: list<list<splines>> / None
        """

        self._check_custom_parameter(
            contact_length, "contact length", self._CONTACT_LENGTH_BOUNDS
        )
        # Process input
        parameters, n_derivatives, derivatives = self._process_input(
            parameters=parameters,
            parameter_sensitivities=parameter_sensitivities,
        )

        v_h_void = parameters[0, 0]

        if closure is not None:
            return self._closing_tile(
                parameters=parameters,
                parameter_sensitivities=parameter_sensitivities,
                closure=closure,
                contact_length=contact_length,
                **kwargs,
            )

        splines = []
        for i_derivative in range(n_derivatives + 1):
            if i_derivative == 0:
                v_zero = 0.0
                v_one_half = 0.5
                v_one = 1.0
                v_outer_c_h = contact_length * 0.5
                v_inner_c_h = contact_length * v_h_void
            else:
                v_zero = 0.0
                v_one_half = 0.0
                v_one = 0.0
                v_h_void = parameter_sensitivities[0, 0, i_derivative - 1]
                v_outer_c_h = 0.0
                v_inner_c_h = contact_length * v_h_void

            spline_list = []

            # set points:
            right = _np.array(
                [
                    [v_h_void + v_one_half, -v_inner_c_h + v_one_half, v_zero],
                    [v_one, -v_outer_c_h + v_one_half, v_zero],
                    [v_h_void + v_one_half, v_inner_c_h + v_one_half, v_zero],
                    [v_one, v_outer_c_h + v_one_half, v_zero],
                    [v_h_void + v_one_half, -v_inner_c_h + v_one_half, v_one],
                    [v_one, -v_outer_c_h + v_one_half, v_one],
                    [v_h_void + v_one_half, v_inner_c_h + v_one_half, v_one],
                    [v_one, v_outer_c_h + v_one_half, v_one],
                ]
            )

            right_top = _np.array(
                [
                    [v_h_void + v_one_half, v_inner_c_h + v_one_half, v_zero],
                    [v_one, v_outer_c_h + v_one_half, v_zero],
                    [v_inner_c_h + v_one_half, v_h_void + v_one_half, v_zero],
                    [v_outer_c_h + v_one_half, v_one, v_zero],
                    [v_h_void + v_one_half, v_inner_c_h + v_one_half, v_one],
                    [v_one, v_outer_c_h + v_one_half, v_one],
                    [v_inner_c_h + v_one_half, v_h_void + v_one_half, v_one],
                    [v_outer_c_h + v_one_half, v_one, v_one],
                ]
            )

            top = _np.array(
                [
                    [v_inner_c_h + v_one_half, v_h_void + v_one_half, v_zero],
                    [v_outer_c_h + v_one_half, v_one, v_zero],
                    [-v_inner_c_h + v_one_half, v_h_void + v_one_half, v_zero],
                    [-v_outer_c_h + v_one_half, v_one, v_zero],
                    [v_inner_c_h + v_one_half, v_h_void + v_one_half, v_one],
                    [v_outer_c_h + v_one_half, v_one, v_one],
                    [-v_inner_c_h + v_one_half, v_h_void + v_one_half, v_one],
                    [-v_outer_c_h + v_one_half, v_one, v_one],
                ]
            )

            bottom_left = _np.array(
                [
                    [
                        -v_h_void + v_one_half,
                        -v_inner_c_h + v_one_half,
                        v_zero,
                    ],
                    [v_zero, -v_outer_c_h + v_one_half, v_zero],
                    [
                        -v_inner_c_h + v_one_half,
                        -v_h_void + v_one_half,
                        v_zero,
                    ],
                    [-v_outer_c_h + v_one_half, v_zero, v_zero],
                    [-v_h_void + v_one_half, -v_inner_c_h + v_one_half, v_one],
                    [v_zero, -v_outer_c_h + v_one_half, v_one],
                    [-v_inner_c_h + v_one_half, -v_h_void + v_one_half, v_one],
                    [-v_outer_c_h + v_one_half, v_zero, v_one],
                ]
            )

            left = _np.array(
                [
                    [v_zero, -v_outer_c_h + v_one_half, v_zero],
                    [
                        -v_h_void + v_one_half,
                        -v_inner_c_h + v_one_half,
                        v_zero,
                    ],
                    [v_zero, v_outer_c_h + v_one_half, v_zero],
                    [-v_h_void + v_one_half, v_inner_c_h + v_one_half, v_zero],
                    [v_zero, -v_outer_c_h + v_one_half, v_one],
                    [-v_h_void + v_one_half, -v_inner_c_h + v_one_half, v_one],
                    [v_zero, v_outer_c_h + v_one_half, v_one],
                    [-v_h_void + v_one_half, v_inner_c_h + v_one_half, v_one],
                ]
            )

            top_left = _np.array(
                [
                    [v_zero, v_outer_c_h + v_one_half, v_zero],
                    [-v_h_void + v_one_half, v_inner_c_h + v_one_half, v_zero],
                    [-v_outer_c_h + v_one_half, v_one, v_zero],
                    [-v_inner_c_h + v_one_half, v_h_void + v_one_half, v_zero],
                    [v_zero, v_outer_c_h + v_one_half, v_one],
                    [-v_h_void + v_one_half, v_inner_c_h + v_one_half, v_one],
                    [-v_outer_c_h + v_one_half, v_one, v_one],
                    [-v_inner_c_h + v_one_half, v_h_void + v_one_half, v_one],
                ]
            )

            bottom = _np.array(
                [
                    [v_outer_c_h + v_one_half, v_zero, v_zero],
                    [v_inner_c_h + v_one_half, -v_h_void + v_one_half, v_zero],
                    [-v_outer_c_h + v_one_half, v_zero, v_zero],
                    [
                        -v_inner_c_h + v_one_half,
                        -v_h_void + v_one_half,
                        v_zero,
                    ],
                    [v_outer_c_h + v_one_half, v_zero, v_one],
                    [v_inner_c_h + v_one_half, -v_h_void + v_one_half, v_one],
                    [-v_outer_c_h + v_one_half, v_zero, v_one],
                    [-v_inner_c_h + v_one_half, -v_h_void + v_one_half, v_one],
                ]
            )

            bottom_right = _np.array(
                [
                    [v_inner_c_h + v_one_half, -v_h_void + v_one_half, v_zero],
                    [v_outer_c_h + v_one_half, v_zero, v_zero],
                    [v_h_void + v_one_half, -v_inner_c_h + v_one_half, v_zero],
                    [v_one, -v_outer_c_h + v_one_half, v_zero],
                    [v_inner_c_h + v_one_half, -v_h_void + v_one_half, v_one],
                    [v_outer_c_h + v_one_half, v_zero, v_one],
                    [v_h_void + v_one_half, -v_inner_c_h + v_one_half, v_one],
                    [v_one, -v_outer_c_h + v_one_half, v_one],
                ]
            )

            for control_points in [
                right,
                right_top,
                top,
                bottom_left,
                left,
                top_left,
                bottom,
                bottom_right,
            ]:
                spline_list.append(
                    _Bezier(degrees=[1, 1, 1], control_points=control_points)
                )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        return (splines, derivatives)

    def _closing_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        contact_length=0.2,
        closure=None,
    ):
        """Create a closing tile to match with closed surface.

        Parameters
        ----------
        parameters: np.ndarray(1, 1)
          One evaluation point with one parameter is used. This parameter
          specifies the distance from the center to the inner edge, where
          the value must be between 0.01 and 0.49.
        parameter_sensitivities: np.ndarray
          Describes the parameter sensitivities with respect to some design
          variable. In case the design variables directly apply to the
          parameter itself, they evaluate as delta_ij
        closure : int
          parametric dimension that needs to be closed. Positive values mean
          that minimum parametric dimension is requested. That means,
          i.e. -2 closes the tile at maximum z-coordinate.
          (must currently be either -2 or 2)
        contact_length: float
          the length of the wall that contacts the other microstructure

        Returns
        -------
        spline_list : list
        derivatives: list<list<splines>> / None
        """
        if closure is None:
            raise ValueError("No closing direction given")

        # Process input
        parameters, n_derivatives, derivatives = self._process_input(
            parameters=parameters,
            parameter_sensitivities=parameter_sensitivities,
        )

        v_h_void = parameters[0, 0]

        splines = []
        for i_derivative in range(n_derivatives + 1):
            if i_derivative == 0:
                v_zero = 0.0
                v_one_half = 0.5
                v_one = 1.0
                v_outer_c_h = contact_length * 0.5
                v_inner_c_h = contact_length * parameters[0, 0]
            else:
                v_h_void = parameter_sensitivities[0, 0, i_derivative - 1]
                v_zero, v_one_half, v_one = [0.0, 0.0, 0.0]
                v_outer_c_h = 0.0
                v_inner_c_h = contact_length * v_h_void

            spline_list = []

            if closure == "x_min":
                # set points:
                right = _np.array(
                    [
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, -v_outer_c_h + v_one_half, v_zero],
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, v_outer_c_h + v_one_half, v_zero],
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, -v_outer_c_h + v_one_half, v_one],
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, v_outer_c_h + v_one_half, v_one],
                    ]
                )

                right_top = _np.array(
                    [
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, v_outer_c_h + v_one_half, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_outer_c_h + v_one_half, v_one, v_zero],
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, v_outer_c_h + v_one_half, v_one],
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_outer_c_h + v_one_half, v_one, v_one],
                    ]
                )

                top = _np.array(
                    [
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_outer_c_h + v_one_half, v_one, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [-v_outer_c_h + v_one_half, v_one, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_outer_c_h + v_one_half, v_one, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                        [-v_outer_c_h + v_one_half, v_one, v_one],
                    ]
                )

                bottom_left = _np.array(
                    [
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_zero, v_zero, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [-v_outer_c_h + v_one_half, v_zero, v_zero],
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_zero, v_zero, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                        [-v_outer_c_h + v_one_half, v_zero, v_one],
                    ]
                )

                left = _np.array(
                    [
                        [v_zero, v_zero, v_zero],
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_zero, v_one, v_zero],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_zero, v_zero, v_one],
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_zero, v_one, v_one],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                    ]
                )

                top_left = _np.array(
                    [
                        [v_zero, v_one, v_zero],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [-v_outer_c_h + v_one_half, v_one, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_zero, v_one, v_one],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [-v_outer_c_h + v_one_half, v_one, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                    ]
                )

                bottom = _np.array(
                    [
                        [v_outer_c_h + v_one_half, v_zero, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [-v_outer_c_h + v_one_half, v_zero, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_outer_c_h + v_one_half, v_zero, v_one],
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                        [-v_outer_c_h + v_one_half, v_zero, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                    ]
                )

                bottom_right = _np.array(
                    [
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_outer_c_h + v_one_half, v_zero, v_zero],
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, -v_outer_c_h + v_one_half, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_outer_c_h + v_one_half, v_zero, v_one],
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, -v_outer_c_h + v_one_half, v_one],
                    ]
                )

            elif closure == "x_max":
                right = _np.array(
                    [
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, v_zero, v_zero],
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, v_one, v_zero],
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, v_zero, v_one],
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, v_one, v_one],
                    ]
                )

                right_top = _np.array(
                    [
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, v_one, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_outer_c_h + v_one_half, v_one, v_zero],
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, v_one, v_one],
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_outer_c_h + v_one_half, v_one, v_one],
                    ]
                )

                top = _np.array(
                    [
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_outer_c_h + v_one_half, v_one, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [-v_outer_c_h + v_one_half, v_one, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_outer_c_h + v_one_half, v_one, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                        [-v_outer_c_h + v_one_half, v_one, v_one],
                    ]
                )

                bottom_left = _np.array(
                    [
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_zero, -v_outer_c_h + v_one_half, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [-v_outer_c_h + v_one_half, v_zero, v_zero],
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_zero, -v_outer_c_h + v_one_half, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                        [-v_outer_c_h + v_one_half, v_zero, v_one],
                    ]
                )

                left = _np.array(
                    [
                        [v_zero, -v_outer_c_h + v_one_half, v_zero],
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_zero, v_outer_c_h + v_one_half, v_zero],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_zero, -v_outer_c_h + v_one_half, v_one],
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_zero, v_outer_c_h + v_one_half, v_one],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                    ]
                )

                top_left = _np.array(
                    [
                        [v_zero, v_outer_c_h + v_one_half, v_zero],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [-v_outer_c_h + v_one_half, v_one, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_zero, v_outer_c_h + v_one_half, v_one],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [-v_outer_c_h + v_one_half, v_one, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                    ]
                )

                bottom = _np.array(
                    [
                        [v_outer_c_h + v_one_half, v_zero, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [-v_outer_c_h + v_one_half, v_zero, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_outer_c_h + v_one_half, v_zero, v_one],
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                        [-v_outer_c_h + v_one_half, v_zero, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                    ]
                )

                bottom_right = _np.array(
                    [
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_outer_c_h + v_one_half, v_zero, v_zero],
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, v_zero, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_outer_c_h + v_one_half, v_zero, v_one],
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, v_zero, v_one],
                    ]
                )

            elif closure == "y_min":
                # set points:
                right = _np.array(
                    [
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, -v_outer_c_h + v_one_half, v_zero],
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, v_outer_c_h + v_one_half, v_zero],
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, -v_outer_c_h + v_one_half, v_one],
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, v_outer_c_h + v_one_half, v_one],
                    ]
                )

                right_top = _np.array(
                    [
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, v_outer_c_h + v_one_half, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_outer_c_h + v_one_half, v_one, v_zero],
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, v_outer_c_h + v_one_half, v_one],
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_outer_c_h + v_one_half, v_one, v_one],
                    ]
                )

                top = _np.array(
                    [
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_outer_c_h + v_one_half, v_one, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [-v_outer_c_h + v_one_half, v_one, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_outer_c_h + v_one_half, v_one, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                        [-v_outer_c_h + v_one_half, v_one, v_one],
                    ]
                )

                bottom_left = _np.array(
                    [
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_zero, -v_outer_c_h + v_one_half, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_zero, v_zero, v_zero],
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_zero, -v_outer_c_h + v_one_half, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_zero, v_zero, v_one],
                    ]
                )

                left = _np.array(
                    [
                        [v_zero, -v_outer_c_h + v_one_half, v_zero],
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_zero, v_outer_c_h + v_one_half, v_zero],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_zero, -v_outer_c_h + v_one_half, v_one],
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_zero, v_outer_c_h + v_one_half, v_one],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                    ]
                )

                top_left = _np.array(
                    [
                        [v_zero, v_outer_c_h + v_one_half, v_zero],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [-v_outer_c_h + v_one_half, v_one, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_zero, v_outer_c_h + v_one_half, v_one],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [-v_outer_c_h + v_one_half, v_one, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                    ]
                )

                bottom = _np.array(
                    [
                        [v_one, v_zero, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_zero, v_zero, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_one, v_zero, v_one],
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_zero, v_zero, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                    ]
                )

                bottom_right = _np.array(
                    [
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_one, v_zero, v_zero],
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, -v_outer_c_h + v_one_half, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_one, v_zero, v_one],
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, -v_outer_c_h + v_one_half, v_one],
                    ]
                )

            elif closure == "y_max":
                # set points:
                right = _np.array(
                    [
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, -v_outer_c_h + v_one_half, v_zero],
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, v_outer_c_h + v_one_half, v_zero],
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, -v_outer_c_h + v_one_half, v_one],
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, v_outer_c_h + v_one_half, v_one],
                    ]
                )

                right_top = _np.array(
                    [
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, v_outer_c_h + v_one_half, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_one, v_one, v_zero],
                        [
                            v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, v_outer_c_h + v_one_half, v_one],
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_one, v_one, v_one],
                    ]
                )

                top = _np.array(
                    [
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_one, v_one, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_zero, v_one, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_one, v_one, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_zero, v_one, v_one],
                    ]
                )

                bottom_left = _np.array(
                    [
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_zero, -v_outer_c_h + v_one_half, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [-v_outer_c_h + v_one_half, v_zero, v_zero],
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_zero, -v_outer_c_h + v_one_half, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                        [-v_outer_c_h + v_one_half, v_zero, v_one],
                    ]
                )

                left = _np.array(
                    [
                        [v_zero, -v_outer_c_h + v_one_half, v_zero],
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_zero, v_outer_c_h + v_one_half, v_zero],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_zero, -v_outer_c_h + v_one_half, v_one],
                        [
                            -v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_zero, v_outer_c_h + v_one_half, v_one],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                    ]
                )

                top_left = _np.array(
                    [
                        [v_zero, v_outer_c_h + v_one_half, v_zero],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_zero, v_one, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_zero, v_outer_c_h + v_one_half, v_one],
                        [
                            -v_h_void + v_one_half,
                            v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_zero, v_one, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            v_h_void + v_one_half,
                            v_one,
                        ],
                    ]
                )

                bottom = _np.array(
                    [
                        [v_outer_c_h + v_one_half, v_zero, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [-v_outer_c_h + v_one_half, v_zero, v_zero],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_outer_c_h + v_one_half, v_zero, v_one],
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                        [-v_outer_c_h + v_one_half, v_zero, v_one],
                        [
                            -v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                    ]
                )

                bottom_right = _np.array(
                    [
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_zero,
                        ],
                        [v_outer_c_h + v_one_half, v_zero, v_zero],
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_zero,
                        ],
                        [v_one, -v_outer_c_h + v_one_half, v_zero],
                        [
                            v_inner_c_h + v_one_half,
                            -v_h_void + v_one_half,
                            v_one,
                        ],
                        [v_outer_c_h + v_one_half, v_zero, v_one],
                        [
                            v_h_void + v_one_half,
                            -v_inner_c_h + v_one_half,
                            v_one,
                        ],
                        [v_one, -v_outer_c_h + v_one_half, v_one],
                    ]
                )

            for control_points in [
                right,
                right_top,
                bottom,
                bottom_left,
                left,
                top_left,
                top,
                bottom_right,
            ]:
                spline_list.append(
                    _Bezier(degrees=[1, 1, 1], control_points=control_points)
                )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        return (splines, derivatives)
