import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class Armadillo(_TileBase):
    """
    Tile in the in a shape of a multi sided dice, where every side connects
    to the neighbor tile in the center of the surface.

    .. raw:: html

        <p><a href="../_static/Armadillo.html">Fullscreen</a>.</p>
        <embed type="text/html" width="100%" height="400" src="../_static/Armadillo.html" />

    """  # noqa: E501

    _para_dim = 3
    _dim = 3
    _evaluation_points = _np.array([[0.5, 0.5, 0.5]])
    _n_info_per_eval_point = 1
    _sensitivities_implemented = True
    _closure_directions = [
        "x_min",
        "x_max",
        "y_min",
        "y_max",
        "z_min",
        "z_max",
    ]
    _parameter_bounds = [[0.0, 0.5]]
    _parameters_shape = (1, 1)
    _default_parameter_value = 0.2

    def _closing_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        contact_length=0.3,
        closure=None,
        **kwargs,  # noqa ARG002
    ):
        """Create a closing tile to match with closed surface.

        The tile will look like this.

        Parameters
        ----------
        parameters: np.ndarray
            An evaluation point with one parameter is used. This parameter
            describes the thickness of the wall. The parameters must be a
            two-dimensional np.array, where the value must be between 0.01
            and 0.49
        parameter_sensitivities: np.ndarray
          Describes the parameter sensitivities with respect to some design
          variable. In case the design variables directly apply to the
          parameter itself, they evaluate as delta_ij
        contact_length: float
          the length of the wall that contacts the other microstructure
        closure : str
          parametric dimension that needs to be closed, given in the form
          "x_min", "x_max", etc.

        Results
        -------
        spline_list : list
        derivative_list : list / None
        """
        if closure is None:
            raise ValueError("No closing direction given")

        self._check_custom_parameter(
            contact_length, "contact length", 0.0, 0.99
        )

        parameters, n_derivatives, derivatives = self._process_input(
            parameters=parameters,
            parameter_sensitivities=parameter_sensitivities,
        )

        splines = []
        for i_derivative in range(n_derivatives + 1):
            if i_derivative == 0:
                v_wall_thickness = parameters[0, 0]
                v_zero = 0.0
                v_one_half = 0.5
                v_one = 1.0
                v_half_contact_length = contact_length * 0.5
                v_inner_half_contact_length = contact_length * v_wall_thickness
            else:
                v_wall_thickness = parameter_sensitivities[
                    0, 0, i_derivative - 1
                ]
                v_zero = 0.0
                v_one_half = 0.0
                v_one = 0.0
                v_half_contact_length = 0.0
                v_inner_half_contact_length = contact_length * v_wall_thickness

            spline_list = []

            if closure == "x_min":
                # set points:
                right = _np.array(
                    [
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                connection_front_right = _np.array(
                    [
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                front = _np.array(
                    [
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                connection_back_left = _np.array(
                    [
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_one,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_zero,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                left = _np.array(
                    [
                        [
                            v_zero,
                            v_one,
                            v_one,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_one,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one,
                            v_zero,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_zero,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_front_left = _np.array(
                    [
                        [
                            v_zero,
                            v_one,
                            v_one,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one,
                            v_zero,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                back = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_back_right = _np.array(
                    [
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                bottom = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                top = _np.array(
                    [
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                    ]
                )

                connection_front_bottom = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                connection_front_top = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                    ]
                )

                connection_back_bottom = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_back_top = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_top_right = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one,
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_top_left = _np.array(
                    [
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_zero,
                            v_one,
                            v_one,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_one,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_bottom_left = _np.array(
                    [
                        [
                            v_zero,
                            v_one,
                            v_zero,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                connection_bottom_right = _np.array(
                    [
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

            elif closure == "x_max":
                # set points:
                right = _np.array(
                    [
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_one,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_zero,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one,
                            v_one,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one,
                            v_zero,
                        ],
                    ]
                )

                connection_front_right = _np.array(
                    [
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one,
                            v_one,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one,
                            v_zero,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                front = _np.array(
                    [
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                connection_back_left = _np.array(
                    [
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                left = _np.array(
                    [
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_front_left = _np.array(
                    [
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                back = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_back_right = _np.array(
                    [
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_one,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_zero,
                        ],
                    ]
                )

                bottom = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                top = _np.array(
                    [
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                    ]
                )

                connection_front_bottom = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                connection_front_top = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                    ]
                )

                connection_back_bottom = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_back_top = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_top_right = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one,
                            v_one,
                            v_one,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_one,
                        ],
                        [
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_top_left = _np.array(
                    [
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_zero,
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_bottom_left = _np.array(
                    [
                        [
                            v_zero,
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                connection_bottom_right = _np.array(
                    [
                        [
                            v_one,
                            v_one,
                            v_zero,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_zero,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

            elif closure == "y_min":
                # set points:
                right = _np.array(
                    [
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                connection_front_right = _np.array(
                    [
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                front = _np.array(
                    [
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                connection_back_left = _np.array(
                    [
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_one,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_zero,
                        ],
                    ]
                )

                left = _np.array(
                    [
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_front_left = _np.array(
                    [
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                back = _np.array(
                    [
                        [
                            v_one,
                            v_zero,
                            v_one,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_zero,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_one,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_zero,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_back_right = _np.array(
                    [
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_one,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_zero,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                bottom = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                top = _np.array(
                    [
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                    ]
                )

                connection_front_bottom = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                connection_front_top = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                    ]
                )

                connection_back_bottom = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_zero,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_zero,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_back_top = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_one,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_one,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_top_right = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one,
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_top_left = _np.array(
                    [
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_zero,
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_bottom_left = _np.array(
                    [
                        [
                            v_zero,
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                connection_bottom_right = _np.array(
                    [
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

            elif closure == "y_max":
                # set points:
                # set points:
                right = _np.array(
                    [
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                connection_front_right = _np.array(
                    [
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one,
                            v_one,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one,
                            v_zero,
                        ],
                    ]
                )

                front = _np.array(
                    [
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one,
                            v_one,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one,
                            v_zero,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one,
                            v_one,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one,
                            v_zero,
                        ],
                    ]
                )

                connection_back_left = _np.array(
                    [
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                left = _np.array(
                    [
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_top_left = _np.array(
                    [
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one,
                            v_one,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one,
                            v_zero,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                back = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_back_right = _np.array(
                    [
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                bottom = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                top = _np.array(
                    [
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                    ]
                )

                connection_front_bottom = _np.array(
                    [
                        [
                            v_one,
                            v_one,
                            v_zero,
                        ],
                        [
                            v_zero,
                            v_one,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                connection_front_top = _np.array(
                    [
                        [
                            v_one,
                            v_one,
                            v_one,
                        ],
                        [
                            v_zero,
                            v_one,
                            v_one,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                    ]
                )

                connection_back_bottom = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_back_top = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_top_right = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one,
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_front_left = _np.array(
                    [
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_zero,
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_bottom_left = _np.array(
                    [
                        [
                            v_zero,
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                connection_bottom_right = _np.array(
                    [
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

            elif closure == "z_max":
                # set points:
                right = _np.array(
                    [
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                connection_front_right = _np.array(
                    [
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                front = _np.array(
                    [
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                connection_back_left = _np.array(
                    [
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                left = _np.array(
                    [
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_front_left = _np.array(
                    [
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                back = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_back_right = _np.array(
                    [
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                bottom = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                top = _np.array(
                    [
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one,
                            v_one,
                            v_one,
                        ],
                        [
                            v_zero,
                            v_one,
                            v_one,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_one,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_one,
                        ],
                    ]
                )

                connection_front_bottom = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                connection_front_top = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one,
                            v_one,
                        ],
                        [
                            v_zero,
                            v_one,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                    ]
                )

                connection_back_bottom = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_back_top = _np.array(
                    [
                        [
                            v_one,
                            v_zero,
                            v_one,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_top_right = _np.array(
                    [
                        [
                            v_one,
                            v_one,
                            v_one,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one,
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_top_left = _np.array(
                    [
                        [
                            v_zero,
                            v_one,
                            v_one,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_one,
                        ],
                        [
                            v_zero,
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_bottom_left = _np.array(
                    [
                        [
                            v_zero,
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                connection_bottom_right = _np.array(
                    [
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_zero,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

            elif closure == "z_min":
                # set points:
                # set points:
                right = _np.array(
                    [
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                connection_front_right = _np.array(
                    [
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                front = _np.array(
                    [
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                connection_back_left = _np.array(
                    [
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                left = _np.array(
                    [
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_front_left = _np.array(
                    [
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                back = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_back_right = _np.array(
                    [
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            -v_wall_thickness + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                    ]
                )

                bottom = _np.array(
                    [
                        [
                            v_one,
                            v_one,
                            v_zero,
                        ],
                        [
                            v_zero,
                            v_one,
                            v_zero,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_zero,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_zero,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                top = _np.array(
                    [
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                    ]
                )

                connection_front_bottom = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one,
                            v_zero,
                        ],
                        [
                            v_zero,
                            v_one,
                            v_zero,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                connection_front_top = _np.array(
                    [
                        [
                            v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_one,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                    ]
                )

                connection_back_bottom = _np.array(
                    [
                        [
                            v_one,
                            v_zero,
                            v_zero,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_zero,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_zero,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_back_top = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            -v_half_contact_length + v_one_half,
                            v_zero,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_top_right = _np.array(
                    [
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one,
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one_half + v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half + v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_top_left = _np.array(
                    [
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                            v_one,
                        ],
                        [
                            v_zero,
                            v_one_half + v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one_half - v_half_contact_length,
                            v_one_half + v_half_contact_length,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_wall_thickness,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                        ],
                    ]
                )

                connection_bottom_left = _np.array(
                    [
                        [
                            v_zero,
                            v_one_half + v_half_contact_length,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one_half - v_half_contact_length,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_zero,
                            v_one,
                            v_zero,
                        ],
                        [
                            v_zero,
                            v_zero,
                            v_zero,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_wall_thickness,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

                connection_bottom_right = _np.array(
                    [
                        [
                            v_one,
                            v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_one,
                            -v_half_contact_length + v_one_half,
                            v_one_half - v_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_wall_thickness + v_one_half,
                            -v_inner_half_contact_length + v_one_half,
                            v_one_half - v_inner_half_contact_length,
                        ],
                        [
                            v_one,
                            v_one,
                            v_zero,
                        ],
                        [
                            v_one,
                            v_zero,
                            v_zero,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                        [
                            v_one_half + v_inner_half_contact_length,
                            v_one_half - v_inner_half_contact_length,
                            v_one_half - v_wall_thickness,
                        ],
                    ]
                )

            for control_points in [
                right,
                connection_front_right,
                front,
                connection_back_left,
                left,
                connection_front_left,
                back,
                connection_back_right,
                bottom,
                top,
                connection_front_bottom,
                connection_front_top,
                connection_back_bottom,
                connection_back_top,
                connection_top_right,
                connection_top_left,
                connection_bottom_left,
                connection_bottom_right,
            ]:
                spline_list.append(
                    _Bezier(degrees=[1, 1, 1], control_points=control_points)
                )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        return (splines, derivatives)

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,  # TODO
        contact_length=0.3,
        closure=None,
        **kwargs,  # noqa ARG002
    ):
        """Create a microtile based on the parameters that describe the wall
        thicknesses.

        Thickness parameters are used to describe the inner radius of the
        outward facing branches

        Parameters
        ----------
        parameters : np.array
            One evaluation point with one parameter is used. This parameter
            describes the thickness of the wall. The parameters must be a
            two-dimensional np.array, where the value must be between 0.01
            and 0.49
        parameter_sensitivities: np.ndarray
          Describes the parameter sensitivities with respect to some design
          variable. In case the design variables directly apply to the
          parameter itself, they evaluate as delta_ij
        contact_length : float
            the length of the wall that contacts the other microstructure
        closure : str
          parametric dimension that needs to be closed, given in the form
          "x_min", "x_max", etc.

        Returns
        -------
        microtile_list : list(splines)
        derivative_list : list / None
        """

        self._check_custom_parameter(
            contact_length, "contact length", 0.0, 0.99
        )

        parameters, n_derivatives, derivatives = self._process_input(
            parameters=parameters,
            parameter_sensitivities=parameter_sensitivities,
        )

        if closure is not None:
            return self._closing_tile(
                parameters=parameters,
                parameter_sensitivities=parameter_sensitivities,
                contact_length=contact_length,
                closure=closure,
                **kwargs,
            )

        splines = []
        for i_derivative in range(n_derivatives + 1):
            if i_derivative == 0:
                v_wall_thickness = parameters[0, 0]
                v_zero = 0.0
                v_one_half = 0.5
                v_one = 1.0
                v_half_contact_length = contact_length * 0.5
                v_inner_half_contact_length = contact_length * v_wall_thickness
            else:
                v_wall_thickness = parameter_sensitivities[
                    0, 0, i_derivative - 1
                ]
                v_zero = 0.0
                v_one_half = 0.0
                v_one = 0.0
                v_half_contact_length = 0.0
                v_inner_half_contact_length = contact_length * v_wall_thickness

            spline_list = []

            # set points:
            right = _np.array(
                [
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            connection_front_right = _np.array(
                [
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            front = _np.array(
                [
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            connection_back_left = _np.array(
                [
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            left = _np.array(
                [
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_front_left = _np.array(
                [
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            back = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_back_right = _np.array(
                [
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            bottom = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            top = _np.array(
                [
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                ]
            )

            connection_front_bottom = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            connection_front_top = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                ]
            )

            connection_back_bottom = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_back_top = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_top_right = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one,
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_top_left = _np.array(
                [
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_bottom_left = _np.array(
                [
                    [
                        v_zero,
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            connection_bottom_right = _np.array(
                [
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            for control_points in [
                right,
                connection_front_right,
                front,
                connection_back_left,
                left,
                connection_front_left,
                back,
                connection_back_right,
                bottom,
                top,
                connection_front_bottom,
                connection_front_top,
                connection_back_bottom,
                connection_back_top,
                connection_top_right,
                connection_top_left,
                connection_bottom_left,
                connection_bottom_right,
            ]:
                spline_list.append(
                    _Bezier(degrees=[1, 1, 1], control_points=control_points)
                )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        return (splines, derivatives)
