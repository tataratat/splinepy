import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class HollowOctagon(_TileBase):
    """Simple tile - looks like a nut

    .. raw:: html

        <p><a href="../_static/HollowOctagon.html">Fullscreen</a>.</p>
        <embed type="text/html" width="100%" height="400" src="../_static/HollowOctagon.html" />

    """  # noqa: E501

    _dim = 2
    _para_dim = 2
    _evaluation_points = _np.array([[0.5, 0.5]])
    _n_info_per_eval_point = 1
    _sensitivities_implemented = True
    _closure_directions = ["x_min", "x_max", "y_min", "y_max"]
    _parameter_bounds = [[0.0, 0.5]]
    _parameters_shape = (1, 1)
    _default_parameter_value = 0.2

    _CONTACT_LENGTH_BOUNDS = [0.0, 0.99]

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
        contact_length: float
          the length of the wall that contacts the other microstructure
        closure : str
          parametric dimension that needs to be closed, given in the form
          "x_min", "x_max", etc.

        Results
        -------
        spline_list : list
        """
        if closure is None:
            raise ValueError("No closing direction given")

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
                v_inner_c_h = contact_length * v_h_void
            else:
                v_zero = 0.0
                v_one_half = 0.0
                v_one = 0.0
                v_h_void = parameter_sensitivities[0, 0, i_derivative - 1]
                v_outer_c_h = 0.0
                v_inner_c_h = contact_length * v_h_void

            spline_list = []

            if closure == "x_min":
                # set points:
                right = _np.array(
                    [
                        [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_one, -v_outer_c_h + v_one_half],
                        [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                        [v_one, v_outer_c_h + v_one_half],
                    ]
                )

                right_top = _np.array(
                    [
                        [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                        [v_one, v_outer_c_h + v_one_half],
                        [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                        [v_outer_c_h + v_one_half, v_one],
                    ]
                )

                top = _np.array(
                    [
                        [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                        [v_outer_c_h + v_one_half, v_one],
                        [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                        [-v_outer_c_h + v_one_half, v_one],
                    ]
                )

                bottom_left = _np.array(
                    [
                        [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_zero, v_zero],
                        [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                        [-v_outer_c_h + v_one_half, v_zero],
                    ]
                )

                left = _np.array(
                    [
                        [v_zero, v_zero],
                        [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_zero, v_one],
                        [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    ]
                )

                top_left = _np.array(
                    [
                        [v_zero, v_one],
                        [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                        [-v_outer_c_h + v_one_half, v_one],
                        [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    ]
                )

                bottom = _np.array(
                    [
                        [v_outer_c_h + v_one_half, v_zero],
                        [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                        [-v_outer_c_h + v_one_half, v_zero],
                        [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    ]
                )

                bottom_right = _np.array(
                    [
                        [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                        [v_outer_c_h + v_one_half, v_zero],
                        [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_one, -v_outer_c_h + v_one_half],
                    ]
                )

            elif closure == "x_max":
                right = _np.array(
                    [
                        [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_one, v_zero],
                        [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                        [v_one, v_one],
                    ]
                )

                right_top = _np.array(
                    [
                        [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                        [v_one, v_one],
                        [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                        [v_outer_c_h + v_one_half, v_one],
                    ]
                )

                top = _np.array(
                    [
                        [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                        [v_outer_c_h + v_one_half, v_one],
                        [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                        [-v_outer_c_h + v_one_half, v_one],
                    ]
                )

                bottom_left = _np.array(
                    [
                        [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_zero, -v_outer_c_h + v_one_half],
                        [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                        [-v_outer_c_h + v_one_half, v_zero],
                    ]
                )

                left = _np.array(
                    [
                        [v_zero, -v_outer_c_h + v_one_half],
                        [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_zero, v_outer_c_h + v_one_half],
                        [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    ]
                )

                top_left = _np.array(
                    [
                        [v_zero, v_outer_c_h + v_one_half],
                        [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                        [-v_outer_c_h + v_one_half, v_one],
                        [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    ]
                )

                bottom = _np.array(
                    [
                        [v_outer_c_h + v_one_half, v_zero],
                        [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                        [-v_outer_c_h + v_one_half, v_zero],
                        [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    ]
                )

                bottom_right = _np.array(
                    [
                        [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                        [v_outer_c_h + v_one_half, v_zero],
                        [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_one, v_zero],
                    ]
                )

            elif closure == "y_min":
                # set points:
                right = _np.array(
                    [
                        [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_one, -v_outer_c_h + v_one_half],
                        [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                        [v_one, v_outer_c_h + v_one_half],
                    ]
                )

                right_top = _np.array(
                    [
                        [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                        [v_one, v_outer_c_h + v_one_half],
                        [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                        [v_outer_c_h + v_one_half, v_one],
                    ]
                )

                top = _np.array(
                    [
                        [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                        [v_outer_c_h + v_one_half, v_one],
                        [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                        [-v_outer_c_h + v_one_half, v_one],
                    ]
                )

                bottom_left = _np.array(
                    [
                        [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_zero, -v_outer_c_h + v_one_half],
                        [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                        [v_zero, v_zero],
                    ]
                )

                left = _np.array(
                    [
                        [v_zero, -v_outer_c_h + v_one_half],
                        [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_zero, v_outer_c_h + v_one_half],
                        [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    ]
                )

                top_left = _np.array(
                    [
                        [v_zero, v_outer_c_h + v_one_half],
                        [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                        [-v_outer_c_h + v_one_half, v_one],
                        [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    ]
                )

                bottom = _np.array(
                    [
                        [v_one, v_zero],
                        [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                        [v_zero, v_zero],
                        [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    ]
                )

                bottom_right = _np.array(
                    [
                        [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                        [v_one, v_zero],
                        [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_one, -v_outer_c_h + v_one_half],
                    ]
                )

            elif closure == "y_max":
                # set points:
                right = _np.array(
                    [
                        [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_one, -v_outer_c_h + v_one_half],
                        [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                        [v_one, v_outer_c_h + v_one_half],
                    ]
                )

                right_top = _np.array(
                    [
                        [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                        [v_one, v_outer_c_h + v_one_half],
                        [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                        [v_one, v_one],
                    ]
                )

                top = _np.array(
                    [
                        [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                        [v_one, v_one],
                        [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                        [v_zero, v_one],
                    ]
                )

                bottom_left = _np.array(
                    [
                        [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_zero, -v_outer_c_h + v_one_half],
                        [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                        [-v_outer_c_h + v_one_half, v_zero],
                    ]
                )

                left = _np.array(
                    [
                        [v_zero, -v_outer_c_h + v_one_half],
                        [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_zero, v_outer_c_h + v_one_half],
                        [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    ]
                )

                top_left = _np.array(
                    [
                        [v_zero, v_outer_c_h + v_one_half],
                        [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                        [v_zero, v_one],
                        [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    ]
                )

                bottom = _np.array(
                    [
                        [v_outer_c_h + v_one_half, v_zero],
                        [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                        [-v_outer_c_h + v_one_half, v_zero],
                        [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    ]
                )

                bottom_right = _np.array(
                    [
                        [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                        [v_outer_c_h + v_one_half, v_zero],
                        [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                        [v_one, -v_outer_c_h + v_one_half],
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
                    _Bezier(degrees=[1, 1], control_points=control_points)
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
          the value must be between non-inclusive (0, 0.5).
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
        derivatives: list<list<splines>> / None
        """
        self._check_custom_parameter(
            contact_length, "contact length", self._CONTACT_LENGTH_BOUNDS
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
            )

        splines = []
        for i_derivative in range(n_derivatives + 1):
            if i_derivative == 0:
                v_zero = 0.0
                v_one_half = 0.5
                v_one = 1.0
                v_h_void = parameters[0, 0]
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
                    [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_one, -v_outer_c_h + v_one_half],
                    [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [v_one, v_outer_c_h + v_one_half],
                ]
            )

            right_top = _np.array(
                [
                    [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [v_one, v_outer_c_h + v_one_half],
                    [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [v_outer_c_h + v_one_half, v_one],
                ]
            )

            top = _np.array(
                [
                    [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [v_outer_c_h + v_one_half, v_one],
                    [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [-v_outer_c_h + v_one_half, v_one],
                ]
            )

            bottom_left = _np.array(
                [
                    [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_zero, -v_outer_c_h + v_one_half],
                    [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [-v_outer_c_h + v_one_half, v_zero],
                ]
            )

            left = _np.array(
                [
                    [v_zero, -v_outer_c_h + v_one_half],
                    [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_zero, v_outer_c_h + v_one_half],
                    [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                ]
            )

            top_left = _np.array(
                [
                    [v_zero, v_outer_c_h + v_one_half],
                    [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [-v_outer_c_h + v_one_half, v_one],
                    [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                ]
            )

            bottom = _np.array(
                [
                    [v_outer_c_h + v_one_half, v_zero],
                    [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [-v_outer_c_h + v_one_half, v_zero],
                    [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                ]
            )

            bottom_right = _np.array(
                [
                    [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [v_outer_c_h + v_one_half, v_zero],
                    [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_one, -v_outer_c_h + v_one_half],
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
                    _Bezier(degrees=[1, 1], control_points=control_points)
                )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        return (splines, derivatives)
