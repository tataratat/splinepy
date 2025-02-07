import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class SMX2DInverse(_TileBase):
    """Side view of the inverse of Sulzer SMX tile.

    Parameters:
        - Thickness of bar at NE-corner
        - Thickness of bar at SE-corner
        - Thickness of bar at SW-corner
        - Thickness of bar at NW-corner

    Ideas for 3D counterpart:
        - Two _n_info_per_eval_point: corresponding to bar width (problem: self
        intersections of bars possible)
    """

    _dim = 2
    _para_dim = 2
    _evaluation_points = _np.array(
        [[1.0, 1.0], [1.0, 0.0], [0.0, 0.0], [0.0, 1.0]]
    )
    _n_info_per_eval_point = 1
    _parameter_bounds = [[0.0, 0.5], [0.0, 0.5], [0.0, 0.5], [0.0, 0.5]]
    _parameters_shape = (4, 1)
    _sensitivities_implemented = False
    _default_parameter_value = 0.2

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        closure=None,
        **kwargs,  # noqa ARG002
    ):
        """Create a tile based on the parameters which describe the thicknesses of
        the crossbar at the corners

        Parameters
        ---------
        parameters: _np.ndarray
            The design parameters of the microtile. The parameters are as follows:
                1. Thickness of bar at NE-corner
                2. Thickness of bar at SE-corner
                3. Thickness of bar at SW-corner
                4. Thickness of bar at NW-corner
        parameter_sensitivities: _np.ndarray
            The sensitivities of the parameters
        closure: str(optional)
            Parametric dimension which needs to be closed. Currently not implemented.

        Returns
        ------------
        microtile_list: list(splines)
        """

        parameters, n_derivatives, derivatives = self._process_input(
            parameters=parameters,
            parameter_sensitivities=parameter_sensitivities,
        )

        bar_thickness = 0.2

        splines = []
        for i_derivative in range(n_derivatives + 1):
            if i_derivative == 0:
                v_zero = 0.0
                v_one_half = 0.5
                v_one = 1.0

                th_ne, th_se, th_sw, th_nw = parameters.ravel()
            else:
                raise NotImplementedError(
                    "Derivatives not implemented for this tile"
                )

            # Define crossing points of bars, analytically computed beforehand
            crossing_s_x = (
                -th_ne * th_se * th_sw
                + th_ne * th_sw
                + th_nw * th_se * th_sw
                - th_nw * th_se
                - th_nw * th_sw
                + th_nw
                + th_se
                - 1
            ) / (
                -th_ne * th_se
                + th_ne
                - th_nw * th_sw
                + th_nw
                + th_se
                + th_sw
                - 2
            )
            crossing_s_y = (
                -th_ne * th_nw * th_se
                - th_ne * th_nw * th_sw
                + th_ne * th_nw
                + th_ne * th_se
                + th_ne * th_sw
                - th_ne
                + th_nw * th_se
                + th_nw * th_sw
                - th_nw
                - th_se
                - th_sw
                + 1
            ) / (
                th_ne * th_se
                - th_ne
                + th_nw * th_sw
                - th_nw
                - th_se
                - th_sw
                + 2
            )
            crossing_w_x = (
                th_ne * th_nw * th_se
                - th_ne * th_nw
                + th_ne * th_se * th_sw
                - th_ne * th_se
                - th_ne * th_sw
                + th_ne
                - th_nw * th_se
                + th_nw
                - th_se * th_sw
                + th_se
                + th_sw
                - 1
            ) / (
                -th_ne * th_nw
                + th_ne
                + th_nw
                - th_se * th_sw
                + th_se
                + th_sw
                - 2
            )
            crossing_w_y = (
                th_ne * th_nw * th_sw
                - th_ne * th_sw
                - th_nw * th_se * th_sw
                + th_nw * th_se
                - th_nw
                + th_se * th_sw
                - th_se
                + 1
            ) / (
                th_ne * th_nw
                - th_ne
                - th_nw
                + th_se * th_sw
                - th_se
                - th_sw
                + 2
            )
            crossing_e_x = (
                -th_ne * th_nw * th_sw
                + th_ne * th_sw
                - th_nw * th_se * th_sw
                + th_nw * th_se
                + th_nw * th_sw
                - 1
            ) / (
                -th_ne * th_nw
                + th_ne
                + th_nw
                - th_se * th_sw
                + th_se
                + th_sw
                - 2
            )
            crossing_e_y = (
                th_ne * th_nw * th_se
                - th_ne * th_se * th_sw
                + th_ne * th_sw
                - th_ne
                - th_nw * th_se
                + th_se * th_sw
                - th_sw
                + 1
            ) / (
                th_ne * th_nw
                - th_ne
                - th_nw
                + th_se * th_sw
                - th_se
                - th_sw
                + 2
            )
            crossing_n_x = (
                -th_ne * th_nw * th_se
                + th_ne * th_nw * th_sw
                - th_ne * th_sw
                + th_ne
                + th_nw * th_se
                - th_nw * th_sw
                + th_sw
                - 1
            ) / (
                -th_ne * th_se
                + th_ne
                - th_nw * th_sw
                + th_nw
                + th_se
                + th_sw
                - 2
            )
            crossing_n_y = (
                th_ne * th_se * th_sw
                - th_ne * th_sw
                + th_nw * th_se * th_sw
                - th_nw * th_se
                - th_se * th_sw
                + 1
            ) / (
                th_ne * th_se
                - th_ne
                + th_nw * th_sw
                - th_nw
                - th_se
                - th_sw
                + 2
            )

            spline_list = []

            # Define control points
            sw_bar_bottom = _np.array(
                [
                    [th_sw, v_zero],
                    [crossing_s_x, crossing_s_y],
                    [th_sw / 2, th_sw / 2],
                    [
                        (crossing_w_x + crossing_s_x) / 2,
                        (crossing_w_y + crossing_s_y) / 2,
                    ],
                ]
            )

            sw_bar_top = _np.array(
                [
                    [th_sw / 2, th_sw / 2],
                    [
                        (crossing_w_x + crossing_s_x) / 2,
                        (crossing_w_y + crossing_s_y) / 2,
                    ],
                    [v_zero, th_sw],
                    [crossing_w_x, crossing_w_y],
                ]
            )

            se_bar_bottom = _np.array(
                [
                    [crossing_s_x, crossing_s_y],
                    [v_one - th_se, v_zero],
                    [
                        (crossing_s_x + crossing_e_x) / 2,
                        (crossing_s_y + crossing_e_y) / 2,
                    ],
                    [v_one - th_se / 2, th_se / 2],
                ]
            )

            se_bar_top = _np.array(
                [
                    [
                        (crossing_s_x + crossing_e_x) / 2,
                        (crossing_s_y + crossing_e_y) / 2,
                    ],
                    [v_one - th_se / 2, th_se / 2],
                    [crossing_e_x, crossing_e_y],
                    [v_one, th_se],
                ]
            )

            nw_bar_bottom = _np.array(
                [
                    [v_zero, v_one - th_nw],
                    [crossing_w_x, crossing_w_y],
                    [th_nw / 2, v_one - th_nw / 2],
                    [
                        (crossing_w_x + crossing_n_x) / 2,
                        (crossing_w_y + crossing_n_y) / 2,
                    ],
                ]
            )

            nw_bar_top = _np.array(
                [
                    [th_nw / 2, v_one - th_nw / 2],
                    [
                        (crossing_w_x + crossing_n_x) / 2,
                        (crossing_w_y + crossing_n_y) / 2,
                    ],
                    [th_nw, v_one],
                    [crossing_n_x, crossing_n_y],
                ]
            )

            ne_bar_bottom = _np.array(
                [
                    [crossing_e_x, crossing_e_y],
                    [v_one, v_one - th_ne],
                    [
                        (crossing_e_x + crossing_n_x) / 2,
                        (crossing_e_y + crossing_n_y) / 2,
                    ],
                    [v_one - th_ne / 2, v_one - th_ne / 2],
                ]
            )

            ne_bar_top = _np.array(
                [
                    [
                        (crossing_e_x + crossing_n_x) / 2,
                        (crossing_e_y + crossing_n_y) / 2,
                    ],
                    [v_one - th_ne / 2, v_one - th_ne / 2],
                    [crossing_n_x, crossing_n_y],
                    [v_one - th_ne, v_one],
                ]
            )

            # Define control points for middle part
            middle_s = _np.array(
                [
                    [crossing_s_x, crossing_s_y],
                    [
                        (crossing_s_x + crossing_e_x) / 2,
                        (crossing_s_y + crossing_e_y) / 2,
                    ],
                    [
                        (crossing_s_x + crossing_w_x) / 2,
                        (crossing_s_y + crossing_w_y) / 2,
                    ],
                    [v_one_half, v_one_half],
                ]
            )

            middle_w = _np.array(
                [
                    [
                        (crossing_s_x + crossing_w_x) / 2,
                        (crossing_s_y + crossing_w_y) / 2,
                    ],
                    [v_one_half, v_one_half],
                    [crossing_w_x, crossing_w_y],
                    [
                        (crossing_n_x + crossing_w_x) / 2,
                        (crossing_n_y + crossing_w_y) / 2,
                    ],
                ]
            )

            middle_e = _np.array(
                [
                    [
                        (crossing_s_x + crossing_e_x) / 2,
                        (crossing_s_y + crossing_e_y) / 2,
                    ],
                    [crossing_e_x, crossing_e_y],
                    [v_one_half, v_one_half],
                    [
                        (crossing_n_x + crossing_e_x) / 2,
                        (crossing_n_y + crossing_e_y) / 2,
                    ],
                ]
            )

            middle_n = _np.array(
                [
                    [v_one_half, v_one_half],
                    [
                        (crossing_n_x + crossing_e_x) / 2,
                        (crossing_n_y + crossing_e_y) / 2,
                    ],
                    [
                        (crossing_n_x + crossing_w_x) / 2,
                        (crossing_n_y + crossing_w_y) / 2,
                    ],
                    [crossing_n_x, crossing_n_y],
                ]
            )

            for control_points in [
                sw_bar_bottom,
                sw_bar_top,
                se_bar_bottom,
                se_bar_top,
                nw_bar_top,
                nw_bar_bottom,
                ne_bar_bottom,
                ne_bar_top,
                middle_s,
                middle_w,
                middle_e,
                middle_n,
            ]:
                spline_list.append(
                    _Bezier(degrees=[1, 1], control_points=control_points)
                )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        return (splines, derivatives)
