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
        closure=None,  # noqa ARG002
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

        splines = []
        for i_derivative in range(n_derivatives + 1):
            if i_derivative == 0:
                v_zero = 0.0
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
            # Bars
            sw_bar_bottom = _np.array(
                [
                    [th_sw, v_zero],
                    [(th_sw + crossing_s_x) / 2, crossing_s_y / 2],
                    [v_zero, th_sw],
                    [crossing_w_x / 2, (th_sw + crossing_w_y) / 2],
                ]
            )

            sw_bar_top = _np.array(
                [
                    [(th_sw + crossing_s_x) / 2, crossing_s_y / 2],
                    [crossing_s_x, crossing_s_y],
                    [crossing_w_x / 2, (th_sw + crossing_w_y) / 2],
                    [crossing_w_x, crossing_w_y],
                ]
            )

            se_bar_bottom = _np.array(
                [
                    [v_one - th_se, v_zero],
                    [v_one, th_se],
                    [(v_one - th_se + crossing_s_x) / 2, crossing_s_y / 2],
                    [(crossing_e_x + v_one) / 2, (th_se + crossing_e_y) / 2],
                ]
            )

            se_bar_top = _np.array(
                [
                    [(v_one - th_se + crossing_s_x) / 2, crossing_s_y / 2],
                    [(crossing_e_x + v_one) / 2, (th_se + crossing_e_y) / 2],
                    [crossing_s_x, crossing_s_y],
                    [crossing_e_x, crossing_e_y],
                ]
            )

            nw_bar_bottom = _np.array(
                [
                    [crossing_w_x, crossing_w_y],
                    [crossing_n_x, crossing_n_y],
                    [crossing_w_x / 2, (crossing_w_y + v_one - th_nw) / 2],
                    [
                        (th_nw + crossing_n_x) / 2,
                        (v_one + crossing_n_y) / 2,
                    ],
                ]
            )

            nw_bar_top = _np.array(
                [
                    [crossing_w_x / 2, (crossing_w_y + v_one - th_nw) / 2],
                    [
                        (th_nw + crossing_n_x) / 2,
                        (v_one + crossing_n_y) / 2,
                    ],
                    [v_zero, v_one - th_nw],
                    [th_nw, v_one],
                ]
            )

            ne_bar_bottom = _np.array(
                [
                    [crossing_n_x, crossing_n_y],
                    [crossing_e_x, crossing_e_y],
                    [
                        (crossing_n_x + v_one - th_ne) / 2,
                        (crossing_n_y + v_one) / 2,
                    ],
                    [
                        (crossing_e_x + v_one) / 2,
                        (crossing_e_y + v_one - th_ne) / 2,
                    ],
                ]
            )

            ne_bar_top = _np.array(
                [
                    [
                        (crossing_n_x + v_one - th_ne) / 2,
                        (crossing_n_y + v_one) / 2,
                    ],
                    [
                        (crossing_e_x + v_one) / 2,
                        (crossing_e_y + v_one - th_ne) / 2,
                    ],
                    [v_one - th_ne, v_one],
                    [v_one, v_one - th_ne],
                ]
            )

            # Define the patches for the triangle
            # Define the middle points of the triangles, not as the triangle's centroid
            # but rather as a point connecting the triangle's innermost point to the
            # connecting point at the tile's interface

            # Connecting points at the interface
            connecting_interface_point_s_x = (
                th_sw + (v_one - th_sw - th_se) / 2
            )
            connecting_interface_point_w_y = (
                th_sw + (v_one - th_sw - th_nw) / 2
            )
            connecting_interface_point_e_y = (
                th_se + (v_one - th_se - th_ne) / 2
            )
            connecting_interface_point_n_x = (
                th_nw + (v_one - th_nw - th_ne) / 2
            )

            # Triangle middle points
            triangle_s_middle = [
                (2 * connecting_interface_point_s_x + crossing_s_x) / 3,
                crossing_s_y / 3,
            ]
            triangle_w_middle = [
                crossing_w_x / 3,
                (2 * connecting_interface_point_w_y + crossing_w_y) / 3,
            ]
            triangle_e_middle = [
                (2 * v_one + crossing_e_x) / 3,
                (2 * connecting_interface_point_e_y + crossing_e_y) / 3,
            ]
            triangle_n_middle = [
                (2 * connecting_interface_point_n_x + crossing_n_x) / 3,
                (2 * v_one + crossing_n_y) / 3,
            ]

            # South triangle
            triangle_s_sw = _np.array(
                [
                    [th_sw, v_zero],
                    [connecting_interface_point_s_x, v_zero],
                    [(th_sw + crossing_s_x) / 2, crossing_s_y / 2],
                    triangle_s_middle,
                ]
            )

            triangle_s_se = _np.array(
                [
                    [connecting_interface_point_s_x, v_zero],
                    [v_one - th_se, v_zero],
                    triangle_s_middle,
                    [(crossing_s_x + v_one - th_se) / 2, crossing_s_y / 2],
                ]
            )

            triangle_s_n = _np.array(
                [
                    [(th_sw + crossing_s_x) / 2, crossing_s_y / 2],
                    triangle_s_middle,
                    [crossing_s_x, crossing_s_y],
                    [(crossing_s_x + v_one - th_se) / 2, crossing_s_y / 2],
                ]
            )

            # West triangle
            triangle_w_s = _np.array(
                [
                    [v_zero, th_sw],
                    [crossing_w_x / 2, (th_sw + crossing_w_y) / 2],
                    [v_zero, connecting_interface_point_w_y],
                    triangle_w_middle,
                ]
            )

            triangle_w_n = _np.array(
                [
                    [v_zero, connecting_interface_point_w_y],
                    triangle_w_middle,
                    [v_zero, v_one - th_nw],
                    [crossing_w_x / 2, (v_one - th_nw + crossing_w_y) / 2],
                ]
            )

            triangle_w_e = _np.array(
                [
                    [crossing_w_x / 2, (th_sw + crossing_w_y) / 2],
                    [crossing_w_x, crossing_w_y],
                    triangle_w_middle,
                    [crossing_w_x / 2, (v_one - th_nw + crossing_w_y) / 2],
                ]
            )

            # East triangle
            triangle_e_s = _np.array(
                [
                    [(v_one + crossing_e_x) / 2, (th_se + crossing_e_y) / 2],
                    [v_one, th_se],
                    triangle_e_middle,
                    [v_one, connecting_interface_point_e_y],
                ]
            )

            triangle_e_n = _np.array(
                [
                    triangle_e_middle,
                    [v_one, connecting_interface_point_e_y],
                    [
                        (v_one + crossing_e_x) / 2,
                        (crossing_e_y + v_one - th_ne) / 2,
                    ],
                    [v_one, v_one - th_ne],
                ]
            )

            triangle_e_w = _np.array(
                [
                    [crossing_e_x, crossing_e_y],
                    [(crossing_e_x + v_one) / 2, (crossing_e_y + th_se) / 2],
                    [
                        (v_one + crossing_e_x) / 2,
                        (crossing_e_y + v_one - th_ne) / 2,
                    ],
                    triangle_e_middle,
                ]
            )

            # North triangle
            triangle_n_s = _np.array(
                [
                    [crossing_n_x, crossing_n_y],
                    [
                        (crossing_n_x + v_one - th_ne) / 2,
                        (crossing_n_y + v_one) / 2,
                    ],
                    [(crossing_n_x + th_nw) / 2, (crossing_n_y + v_one) / 2],
                    triangle_n_middle,
                ]
            )

            triangle_n_nw = _np.array(
                [
                    [(crossing_n_x + th_nw) / 2, (crossing_n_y + v_one) / 2],
                    triangle_n_middle,
                    [th_nw, v_one],
                    [connecting_interface_point_n_x, v_one],
                ]
            )

            triangle_n_ne = _np.array(
                [
                    [
                        (crossing_n_x + v_one - th_ne) / 2,
                        (crossing_n_y + v_one) / 2,
                    ],
                    [v_one - th_ne, v_one],
                    triangle_n_middle,
                    [connecting_interface_point_n_x, v_one],
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
                triangle_s_sw,
                triangle_s_se,
                triangle_s_n,
                triangle_w_s,
                triangle_w_n,
                triangle_w_e,
                triangle_e_s,
                triangle_e_n,
                triangle_e_w,
                triangle_n_s,
                triangle_n_nw,
                triangle_n_ne,
            ]:
                spline_list.append(
                    _Bezier(degrees=[1, 1], control_points=control_points)
                )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        return (splines, derivatives)
