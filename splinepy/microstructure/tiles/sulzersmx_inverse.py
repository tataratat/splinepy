import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class SulzerSMXInverse(_TileBase):
    """Inverse of element used for Sulzer SMX static mixers. Using their nomenclature,
    the design parameters are: Nx=2 (number of cross bars), Np=1 (number of parallel
    cross bars)

    Ideas for parameters:
        - Thickness of one crossbar
        - Curvature of one crossbar:
            - just one curve
            - maybe multiple curves (wave-like)
        - If just one curve parameter, may set _n_info_per_eval_point to 2
    """

    _dim = 3
    _para_dim = 3
    _evaluation_points = _np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
        ]
    )
    _n_info_per_eval_point = 1
    _parameter_bounds = [[0.0, 0.5]] * 8
    _parameters_shape = (8, 1)
    _sensitivities_implemented = False

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        closure=None,
        **kwargs,  # noqa ARG002
    ):
        """Create a tile based on the parameters which describe the thickness of one
        crossbar.

        Parameters
        -----------
        parameters: np.ndarray
            The design parameters of the microtile. The parameters are as follows:
                1. Thickness of front crossbar
                2. Thickness of back crossbar
        parameter_sensitivities: np.ndarray
            The sensitivities of the parameters
        closure: str(optional)
            Parametric dimension which needs to be closed. Currently not implemented.

        Returns
        ----------
        microtile_list: list(splines)
        """

        if (
            parameters is not None
            and parameters.shape != self._parameters_shape
        ):
            raise ValueError(
                "Parameters are not in the correct shape. Needs "
                + str(self._parameters_shape)
            )

        # Set parameters to default values if nothing is given
        if parameters is None:
            self._logd("Setting cross bar thickness to default 0.2")
            parameters = 0.2 * _np.ones(self._parameters_shape)

        # TODO: Check if parameters are in allowed range

        self.check_params(parameters)

        if parameter_sensitivities is None:
            n_derivatives = 0
            derivatives = None
        else:
            self.check_param_derivatives(parameter_sensitivities)

            n_derivatives = parameter_sensitivities.shape[2]
            derivatives = []
            raise NotImplementedError(
                "Parameter sensitivities not implemented yet!"
            )

        if closure is not None:
            raise NotImplementedError("Closure not implemented!")

        def get_crossing_points(th_ne, th_se, th_sw, th_nw):
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

            return (
                crossing_s_x,
                crossing_s_y,
                crossing_w_x,
                crossing_w_y,
                crossing_e_x,
                crossing_e_y,
                crossing_n_x,
                crossing_n_y,
            )

        def create_triangle_yz_cps(v_zero, v_one, th_ne, th_se, th_sw, th_nw):
            # Get crossing points
            (
                crossing_s_x,
                crossing_s_y,
                crossing_w_x,
                crossing_w_y,
                crossing_e_x,
                crossing_e_y,
                crossing_n_x,
                crossing_n_y,
            ) = get_crossing_points(th_ne, th_se, th_sw, th_nw)
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

            return (
                crossing_s_x,
                crossing_s_y,
                crossing_w_x,
                crossing_w_y,
                crossing_e_x,
                crossing_e_y,
                crossing_n_x,
                crossing_n_y,
                connecting_interface_point_s_x,
                connecting_interface_point_w_y,
                connecting_interface_point_e_y,
                connecting_interface_point_n_x,
                triangle_s_middle,
                triangle_w_middle,
                triangle_e_middle,
                triangle_n_middle,
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
            )

        def create_bar_yz_cps(
            v_zero,
            v_one,
            th_ne,
            th_se,
            th_sw,
            th_nw,
            crossing_s_x,
            crossing_s_y,
            crossing_w_x,
            crossing_w_y,
            crossing_e_x,
            crossing_e_y,
            crossing_n_x,
            crossing_n_y,
        ):
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

            return (
                sw_bar_bottom,
                sw_bar_top,
                se_bar_bottom,
                se_bar_top,
                nw_bar_top,
                nw_bar_bottom,
                ne_bar_bottom,
                ne_bar_top,
            )

        e = _np.ones((4, 1))
        splines = []
        for i_derivative in range(n_derivatives + 1):
            if i_derivative == 0:
                (
                    th_swf,
                    th_sef,
                    th_nwf,
                    th_nef,
                    th_swb,
                    th_seb,
                    th_nwb,
                    th_neb,
                ) = parameters.flatten()

                v_zero = 0.0
                v_one_half = 0.5
                v_one = 1.0

            else:
                raise NotImplementedError("Derivatives not yet implemented")

            spline_list = []

            # Depth/x control points
            x_cps_front = _np.vstack((v_zero * e, v_one_half * e))
            x_cps_back = _np.vstack((v_one_half * e, v_one * e))

            # Create triangle pieces and auxiulary values
            triangle_pieces_front = create_triangle_yz_cps(
                v_zero, v_one, th_nef, th_sef, th_swf, th_nwf
            )
            triangle_pieces_back = create_triangle_yz_cps(
                v_zero, v_one, th_neb, th_seb, th_swb, th_nwb
            )

            # Go through the triangle pieces and create patches front to back
            for triangle_yz_cps_front, triangle_yz_cps_back in zip(
                triangle_pieces_front[16:], triangle_pieces_back[16:]
            ):
                # Cps are defined in xy-plane, but physical plane is zy-plane
                triangle_yz_cps_front = _np.fliplr(triangle_yz_cps_front)
                triangle_yz_cps_back = _np.fliplr(triangle_yz_cps_back)
                # Yz-control points at the middle yz-plane of the tile
                middle_yz_cps = (
                    triangle_yz_cps_front + triangle_yz_cps_back
                ) * 0.5
                # Go through front and back and create patches
                for x_cps, yz_cps in zip(
                    [x_cps_front, x_cps_back],
                    [
                        _np.vstack((triangle_yz_cps_front, middle_yz_cps)),
                        _np.vstack((middle_yz_cps, triangle_yz_cps_back)),
                    ],
                ):
                    control_points = _np.hstack((x_cps, yz_cps))
                    spline_list.append(
                        _Bezier(
                            degrees=[1, 1, 1], control_points=control_points
                        )
                    )

            # Add patches for the bars
            # Get bar pieces
            bar_pieces_front = create_bar_yz_cps(
                v_zero,
                v_one,
                th_nef,
                th_sef,
                th_swf,
                th_nwf,
                *triangle_pieces_front[:8],
            )
            bar_pieces_back = create_bar_yz_cps(
                v_zero,
                v_one,
                th_neb,
                th_seb,
                th_swb,
                th_nwb,
                *triangle_pieces_back[:8],
            )

            # Define whether to put patches at front or at back. The list corresponds
            # to the return values in create_bar_yz_cps
            bar_at_front_list = [
                False,
                False,
                True,
                True,
                True,
                True,
                False,
                False,
            ]

            for bar_yz_front, bar_yz_back, bar_at_front in zip(
                bar_pieces_front, bar_pieces_back, bar_at_front_list
            ):
                # Control points are defined in xy-plane, but physical points should
                # lie in zy-plane
                bar_yz_front = _np.fliplr(bar_yz_front)
                bar_yz_back = _np.fliplr(bar_yz_back)
                middle_yz_cps = (bar_yz_front + bar_yz_back) * 0.5
                if bar_at_front:
                    x_cps = x_cps_front
                    yz_cps = _np.vstack((bar_yz_front, middle_yz_cps))
                else:
                    x_cps = x_cps_back
                    yz_cps = _np.vstack((middle_yz_cps, bar_yz_back))
                control_points = _np.hstack((x_cps, yz_cps))
                spline_list.append(
                    _Bezier(degrees=[1, 1, 1], control_points=control_points)
                )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        return (splines, derivatives)
