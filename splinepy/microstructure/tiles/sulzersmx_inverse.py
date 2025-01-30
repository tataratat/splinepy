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
    _evaluation_points = _np.array([[0.5, 0.5, 0.5], [0.5, 0.5, 0.5]])
    _n_info_per_eval_point = 1
    _parameter_bounds = [[0.0, 0.5], [0.0, 0.5]]
    _parameters_shape = (2, 1)
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
            parameters = _np.array([[0.1], [0.1]])
            # TODO: Only simulatable if both parameters are the same

        # Check if parameters are in allowed range
        if _np.sum(
            (parameters.flatten() > _np.array(self._parameter_bounds)[:, 0])
            & (parameters.flatten() < _np.array(self._parameter_bounds)[:, 1])
        ) < _np.prod(self._parameters_shape):
            raise ValueError("Given parameters are not within allowed range!")

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

        # Helper array for defining x-component of control points
        e = _np.ones((4, 1))

        splines = []
        for i_derivative in range(n_derivatives + 1):
            if i_derivative == 0:
                [thickness_front, thickness_back] = parameters.flatten()

                v_zero = 0.0
                v_one_half = 0.5
                v_one = 1.0
            else:
                raise NotImplementedError("Derivatives not yet implemented")
                [thickness_front, thickness_back] = (
                    parameter_sensitivities.flatten()
                )

                v_zero, v_one = [0.0] * 2

            # Auxiliary values
            triangle_height = 0.5 * (v_one - thickness_front - thickness_back)

            spline_list = []

            # Define control points of crossbars
            xs_front = _np.vstack((v_zero * e, v_one_half * e))
            xs_back = _np.vstack((v_one_half * e, v_one * e))

            # Define the y- and z-coordinates of the triangle on the bottom
            triangle_bottom_sw = _np.array(
                [
                    [v_zero, thickness_front],
                    [v_zero, v_one_half],
                    [
                        triangle_height * 0.5,
                        v_one_half - triangle_height * 0.5,
                    ],
                    [triangle_height / 3, v_one_half],
                ]
            )

            triangle_bottom_se = _np.array(
                [
                    [v_zero, v_one_half],
                    [v_zero, v_one - thickness_back],
                    [triangle_height / 3, v_one_half],
                    [
                        triangle_height * 0.5,
                        v_one_half + triangle_height * 0.5,
                    ],
                ]
            )
            triangle_bottom_n = _np.array(
                [
                    [
                        triangle_height * 0.5,
                        v_one_half - triangle_height * 0.5,
                    ],
                    [triangle_height / 3, v_one_half],
                    [triangle_height, v_one_half],
                    [
                        triangle_height * 0.5,
                        v_one_half + triangle_height * 0.5,
                    ],
                ]
            )

            triangle_right_s = _np.array(
                [
                    [
                        v_one_half - triangle_height * 0.5,
                        v_one - triangle_height * 0.5,
                    ],
                    [thickness_back, v_one],
                    [v_one_half, v_one - triangle_height / 3],
                    [v_one_half, v_one],
                ]
            )

            triangle_right_n = _np.array(
                [
                    [v_one_half, v_one - triangle_height / 3],
                    [v_one_half, v_one],
                    [
                        v_one_half + triangle_height * 0.5,
                        v_one - triangle_height * 0.5,
                    ],
                    [v_one - thickness_front, v_one],
                ]
            )

            triangle_right_w = _np.array(
                [
                    [v_one_half, v_one - triangle_height],
                    [
                        v_one_half - triangle_height * 0.5,
                        v_one - triangle_height * 0.5,
                    ],
                    [
                        v_one_half + triangle_height * 0.5,
                        v_one - triangle_height * 0.5,
                    ],
                    [v_one_half, v_one - triangle_height / 3],
                ]
            )

            triangle_top_nw, triangle_top_ne, triangle_top_s = (
                _np.hstack(
                    (
                        v_one - cps[:, 0].reshape(-1, 1),
                        cps[:, 1].reshape(-1, 1),
                    )
                )
                for cps in [
                    triangle_bottom_sw,
                    triangle_bottom_se,
                    triangle_bottom_n,
                ]
            )

            triangle_left_s, triangle_left_n, triangle_left_e = (
                _np.hstack(
                    (
                        cps[:, 0].reshape(-1, 1),
                        v_one - cps[:, 1].reshape(-1, 1),
                    )
                )
                for cps in [
                    triangle_right_s,
                    triangle_right_n,
                    triangle_right_w,
                ]
            )

            # Create patches which have front as well as back
            for yz_cps in [
                triangle_bottom_sw,
                triangle_bottom_se,
                triangle_bottom_n,
                triangle_top_nw,
                triangle_top_ne,
                triangle_top_s,
                triangle_right_s,
                triangle_right_n,
                triangle_right_w,
                triangle_left_s,
                triangle_left_n,
                triangle_left_e,
            ]:
                for xs in [xs_front, xs_back]:
                    control_points = _np.hstack((xs, _np.tile(yz_cps, [2, 1])))
                    spline_list.append(
                        _Bezier(
                            degrees=[1, 1, 1], control_points=control_points
                        )
                    )

            # Create the bars
            front_bottom_bar_s = _np.array(
                [
                    [
                        triangle_height * 0.5,
                        v_one - thickness_back - triangle_height * 0.5,
                    ],
                    [v_zero, v_one - thickness_back],
                    [
                        thickness_back + triangle_height * 0.5,
                        v_one - triangle_height * 0.5,
                    ],
                    [thickness_back, v_one],
                ]
            )

            front_bottom_bar_n = _np.array(
                [
                    [triangle_height, v_one_half],
                    [
                        triangle_height * 0.5,
                        v_one - thickness_back - triangle_height * 0.5,
                    ],
                    [v_one_half, v_one - triangle_height],
                    [
                        thickness_back + triangle_height * 0.5,
                        v_one - triangle_height * 0.5,
                    ],
                ]
            )

            # Helper functions for control point translation
            def identity_function(cps):
                return cps

            def opposite_corner_function(cps):
                return v_one - cps

            # Append bar patches to spline list
            for yz_cps in [front_bottom_bar_s, front_bottom_bar_n]:
                # 1: create front bottom, 2: create front top
                for translation_function in [
                    identity_function,
                    opposite_corner_function,
                ]:
                    yz_cps_translated = translation_function(yz_cps)
                    control_points = _np.hstack(
                        (xs_front, _np.tile(yz_cps_translated, [2, 1]))
                    )
                    spline_list.append(
                        _Bezier(
                            degrees=[1, 1, 1], control_points=control_points
                        )
                    )

            # Define bar control points
            back_bottom_bar_s = _np.array(
                [
                    [v_zero, thickness_front],
                    [
                        triangle_height * 0.5,
                        v_one_half - triangle_height * 0.5,
                    ],
                    [thickness_front, v_zero],
                    [
                        v_one_half - triangle_height * 0.5,
                        triangle_height * 0.5,
                    ],
                ]
            )

            back_bottom_bar_n = _np.array(
                [
                    [
                        triangle_height * 0.5,
                        v_one_half - triangle_height * 0.5,
                    ],
                    [triangle_height, v_one_half],
                    [
                        v_one_half - triangle_height * 0.5,
                        triangle_height * 0.5,
                    ],
                    [v_one_half, triangle_height],
                ]
            )

            # Add back bar patches to spline list
            for yz_cps in [back_bottom_bar_s, back_bottom_bar_n]:
                # 1: create back bottom, 2: create back top bar
                for translation_function in [
                    identity_function,
                    opposite_corner_function,
                ]:
                    yz_cps_translated = translation_function(yz_cps)
                    control_points = _np.hstack(
                        (xs_back, _np.tile(yz_cps_translated, [2, 1]))
                    )
                    spline_list.append(
                        _Bezier(
                            degrees=[1, 1, 1], control_points=control_points
                        )
                    )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        return (splines, derivatives)
