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
    _evaluation_points = _np.array([[0.25, 0.5, 0.5], [0.75, 0.5, 0.5]])
    _n_info_per_eval_point = 1
    _parameter_bounds = [[0.0, 1.0], [0.0, 1.0]]
    _parameters_shape = (2, 1)

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
            parameters = _np.array([[0.2], [0.2]])

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

            spline_list = []

            # Define control points of crossbars
            front_bar_bottom = _np.array(
                [
                    [v_zero, v_zero, thickness_front],
                    [v_one_half, v_zero, thickness_front],
                    [
                        v_zero,
                        (v_one - thickness_front) / 2,
                        (v_one + thickness_front) / 2,
                    ],
                    [
                        v_one_half,
                        (v_one - thickness_front) / 2,
                        (v_one + thickness_front) / 2,
                    ],
                    [v_zero, v_zero, v_one],
                    [v_one_half, v_zero, v_one],
                    [v_zero, v_one - thickness_front, v_one],
                    [v_one_half, v_one - thickness_front, v_one],
                ]
            )

            front_bar_top = _np.array(
                [
                    [v_zero, thickness_front, v_zero],
                    [v_one_half, thickness_front, v_zero],
                    [v_zero, v_one, v_zero],
                    [v_one_half, v_one, v_zero],
                    [
                        v_zero,
                        (v_one + thickness_front) / 2,
                        (v_one - thickness_front) / 2,
                    ],
                    [
                        v_one_half,
                        (v_one + thickness_front) / 2,
                        (v_one - thickness_front) / 2,
                    ],
                    [v_zero, v_one, v_one - thickness_front],
                    [v_one_half, v_one, v_one - thickness_front],
                ]
            )

            back_bar_bottom = _np.array(
                [
                    [v_one_half, v_zero, v_zero],
                    [v_one, v_zero, v_zero],
                    [v_one_half, v_one - thickness_back, v_zero],
                    [v_one, v_one - thickness_back, v_zero],
                    [v_one_half, v_zero, v_one - thickness_back],
                    [v_one, v_zero, v_one - thickness_back],
                    [
                        v_one_half,
                        (v_one - thickness_back) / 2,
                        (v_one - thickness_back) / 2,
                    ],
                    [
                        v_one,
                        (v_one - thickness_back) / 2,
                        (v_one - thickness_back) / 2,
                    ],
                ]
            )

            back_bar_top = _np.array(
                [
                    [
                        v_one_half,
                        (v_one + thickness_back) / 2,
                        (v_one + thickness_back) / 2,
                    ],
                    [
                        v_one,
                        (v_one + thickness_back) / 2,
                        (v_one + thickness_back) / 2,
                    ],
                    [v_one_half, v_one, thickness_back],
                    [v_one, v_one, thickness_back],
                    [v_one_half, thickness_back, v_one],
                    [v_one, thickness_back, v_one],
                    [v_one_half, v_one, v_one],
                    [v_one, v_one, v_one],
                ]
            )

            # Append all the patches into list
            for control_points in [
                front_bar_bottom,
                front_bar_top,
                back_bar_bottom,
                back_bar_top,
            ]:
                spline_list.append(
                    _Bezier(degrees=[1, 1, 1], control_points=control_points)
                )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        return (splines, derivatives)
