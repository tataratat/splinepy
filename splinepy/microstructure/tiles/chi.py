import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class Chi(_TileBase):
    def __init__(self):
        """Chi"""
        self._dim = 2
        self._para_dim = 1
        self._evaluation_points = _np.array([[0.5, 0.5]])
        self._n_info_per_eval_point = 1

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        **kwargs,  # noqa ARG002
    ):
        """Create a microtile based on the parameter that describes a
        drilled cross around the center

        Parameters
        ----------
        parameters : np.array
          only first entry is used, defines the angle of the branches in
          radian, where zero results in a cross
        parameter_sensitivities: np.ndarray
          correlates with the drill angle

        Returns
        -------
        microtile_list : list(splines)
        """

        # set to default if nothing is given
        if parameters is None:
            self._logd(
                "Tile request is not parametrized, setting default Pi/8"
            )
            parameters = _np.array([[_np.pi / 8]])
        else:
            if not (
                _np.all(parameters >= -_np.pi * 0.5)
                and _np.all(parameters < _np.pi * 0.5)
            ):
                raise ValueError("The parameter must be in -Pi/2 and Pi/2")
            pass
        self.check_params(parameters)

        # Check if user requests derivative splines
        if self.check_param_derivatives(parameter_sensitivities):
            n_derivatives = parameter_sensitivities.shape[2]
            derivatives = []
        else:
            n_derivatives = 0
            derivatives = None

        splines = []
        for i_derivative in range(n_derivatives + 1):
            # Constant auxiliary values
            if i_derivative == 0:
                alpha = parameters[0, 0] + _np.pi / 4
                v_one_half = 0.5
                v_zero = 0.0
                r = _np.sqrt(0.125)
                s = r * _np.sin(alpha)
                c = r * _np.cos(alpha)
            else:
                alpha = parameter_sensitivities[0, 0, i_derivative - 1]
                v_one_half = 0.0
                v_zero = 0.0
                r = _np.sqrt(0.125)
                s = r * _np.cos(alpha)
                c = -r * _np.sin(alpha)

            # Init return value
            spline_list = []
            # 1
            spline_list.append(
                _Bezier(
                    degrees=[2],
                    control_points=_np.array(
                        [
                            [v_zero, v_zero],
                            [-s, -c],
                            [-v_one_half, -v_one_half],
                        ]
                    )
                    + v_one_half,
                )
            )
            # 2
            spline_list.append(
                _Bezier(
                    degrees=[2],
                    control_points=_np.array(
                        [
                            [v_zero, v_zero],
                            [c, -s],
                            [v_one_half, -v_one_half],
                        ]
                    )
                    + v_one_half,
                )
            )
            # 3
            spline_list.append(
                _Bezier(
                    degrees=[2],
                    control_points=_np.array(
                        [
                            [v_zero, v_zero],
                            [s, c],
                            [v_one_half, v_one_half],
                        ]
                    )
                    + v_one_half,
                )
            )
            # 4
            spline_list.append(
                _Bezier(
                    degrees=[2],
                    control_points=_np.array(
                        [
                            [v_zero, v_zero],
                            [-c, s],
                            [-v_one_half, v_one_half],
                        ]
                    )
                    + v_one_half,
                )
            )
            # Pass to output
            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        # Return results
        return (splines, derivatives)
