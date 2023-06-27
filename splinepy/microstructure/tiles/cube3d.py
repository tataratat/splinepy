import numpy as np

from splinepy.bezier import Bezier
from splinepy.microstructure.tiles.tilebase import TileBase


class Cube3D(TileBase):
    def __init__(self):
        """Simple tile - looks like a nut"""
        self._dim = 3
        self._evaluation_points = np.array(
            [
                [0.0, 0.5, 0.5],
                [1.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 1.0, 0.5],
                [0.5, 0.5, 0.0],
                [0.5, 0.5, 1.0],
                [0.5, 0.5, 0.5],
            ]
        )
        self._n_info_per_eval_point = 1

    def create_tile(
        self, parameters=None, parameter_sensitivities=None, **kwargs
    ):
        """Create a microtile based on the parameters that describe the strut
        thicknesses.


        Parameters
        ----------
        parameters :
          specifies the hole-thickness on the tile interfaces. Last parameter
          is used for center dimensions
        parameter_sensitivities: np.ndarray
          Describes the parameter sensitivities with respect to some design
          variable. In case the design variables directly apply to the
          parameter itself, they evaluate to delta_ij

        Returns
        -------
        microtile_list : list(splines)
        """

        if parameters is None:
            self._logd("Setting parameters to default value (0.2)")
            parameters = (
                np.ones(
                    (len(self._evaluation_points), self._n_info_per_eval_point)
                ).reshape(-1, 1)
                * 0.2
            )

        self.check_params(parameters)

        if not (
            np.all(parameters[:, :2] > 0.0) and np.all(parameters[:, :2] < 0.5)
        ):
            raise ValueError("The wall thickness must be in (0.0 and 0.5)")

        if not np.all(
            (parameters[:, 2:] > np.deg2rad(-30))
            and (parameters[:, 2:] < np.deg2rad(30))
        ):
            raise ValueError("Rotation is only allowed between +-30 deg")

        if self.check_param_derivatives(parameter_sensitivities):
            n_derivatives = parameter_sensitivities.shape[2]
        else:
            n_derivatives = 0

        derivatives = []
        splines = []

        for i_derivative in range(n_derivatives + 1):
            spline_list = []
            # Constant auxiliary values
            if i_derivative == 0:
                v_zero = 0.0
                v_one = 1.0
                [
                    x_min,
                    x_max,
                    y_min,
                    y_max,
                    z_min,
                    z_max,
                    center,
                ] = parameters[:, 0].flatten()
            else:
                v_zero = 0.0
                v_one = 0.0
                [
                    x_min,
                    x_max,
                    y_min,
                    y_max,
                    z_min,
                    z_max,
                    center,
                ] = parameter_sensitivities[:, 0, i_derivative - 1]

            # x_max_z_max
            spline_list.append(
                Bezier(
                    degrees=[1, 1, 1],
                    control_points=np.array(
                        [
                            [v_one - z_max, z_max, v_one],
                            [v_one, v_zero, v_one],
                            [v_one - z_max, v_one - z_max, v_one],
                            [v_one, v_one, v_one],
                            [v_one - center, center, v_one - center],
                            [v_one, x_max, v_one - x_max],
                            [v_one - center, v_one - center, v_one - center],
                            [v_one, v_one - x_max, v_one - x_max],
                        ]
                    ),
                )
            )

            # x_min_z_min
            spline_list.append(
                Bezier(
                    degrees=[1, 1, 1],
                    control_points=np.array(
                        [
                            [v_zero, x_min, x_min],
                            [center, center, center],
                            [v_zero, v_one - x_min, x_min],
                            [center, v_one - center, center],
                            [v_zero, v_zero, v_zero],
                            [z_min, z_min, v_zero],
                            [v_zero, v_one, v_zero],
                            [z_min, v_one - z_min, v_zero],
                        ]
                    ),
                )
            )

            # x_min_z_max
            spline_list.append(
                Bezier(
                    degrees=[1, 1, 1],
                    control_points=np.array(
                        [
                            [v_zero, v_zero, v_one],
                            [z_max, z_max, v_one],
                            [v_zero, v_one, v_one],
                            [z_max, v_one - z_max, v_one],
                            [v_zero, x_min, v_one - x_min],
                            [center, center, v_one - center],
                            [v_zero, v_one - x_min, v_one - x_min],
                            [center, v_one - center, v_one - center],
                        ]
                    ),
                )
            )

            # x_max_z_min
            spline_list.append(
                Bezier(
                    degrees=[1, 1, 1],
                    control_points=np.array(
                        [
                            [v_one - center, center, center],
                            [v_one, x_max, x_max],
                            [v_one - center, v_one - center, center],
                            [v_one, v_one - x_max, x_max],
                            [v_one - z_min, z_min, v_zero],
                            [v_one, v_zero, v_zero],
                            [v_one - z_min, v_one - z_min, v_zero],
                            [v_one, v_one, v_zero],
                        ]
                    ),
                )
            )

            # x_max_y_max
            spline_list.append(
                Bezier(
                    degrees=[1, 1, 1],
                    control_points=np.array(
                        [
                            [v_one - y_max, v_one, y_max],
                            [v_one, v_one, v_zero],
                            [v_one - y_max, v_one, v_one - y_max],
                            [v_one, v_one, v_one],
                            [v_one - center, v_one - center, center],
                            [v_one, v_one - x_max, x_max],
                            [v_one - center, v_one - center, v_one - center],
                            [v_one, v_one - x_max, v_one - x_max],
                        ]
                    ),
                )
            )

            # y_max_z_min
            spline_list.append(
                Bezier(
                    degrees=[1, 1, 1],
                    control_points=np.array(
                        [
                            [v_zero, v_one, v_zero],
                            [v_one, v_one, v_zero],
                            [y_max, v_one, y_max],
                            [v_one - y_max, v_one, y_max],
                            [z_min, v_one - z_min, v_zero],
                            [v_one - z_min, v_one - z_min, v_zero],
                            [center, v_one - center, center],
                            [v_one - center, v_one - center, center],
                        ]
                    ),
                )
            )

            # x_min_y_max
            spline_list.append(
                Bezier(
                    degrees=[1, 1, 1],
                    control_points=np.array(
                        [
                            [v_zero, v_one, v_zero],
                            [y_max, v_one, y_max],
                            [v_zero, v_one, v_one],
                            [y_max, v_one, v_one - y_max],
                            [v_zero, v_one - x_min, x_min],
                            [center, v_one - center, center],
                            [v_zero, v_one - x_min, v_one - x_min],
                            [center, v_one - center, v_one - center],
                        ]
                    ),
                )
            )

            # y_max_z_max
            spline_list.append(
                Bezier(
                    degrees=[1, 1, 1],
                    control_points=np.array(
                        [
                            [y_max, v_one, v_one - y_max],
                            [v_one - y_max, v_one, v_one - y_max],
                            [v_zero, v_one, v_one],
                            [v_one, v_one, v_one],
                            [center, v_one - center, v_one - center],
                            [v_one - center, v_one - center, v_one - center],
                            [z_max, v_one - z_max, v_one],
                            [v_one - z_max, v_one - z_max, v_one],
                        ]
                    ),
                )
            )

            # x_max_y_min
            spline_list.append(
                Bezier(
                    degrees=[1, 1, 1],
                    control_points=np.array(
                        [
                            [v_one - center, center, center],
                            [v_one, x_max, x_max],
                            [v_one - center, center, v_one - center],
                            [v_one, x_max, v_one - x_max],
                            [v_one - y_min, v_zero, y_min],
                            [v_one, v_zero, v_zero],
                            [v_one - y_min, v_zero, v_one - y_min],
                            [v_one, v_zero, v_one],
                        ]
                    ),
                )
            )

            # x_min_y_min
            spline_list.append(
                Bezier(
                    degrees=[1, 1, 1],
                    control_points=np.array(
                        [
                            [v_zero, x_min, x_min],
                            [center, center, center],
                            [v_zero, x_min, v_one - x_min],
                            [center, center, v_one - center],
                            [v_zero, v_zero, v_zero],
                            [y_min, v_zero, y_min],
                            [v_zero, v_zero, v_one],
                            [y_min, v_zero, v_one - y_min],
                        ]
                    ),
                )
            )

            # y_min_z_max
            spline_list.append(
                Bezier(
                    degrees=[1, 1, 1],
                    control_points=np.array(
                        [
                            [center, center, v_one - center],
                            [v_one - center, center, v_one - center],
                            [z_max, z_max, v_one],
                            [v_one - z_max, z_max, v_one],
                            [y_min, v_zero, v_one - y_min],
                            [v_one - y_min, v_zero, v_one - y_min],
                            [v_zero, v_zero, v_one],
                            [v_one, v_zero, v_one],
                        ]
                    ),
                )
            )

            # y_min_z_min
            spline_list.append(
                Bezier(
                    degrees=[1, 1, 1],
                    control_points=np.array(
                        [
                            [z_min, z_min, v_zero],
                            [v_one - z_min, z_min, v_zero],
                            [center, center, center],
                            [v_one - center, center, center],
                            [v_zero, v_zero, v_zero],
                            [v_one, v_zero, v_zero],
                            [y_min, v_zero, y_min],
                            [v_one - y_min, v_zero, y_min],
                        ]
                    ),
                )
            )

            # Pass to output
            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        # Return results
        if i_derivative == 0:
            return splines
        else:
            return (splines, derivatives)
