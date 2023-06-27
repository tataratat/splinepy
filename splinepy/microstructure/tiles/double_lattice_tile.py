import numpy as np

from splinepy.bezier import Bezier
from splinepy.microstructure.tiles.tilebase import TileBase


class DoubleLatticeTile(TileBase):
    def __init__(self):
        """Simple crosstile with linear-quadratic branches and a trilinear
        center spline."""
        self._dim = 2
        self._evaluation_points = np.array(
            [
                [0.5, 0.5],
            ]
        )
        self._n_info_per_eval_point = 1

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        contact_length=0.5,
        **kwargs,
    ):
        """Create a microtile based on the parameters that describe the branch
        thicknesses.

        Thickness parameters are used to describe the inner radius of the
        outward facing branches

        Parameters
        ----------
        parameters : np.array
          only first entry is used, defines the internal radii of the
          branches
        parameter_sensitivities: np.ndarray
          correlates with thickness of branches and entouring wall
        contact_length : double
          required for conformity between tiles, sets the length of the center
          block on the tiles boundary

        Returns
        -------
        microtile_list : list(splines)
        """

        if not isinstance(contact_length, float):
            raise ValueError("Invalid Type")
        if not ((contact_length > 0.0) and (contact_length < 1.0)):
            raise ValueError("Contact length must be in (0.,1.)")

        # set to default if nothing is given
        if parameters is None:
            self._logd("Tile request is not parametrized, setting default 0.2")
            parameters = np.ones((1, 1)) * 0.1
        else:
            if not (np.all(parameters > 0) and np.all(parameters < 0.25)):
                raise ValueError("The parameter must be in 0.01 and 0.25")
            pass
        self.check_params(parameters)

        # Check if user requests derivative splines
        if self.check_param_derivatives(parameter_sensitivities):
            n_derivatives = parameter_sensitivities.shape[2]
        else:
            n_derivatives = 0

        derivatives = []
        splines = []
        for i_derivative in range(n_derivatives + 1):
            # Constant auxiliary values
            if i_derivative == 0:
                cl = contact_length
                pp = parameters[0, 0]  # parameters.shape == [1]
                v_one_half = 0.5
                v_one = 1.0
                v_zero = 0.0
            else:
                cl = 0.0
                pp = parameter_sensitivities[0, 0, i_derivative - 1]
                v_one_half = 0.0
                v_one = 0.0
                v_zero = 0.0

            # Set variables
            a01 = v_zero
            a02 = pp
            a03 = 2 * pp
            a04 = (v_one - cl) * 0.5
            a05 = v_one_half - pp
            a06 = v_one_half
            a07 = v_one_half + pp
            a08 = (v_one + cl) * 0.5
            a09 = v_one - 2 * pp
            a10 = v_one - pp
            a11 = v_one

            # Init return value
            spline_list = []

            # 1
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a01, a01],
                        [a02, a02],
                        [a01, a04],
                        [a02, a03],
                    ],
                )
            )

            # 2
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a01, a01],
                        [a04, a01],
                        [a02, a02],
                        [a03, a02],
                    ],
                )
            )

            # 3
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a04, a01],
                        [a08, a01],
                        [a03, a02],
                        [a09, a02],
                    ],
                )
            )

            # 4
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a08, a01],
                        [a11, a01],
                        [a09, a02],
                        [a10, a02],
                    ],
                ),
            )

            # 5
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a10, a02],
                        [a11, a01],
                        [a10, a03],
                        [a11, a04],
                    ],
                )
            )

            # 6
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a11, a04],
                        [a11, a08],
                        [a10, a03],
                        [a10, a09],
                    ],
                )
            )

            # 7
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a10, a09],
                        [a11, a08],
                        [a10, a10],
                        [a11, a11],
                    ],
                )
            )
            # 8
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a09, a10],
                        [a10, a10],
                        [a08, a11],
                        [a11, a11],
                    ],
                )
            )

            # 9
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a03, a10],
                        [a09, a10],
                        [a04, a11],
                        [a08, a11],
                    ],
                )
            )

            # 10
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a02, a10],
                        [a03, a10],
                        [a01, a11],
                        [a04, a11],
                    ],
                )
            )

            # 11
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a01, a08],
                        [a02, a09],
                        [a01, a11],
                        [a02, a10],
                    ],
                )
            )

            # 12
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a01, a04],
                        [a02, a03],
                        [a01, a08],
                        [a02, a09],
                    ],
                )
            )

            # 13
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a02, a09],
                        [a05, a06],
                        [a02, a10],
                        [a06, a06],
                    ],
                )
            )

            # 14
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a02, a10],
                        [a06, a06],
                        [a03, a10],
                        [a06, a07],
                    ],
                )
            )

            # 15
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a02, a02],
                        [a06, a06],
                        [a02, a03],
                        [a05, a06],
                    ],
                )
            )

            # 16
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a02, a02],
                        [a03, a02],
                        [a06, a06],
                        [a06, a05],
                    ],
                )
            )

            # 17
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a09, a02],
                        [a10, a02],
                        [a06, a05],
                        [a06, a06],
                    ],
                )
            )

            # 18
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a06, a06],
                        [a10, a02],
                        [a07, a06],
                        [a10, a03],
                    ],
                )
            )

            # 19
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a06, a06],
                        [a07, a06],
                        [a10, a10],
                        [a10, a09],
                    ],
                )
            )

            # 20
            spline_list.append(
                Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [a06, a06],
                        [a10, a10],
                        [a06, a07],
                        [a09, a10],
                    ],
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
