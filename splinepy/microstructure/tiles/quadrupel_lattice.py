import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class QuadrupelLattice(_TileBase):
    def __init__(self):
        """
        Lattice base cell, consisting of a rectangle with two diagonals in the
        center with all directions having a different thickness associated to them
        """
        self._dim = 2
        self._para_dim = 2
        self._evaluation_points = _np.array(
            [
                [0.5, 0.5],
            ]
        )
        self._n_info_per_eval_point = 4

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        contact_length=0.5,
        **kwargs,  # noqa ARG002
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
        min_thickness = 0.1
        max_thickness = min_thickness * (1 + 2**0.5)
        if not isinstance(contact_length, float):
            raise ValueError("Invalid Type")
        if not ((contact_length > 0.0) and (contact_length < 1.0)):
            raise ValueError("Contact length must be in (0.,1.)")

        # set to default if nothing is given
        if parameters is None:
            self._logd("Tile request is not parametrized, setting default 0.2")
            parameters = _np.ones((1, 4)) * 0.1
        else:
            if not (
                _np.all(parameters > min_thickness)
                and _np.all(parameters < max_thickness)
            ):
                raise ValueError(
                    "The parameter must be in 0.01 and 1/(1 + sqrt(2))"
                )
            pass
        self.check_params(parameters)

        # Check if user requests derivative splines
        if self.check_param_derivatives(parameter_sensitivities):
            n_derivatives = parameter_sensitivities.shape[2]
            derivatives = []
        else:
            n_derivatives = 0
            derivatives = None

        # Coefficients
        coeff_sqrt_2 = 2**0.5
        coeff_inv_sqrt_2 = 2 ** (-0.5)
        coeff_one_half = 0.5

        splines = []
        for i_derivative in range(n_derivatives + 1):
            # Constant auxiliary values
            if i_derivative == 0:
                cl = contact_length
                thickness_horizontal = coeff_one_half * parameters[0, 0]
                thickness_vertical = coeff_one_half * parameters[0, 1]
                thickness_diag_up = coeff_one_half * parameters[0, 2]
                thickness_diag_down = coeff_one_half * parameters[0, 3]
                v_one_half = 0.5
                v_one = 1.0
                v_zero = 0.0
            else:
                cl = 0.0
                thickness_horizontal = (
                    coeff_one_half
                    * parameter_sensitivities[0, 0, i_derivative - 1]
                )
                thickness_vertical = (
                    coeff_one_half * parameters[0, 1, i_derivative - 1]
                )
                thickness_diag_up = (
                    coeff_one_half * parameters[0, 2, i_derivative - 1]
                )
                thickness_diag_down = (
                    coeff_one_half * parameters[0, 3, i_derivative - 1]
                )
                v_one_half = 0.0
                v_one = 0.0
                v_zero = 0.0

            # Set variables (with notation x01 = 1 - x11 )
            # Variables in Vertical direction
            v01 = thickness_horizontal
            v02 = thickness_vertical + coeff_sqrt_2 * thickness_diag_up
            v03 = thickness_vertical + coeff_sqrt_2 * thickness_diag_down
            v11 = v_one - v01
            v12 = v_one - v02
            v13 = v_one - v03

            # Variables in Horizontal direction
            h01 = thickness_vertical
            h02 = thickness_horizontal + coeff_sqrt_2 * thickness_diag_up
            h03 = thickness_horizontal + coeff_sqrt_2 * thickness_diag_down
            h11 = v_one - h01
            h12 = v_one - h02
            h13 = v_one - h03

            # Omnidirectional Variables
            c01 = v_zero
            c02 = (v_one - cl) * 0.5
            c03 = v_one_half + coeff_inv_sqrt_2 * (
                thickness_diag_up - thickness_diag_down
            )
            c04 = v_one_half - coeff_inv_sqrt_2 * (
                thickness_diag_up + thickness_diag_down
            )
            c05 = (
                0.5
                * coeff_one_half
                * (thickness_horizontal + thickness_vertical)
            )
            c06 = v_one_half
            c15 = v_one - c05
            c14 = v_one - c04
            c13 = v_one - c03
            c12 = v_one - c02
            c11 = v_one

            # Init return value
            spline_list = []

            # 0
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c01, c01],
                        [c05, c05],
                        [c01, c02],
                        [h01, v02],
                    ],
                )
            )

            # 1
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c01, c01],
                        [c02, c01],
                        [c05, c05],
                        [h02, v01],
                    ],
                )
            )

            # 2
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c02, c01],
                        [c12, c01],
                        [h02, v01],
                        [h13, v01],
                    ],
                )
            )

            # 3
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c12, c01],
                        [c11, c01],
                        [h13, v01],
                        [c15, c05],
                    ],
                ),
            )

            # 4
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c15, c05],
                        [c11, c01],
                        [h11, v03],
                        [c11, c02],
                    ],
                )
            )

            # 5
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [h11, v03],
                        [c11, c02],
                        [h11, v12],
                        [c11, c12],
                    ],
                )
            )

            # 6
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [h11, v12],
                        [c11, c12],
                        [c15, c15],
                        [c11, c11],
                    ],
                )
            )
            # 7
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [h12, v11],
                        [c15, c15],
                        [c12, c11],
                        [c11, c11],
                    ],
                )
            )

            # 8
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [h03, v11],
                        [h12, v11],
                        [c02, c11],
                        [c12, c11],
                    ],
                )
            )

            # 9
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c05, c15],
                        [h03, v11],
                        [c01, c11],
                        [c02, c11],
                    ],
                )
            )

            # 10
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c01, c12],
                        [h01, v13],
                        [c01, c11],
                        [c05, c15],
                    ],
                )
            )

            # 11
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c01, c02],
                        [h01, v02],
                        [c01, c12],
                        [h01, v13],
                    ],
                )
            )

            # 12
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [h01, v02],
                        [c05, c05],
                        [c04, c03],
                        [c06, c06],
                    ],
                )
            )

            # 13
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c05, c05],
                        [h02, v01],
                        [c06, c06],
                        [c03, c04],
                    ],
                )
            )

            # 14
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c03, c04],
                        [h13, v01],
                        [c06, c06],
                        [c15, c05],
                    ],
                )
            )

            # 15
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c06, c06],
                        [c15, c05],
                        [c14, c13],
                        [h11, v03],
                    ],
                )
            )

            # 16
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c06, c06],
                        [c14, c13],
                        [c15, c15],
                        [h11, v12],
                    ],
                )
            )

            # 17
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c06, c06],
                        [c15, c15],
                        [c13, c14],
                        [h12, v11],
                    ],
                )
            )

            # 19
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c05, c15],
                        [c06, c06],
                        [h03, v11],
                        [c13, c14],
                    ],
                )
            )

            # 20
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [h01, v13],
                        [c04, c03],
                        [c05, c15],
                        [c06, c06],
                    ],
                )
            )

            # Pass to output
            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        # Return results
        return (splines, derivatives)
