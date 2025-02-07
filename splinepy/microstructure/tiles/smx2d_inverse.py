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

                # Auxiliary values
                triangle_height = 0.5 * (v_one - bar_thickness - bar_thickness)
            else:
                raise NotImplementedError(
                    "Derivatives not implemented for this tile"
                )

            spline_list = []

            # Create triangle
            triangle_sw = _np.array(
                [
                    [v_zero, bar_thickness],
                    [v_zero, v_one_half],
                    [
                        triangle_height * 0.5,
                        v_one_half - triangle_height * 0.5,
                    ],
                    [triangle_height / 3, v_one_half],
                ]
            )

            triangle_se = _np.array(
                [
                    [v_zero, v_one_half],
                    [v_zero, v_one - bar_thickness],
                    [triangle_height / 3, v_one_half],
                    [
                        triangle_height * 0.5,
                        v_one_half + triangle_height * 0.5,
                    ],
                ]
            )
            triangle_n = _np.array(
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

            # Define bar control points
            sw_bar_bottom = _np.array(
                [
                    [v_zero, bar_thickness],
                    [
                        triangle_height * 0.5,
                        v_one_half - triangle_height * 0.5,
                    ],
                    [bar_thickness, v_zero],
                    [
                        v_one_half - triangle_height * 0.5,
                        triangle_height * 0.5,
                    ],
                ]
            )

            sw_bar_top = _np.array(
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

            # Transformation function to rotate recurrent patterns
            def rotate_function(cps):
                return _np.hstack(
                    (
                        cps[:, 1].reshape(-1, 1),
                        v_one - cps[:, 0].reshape(-1, 1),
                    )
                )

            for patch_yz_cps in [
                triangle_sw,
                triangle_se,
                triangle_n,
                sw_bar_bottom,
                sw_bar_top,
            ]:
                for _ in range(4):
                    patch_yz_cps = rotate_function(patch_yz_cps)
                    spline_list.append(
                        _Bezier(degrees=[1, 1], control_points=patch_yz_cps)
                    )

            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        return (splines, derivatives)
