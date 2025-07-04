import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class HCross(_TileBase):
    """
    Looks like an H with a square in the middle extruded into the depth direction.
    """  # noqa: E501

    _dim = 3
    _para_dim = 3
    _evaluation_points = _np.array([[0.5, 0.5, 0.5]])
    _n_info_per_eval_point = 3
    _sensitivities_implemented = False
    _parameter_bounds = [[0.0, 0.3], [0.0, 0.15], [0.0, 0.3]]
    _parameters_shape = (1, 3)
    _default_parameter_value = _np.array([[0.25, 0.04, 0.15]])

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        **kwargs,  # noqa ARG002
    ):
        """Create a microtile based on the following parameters
            1) square width (sqw): side length of the square in the middle
            2) center roundness (cr): "roundness" of square in the middle
            3) beam width (bw): width of the side beams

        Parameters
        ----------
        parameters : np.array
            following parameters
            1) square width (sqw): side length of the square in the middle
            2) center roundness (cr): "roundness" of square in the middle
            3) beam width (bw): width of the side beams
        parameter_sensitivities: np.ndarray
          Sensitivities of the given parameters

        Returns
        -------
        microtile_list : list(splines)
        """
        parameters, n_derivatives, derivatives = self._process_input(
            parameters=parameters,
            parameter_sensitivities=parameter_sensitivities,
        )

        splines = []
        for i_derivative in range(n_derivatives + 1):
            # Parameter values for the spline
            if i_derivative == 0:
                sqw, cr, bw = parameters.flatten()

                v_zero = 0.0
                v_one_half = 0.5
                v_one = 1.0
            # Parameter values for the derivative
            else:
                v_zero = 0.0
                v_one_half = 0.0
                v_one = 0.0
                raise NotImplementedError(
                    "The derivatives have not been implemented"
                )

            # Init return value
            spline_list = []

            curved_tv3 = _Bezier(
                degrees=[2, 2],
                control_points=[
                    [v_one_half - sqw / 2, v_one_half + sqw / 2, v_zero],
                    [v_one_half, v_one_half + sqw / 2 - cr, v_zero],
                    [v_one_half + sqw / 2, v_one_half + sqw / 2, v_zero],
                    [
                        v_one_half - sqw / 2,
                        v_one_half + sqw / 2,
                        v_one_half - bw,
                    ],
                    [v_one_half, v_one_half + sqw / 2 - cr, v_one_half - bw],
                    [
                        v_one_half + sqw / 2,
                        v_one_half + sqw / 2,
                        v_one_half - bw,
                    ],
                    [v_zero, v_one - bw, v_one_half - bw],
                    [v_one_half, v_one - bw, v_one_half - bw],
                    [v_one, v_one - bw, v_one_half - bw],
                ],
            )

            curved_tv4 = _Bezier(
                degrees=[1, 1],
                control_points=[
                    [v_zero, v_one - bw, v_one_half - bw],
                    [v_one, v_one - bw, v_one_half - bw],
                    [v_zero, v_one, v_one_half - bw],
                    [v_one, v_one, v_one_half - bw],
                ],
            )

            curve_tv_middle = _Bezier(
                degrees=[2, 2],
                control_points=[
                    [v_one_half - sqw / 2, v_one_half - sqw / 2, v_zero],
                    [v_one_half, v_one_half - sqw / 2 + cr, v_zero],
                    [v_one_half + sqw / 2, v_one_half - sqw / 2, v_zero],
                    [v_one_half - sqw / 2 + cr, v_one_half, v_zero],
                    [v_one_half, v_one_half, v_zero],
                    [v_one_half + sqw / 2 - cr, v_one_half, v_zero],
                    [v_one_half - sqw / 2, v_one_half + sqw / 2, v_zero],
                    [v_one_half, v_one_half + sqw / 2 - cr, v_zero],
                    [v_one_half + sqw / 2, v_one_half + sqw / 2, v_zero],
                ],
            )

            curved_tv2 = _Bezier(
                degrees=[2, 2],
                control_points=[
                    [v_zero, bw, v_one_half - bw],
                    [v_one_half, bw, v_one_half - bw],
                    [v_one, bw, v_one_half - bw],
                    [
                        v_one_half - sqw / 2,
                        v_one_half - sqw / 2,
                        v_one_half - bw,
                    ],
                    [v_one_half, v_one_half - sqw / 2 + cr, v_one_half - bw],
                    [
                        v_one_half + sqw / 2,
                        v_one_half - sqw / 2,
                        v_one_half - bw,
                    ],
                    [v_one_half - sqw / 2, v_one_half - sqw / 2, v_zero],
                    [v_one_half, v_one_half - sqw / 2 + cr, v_zero],
                    [v_one_half + sqw / 2, v_one_half - sqw / 2, v_zero],
                ],
            )

            curved_tv1 = _Bezier(
                degrees=[1, 1],
                control_points=[
                    [v_zero, v_zero, v_one_half - bw],
                    [v_one, v_zero, v_one_half - bw],
                    [v_zero, bw, v_one_half - bw],
                    [v_one, bw, v_one_half - bw],
                ],
            )

            for tile_surface in [
                curve_tv_middle,
                curved_tv1,
                curved_tv2,
                curved_tv3,
                curved_tv4,
            ]:
                spline_cps = tile_surface.cps
                middle_cps = spline_cps.copy()
                middle_cps[:, -1] = v_one_half
                new_cps = _np.vstack((spline_cps, middle_cps))
                new_volume = _Bezier(
                    degrees=_np.append(tile_surface.ds, [1]),
                    control_points=new_cps,
                )
                spline_list.append(new_volume)
                # Add mirrored volume
                back_spline_cps = spline_cps.copy()
                back_spline_cps[:, -1] = 1 - back_spline_cps[:, -1]
                new_back_cps = _np.vstack((middle_cps, back_spline_cps))
                new_back_volume = new_volume.copy()
                new_back_volume.cps = new_back_cps
                spline_list.append(new_back_volume)

            # Pass to output
            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        # Return results
        return (splines, derivatives)
