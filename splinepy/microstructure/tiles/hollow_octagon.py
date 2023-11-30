import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class HollowOctagon(_TileBase):
    def __init__(self):
        """Simple tile - looks like a nut"""
        self._dim = 2
        self._para_dim = 2
        self._evaluation_points = _np.array(
            [
                [0.5, 0.5],
            ]
        )
        self._n_info_per_eval_point = 1

    def _closing_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,  # TODO
        contact_length=0.2,
        closure=None,
    ):
        """Create a closing tile to match with closed surface.

        Parameters
        ----------
        parameters: np.ndarray(1, 1)
          One evaluation point with one parameter is used. This parameter
          specifies the distance from the center to the inner edge, where
          the value must be between 0.01 and 0.49.
        parameter_sensitivities: np.ndarray
          Describes the parameter sensitivities with respect to some design
          variable. In case the design variables directly apply to the
          parameter itself, they evaluate as delta_ij
        contact_length: float
          the length of the wall that contacts the other microstructure
        closure : str
          parametric dimension that needs to be closed, given in the form
          "x_min", "x_max", etc.

        Results
        -------
        spline_list : list
        """
        if closure is None:
            raise ValueError("No closing direction given")

        if parameters is None:
            self._log("Tile request is not parametrized, setting default 0.2")
            parameters = _np.array(
                _np.ones(
                    (len(self._evaluation_points), self._n_info_per_eval_point)
                )
                * 0.2
            )

        if not (_np.all(parameters > 0) and _np.all(parameters < 0.5)):
            raise ValueError(
                "The thickness of the wall must be in (0.01 and 0.49)"
            )

        self.check_params(parameters)

        if parameter_sensitivities is not None:
            raise NotImplementedError(
                "Derivatives are not implemented for this tile yet"
            )

        v_h_void = parameters[0, 0]
        if not ((v_h_void > 0.01) and (v_h_void < 0.5)):
            raise ValueError(
                "The thickness of the wall must be in (0.01 and 0.49)"
            )

        spline_list = []
        v_zero = 0.0
        v_one_half = 0.5
        v_one = 1.0
        v_outer_c_h = contact_length * 0.5
        v_inner_c_h = contact_length * parameters[0, 0]

        if closure == "x_min":
            # set points:
            right = _np.array(
                [
                    [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_one, -v_outer_c_h + v_one_half],
                    [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [v_one, v_outer_c_h + v_one_half],
                ]
            )

            right_top = _np.array(
                [
                    [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [v_one, v_outer_c_h + v_one_half],
                    [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [v_outer_c_h + v_one_half, v_one],
                ]
            )

            top = _np.array(
                [
                    [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [v_outer_c_h + v_one_half, v_one],
                    [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [-v_outer_c_h + v_one_half, v_one],
                ]
            )

            bottom_left = _np.array(
                [
                    [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_zero, v_zero],
                    [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [-v_outer_c_h + v_one_half, v_zero],
                ]
            )

            left = _np.array(
                [
                    [v_zero, v_zero],
                    [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_zero, v_one],
                    [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                ]
            )

            top_left = _np.array(
                [
                    [v_zero, v_one],
                    [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [-v_outer_c_h + v_one_half, v_one],
                    [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                ]
            )

            bottom = _np.array(
                [
                    [v_outer_c_h + v_one_half, v_zero],
                    [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [-v_outer_c_h + v_one_half, v_zero],
                    [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                ]
            )

            bottom_right = _np.array(
                [
                    [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [v_outer_c_h + v_one_half, v_zero],
                    [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_one, -v_outer_c_h + v_one_half],
                ]
            )

        elif closure == "x_max":
            right = _np.array(
                [
                    [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_one, v_zero],
                    [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [v_one, v_one],
                ]
            )

            right_top = _np.array(
                [
                    [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [v_one, v_one],
                    [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [v_outer_c_h + v_one_half, v_one],
                ]
            )

            top = _np.array(
                [
                    [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [v_outer_c_h + v_one_half, v_one],
                    [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [-v_outer_c_h + v_one_half, v_one],
                ]
            )

            bottom_left = _np.array(
                [
                    [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_zero, -v_outer_c_h + v_one_half],
                    [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [-v_outer_c_h + v_one_half, v_zero],
                ]
            )

            left = _np.array(
                [
                    [v_zero, -v_outer_c_h + v_one_half],
                    [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_zero, v_outer_c_h + v_one_half],
                    [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                ]
            )

            top_left = _np.array(
                [
                    [v_zero, v_outer_c_h + v_one_half],
                    [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [-v_outer_c_h + v_one_half, v_one],
                    [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                ]
            )

            bottom = _np.array(
                [
                    [v_outer_c_h + v_one_half, v_zero],
                    [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [-v_outer_c_h + v_one_half, v_zero],
                    [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                ]
            )

            bottom_right = _np.array(
                [
                    [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [v_outer_c_h + v_one_half, v_zero],
                    [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_one, v_zero],
                ]
            )

        elif closure == "y_min":
            # set points:
            right = _np.array(
                [
                    [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_one, -v_outer_c_h + v_one_half],
                    [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [v_one, v_outer_c_h + v_one_half],
                ]
            )

            right_top = _np.array(
                [
                    [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [v_one, v_outer_c_h + v_one_half],
                    [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [v_outer_c_h + v_one_half, v_one],
                ]
            )

            top = _np.array(
                [
                    [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [v_outer_c_h + v_one_half, v_one],
                    [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [-v_outer_c_h + v_one_half, v_one],
                ]
            )

            bottom_left = _np.array(
                [
                    [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_zero, -v_outer_c_h + v_one_half],
                    [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [v_zero, v_zero],
                ]
            )

            left = _np.array(
                [
                    [v_zero, -v_outer_c_h + v_one_half],
                    [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_zero, v_outer_c_h + v_one_half],
                    [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                ]
            )

            top_left = _np.array(
                [
                    [v_zero, v_outer_c_h + v_one_half],
                    [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [-v_outer_c_h + v_one_half, v_one],
                    [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                ]
            )

            bottom = _np.array(
                [
                    [v_one, v_zero],
                    [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [v_zero, v_zero],
                    [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                ]
            )

            bottom_right = _np.array(
                [
                    [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [v_one, v_zero],
                    [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_one, -v_outer_c_h + v_one_half],
                ]
            )

        elif closure == "y_max":
            # set points:
            right = _np.array(
                [
                    [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_one, -v_outer_c_h + v_one_half],
                    [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [v_one, v_outer_c_h + v_one_half],
                ]
            )

            right_top = _np.array(
                [
                    [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [v_one, v_outer_c_h + v_one_half],
                    [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [v_one, v_one],
                ]
            )

            top = _np.array(
                [
                    [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [v_one, v_one],
                    [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                    [v_zero, v_one],
                ]
            )

            bottom_left = _np.array(
                [
                    [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_zero, -v_outer_c_h + v_one_half],
                    [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [-v_outer_c_h + v_one_half, v_zero],
                ]
            )

            left = _np.array(
                [
                    [v_zero, -v_outer_c_h + v_one_half],
                    [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_zero, v_outer_c_h + v_one_half],
                    [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                ]
            )

            top_left = _np.array(
                [
                    [v_zero, v_outer_c_h + v_one_half],
                    [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                    [v_zero, v_one],
                    [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                ]
            )

            bottom = _np.array(
                [
                    [v_outer_c_h + v_one_half, v_zero],
                    [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [-v_outer_c_h + v_one_half, v_zero],
                    [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                ]
            )

            bottom_right = _np.array(
                [
                    [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                    [v_outer_c_h + v_one_half, v_zero],
                    [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                    [v_one, -v_outer_c_h + v_one_half],
                ]
            )

        spline_list.append(_Bezier(degrees=[1, 1], control_points=right))

        spline_list.append(_Bezier(degrees=[1, 1], control_points=right_top))

        spline_list.append(_Bezier(degrees=[1, 1], control_points=bottom))

        spline_list.append(_Bezier(degrees=[1, 1], control_points=bottom_left))

        spline_list.append(_Bezier(degrees=[1, 1], control_points=left))

        spline_list.append(_Bezier(degrees=[1, 1], control_points=top_left))

        spline_list.append(_Bezier(degrees=[1, 1], control_points=top))

        spline_list.append(
            _Bezier(degrees=[1, 1], control_points=bottom_right)
        )

        return (spline_list, None)

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,  # TODO
        contact_length=0.2,
        closure=None,
        **kwargs,  # noqa ARG002
    ):
        """Create a microtile based on the parameters that describe the wall
        thicknesses.

        Thickness parameters are used to describe the inner radius of the
        outward facing branches

        Parameters
        ----------
        parameters : np.array(1, 1)
          One evaluation point with one parameter is used. This parameter
          specifies the distance from the center to the inner edge, where
          the value must be between non-inclusive (0, 0.5).
        parameter_sensitivities: np.ndarray
          Describes the parameter sensitivities with respect to some design
          variable. In case the design variables directly apply to the
          parameter itself, they evaluate as delta_ij
        contact_length : float
            the length of the wall that contacts the other microstructure
        closure : str
          parametric dimension that needs to be closed, given in the form
          "x_min", "x_max", etc.

        Returns
        -------
        microtile_list : list(splines)
        derivatives: list<list<splines>> / None
        """

        if not isinstance(contact_length, float):
            raise ValueError("Invalid Type for radius")

        if not ((contact_length > 0) and (contact_length < 0.99)):
            raise ValueError("The length of a side must be in (0.01, 0.99)")

        if parameters is None:
            self._logd("Setting parameters to default values (0.2)")
            parameters = _np.array(
                _np.ones(
                    (len(self._evaluation_points), self._n_info_per_eval_point)
                )
                * 0.2
            )

        if parameter_sensitivities is not None:
            raise NotImplementedError(
                "Derivatives are not implemented for this tile yet"
            )

        self.check_params(parameters)

        v_h_void = parameters[0, 0]
        if not (_np.all(parameters > 0) and _np.all(parameters < 0.5)):
            raise ValueError(
                "The thickness of the wall must be in (0.01 and 0.49)"
            )

        if closure is not None:
            return self._closing_tile(
                parameters=parameters,
                parameter_sensitivities=parameter_sensitivities,  # TODO
                contact_length=contact_length,
                closure=closure,
            )

        v_zero = 0.0
        v_one_half = 0.5
        v_one = 1.0
        v_outer_c_h = contact_length * 0.5
        v_inner_c_h = contact_length * parameters[0, 0]

        spline_list = []

        # set points:
        right = _np.array(
            [
                [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                [v_one, -v_outer_c_h + v_one_half],
                [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                [v_one, v_outer_c_h + v_one_half],
            ]
        )

        right_top = _np.array(
            [
                [v_h_void + v_one_half, v_inner_c_h + v_one_half],
                [v_one, v_outer_c_h + v_one_half],
                [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                [v_outer_c_h + v_one_half, v_one],
            ]
        )

        top = _np.array(
            [
                [v_inner_c_h + v_one_half, v_h_void + v_one_half],
                [v_outer_c_h + v_one_half, v_one],
                [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
                [-v_outer_c_h + v_one_half, v_one],
            ]
        )

        bottom_left = _np.array(
            [
                [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                [v_zero, -v_outer_c_h + v_one_half],
                [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                [-v_outer_c_h + v_one_half, v_zero],
            ]
        )

        left = _np.array(
            [
                [v_zero, -v_outer_c_h + v_one_half],
                [-v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                [v_zero, v_outer_c_h + v_one_half],
                [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
            ]
        )

        top_left = _np.array(
            [
                [v_zero, v_outer_c_h + v_one_half],
                [-v_h_void + v_one_half, v_inner_c_h + v_one_half],
                [-v_outer_c_h + v_one_half, v_one],
                [-v_inner_c_h + v_one_half, v_h_void + v_one_half],
            ]
        )

        bottom = _np.array(
            [
                [v_outer_c_h + v_one_half, v_zero],
                [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                [-v_outer_c_h + v_one_half, v_zero],
                [-v_inner_c_h + v_one_half, -v_h_void + v_one_half],
            ]
        )

        bottom_right = _np.array(
            [
                [v_inner_c_h + v_one_half, -v_h_void + v_one_half],
                [v_outer_c_h + v_one_half, v_zero],
                [v_h_void + v_one_half, -v_inner_c_h + v_one_half],
                [v_one, -v_outer_c_h + v_one_half],
            ]
        )

        spline_list.append(_Bezier(degrees=[1, 1], control_points=right))

        spline_list.append(_Bezier(degrees=[1, 1], control_points=right_top))

        spline_list.append(_Bezier(degrees=[1, 1], control_points=bottom))

        spline_list.append(_Bezier(degrees=[1, 1], control_points=bottom_left))

        spline_list.append(_Bezier(degrees=[1, 1], control_points=left))

        spline_list.append(_Bezier(degrees=[1, 1], control_points=top_left))

        spline_list.append(_Bezier(degrees=[1, 1], control_points=top))

        spline_list.append(
            _Bezier(degrees=[1, 1], control_points=bottom_right)
        )

        return (spline_list, None)
