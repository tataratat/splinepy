import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class Snappy(_TileBase):
    def __init__(self):
        """Snap-through tile consisting of a thin truss and a thick truss that
        collide into each other"""
        self._dim = 2
        self._para_dim = 2
        # Dummy point
        self._evaluation_points = _np.array([[0.5, 0.5]])
        self._n_info_per_eval_point = 1

    def _closing_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,  # TODO
        closure=None,
        contact_length=0.1,
        a=0.1,
        b=0.2,
        r=0.15,
        **kwargs,  # noqa ARG002
    ):
        """Create a closing tile to match with closed surface.

        Parameters
        ----------
        parameters: np.ndarray
          currently unused
        parameter_sensitivities: np.ndarray
          Describes the parameter sensitivities with respect to some design
          variable. In case the design variables directly apply to the
          parameter itself, they evaluate as delta_ij unused
        closure : str
          string specifying the closing dimensions (e.g., x_min)
        contact_length : float
          the length of the wall that contacts neighboring microstructures
        a : float
          height/ thickness of the thinner/upper beam
        b : float
          height/ thickness of the lower/thicker beam
        r : float
          'radius' of the cubic bezier

        Results
        -------
        spline_list : list
        derivatives: list<list<splines>> / None
        """
        if closure is None:
            raise ValueError("No closing direction given")

        self.check_params(parameters)

        if parameter_sensitivities is not None:
            raise NotImplementedError(
                "Derivatives are not implemented for this tile yet"
            )

        spline_list = []
        v_zero = 0.0
        v_one_half = 0.5
        v_one = 1.0
        cl_2 = contact_length * 0.5
        cl_2_inv = 1 - contact_length * 0.5
        a_inv = v_one - a

        if closure == "y_min":
            # set points:
            spline_1 = _np.array(
                [
                    [v_zero, v_zero],
                    [cl_2, v_zero],
                    [v_zero, b],
                    [cl_2, b],
                ]
            )

            spline_list.append(
                _Bezier(degrees=[1, 1], control_points=spline_1)
            )
            spline_2 = _np.array(
                [
                    [cl_2_inv, v_zero],
                    [v_one, v_zero],
                    [cl_2_inv, b],
                    [v_one, b],
                ]
            )

            spline_list.append(
                _Bezier(degrees=[1, 1], control_points=spline_2)
            )
            spline_3 = _np.array(
                [
                    [v_zero, a_inv],
                    [cl_2, a_inv],
                    [v_zero, v_one],
                    [cl_2, v_one],
                ]
            )

            spline_list.append(
                _Bezier(degrees=[1, 1], control_points=spline_3)
            )
            spline_4 = _np.array(
                [
                    [cl_2_inv, a_inv],
                    [v_one, a_inv],
                    [cl_2_inv, v_one],
                    [v_one, v_one],
                ]
            )

            spline_list.append(
                _Bezier(degrees=[1, 1], control_points=spline_4)
            )

            spline_5 = _np.array(
                [
                    [v_one_half - cl_2, v_zero],
                    [v_one_half + cl_2, v_zero],
                    [v_one_half - cl_2, v_one_half],
                    [v_one_half + cl_2, v_one_half],
                ]
            )

            spline_list.append(
                _Bezier(degrees=[1, 1], control_points=spline_5)
            )

            spline_6 = _np.array(
                [
                    [v_one_half - cl_2, v_one_half],
                    [v_one_half + cl_2, v_one_half],
                    [v_one_half - cl_2, v_one_half + a],
                    [v_one_half + cl_2, v_one_half + a],
                ]
            )

            spline_list.append(
                _Bezier(degrees=[1, 1], control_points=spline_6)
            )

            spline_7 = _np.array(
                [
                    [cl_2, v_zero],
                    [cl_2 + r, v_zero],
                    [v_one_half - cl_2 - r, v_zero],
                    [v_one_half - cl_2, v_zero],
                    [cl_2, b],
                    [cl_2 + r, b],
                    [v_one_half - cl_2 - r, v_one_half],
                    [v_one_half - cl_2, v_one_half],
                ]
            )

            spline_list.append(
                _Bezier(degrees=[3, 1], control_points=spline_7)
            )

            spline_8 = _np.array(
                [
                    [cl_2, v_zero],
                    [cl_2 + r, v_zero],
                    [v_one_half - cl_2 - r, v_zero],
                    [v_one_half - cl_2, v_zero],
                    [cl_2, v_one_half],
                    [cl_2 + r, v_one_half],
                    [v_one_half - cl_2 - r, b],
                    [v_one_half - cl_2, b],
                ]
            ) + [v_one_half, v_zero]

            spline_list.append(
                _Bezier(degrees=[3, 1], control_points=spline_8)
            )

            spline_9 = _np.array(
                [
                    [cl_2, a_inv],
                    [cl_2 + r, a_inv],
                    [v_one_half - cl_2 - r, v_one_half],
                    [v_one_half - cl_2, v_one_half],
                    [cl_2, v_one],
                    [cl_2 + r, v_one],
                    [v_one_half - cl_2 - r, v_one_half + a],
                    [v_one_half - cl_2, v_one_half + a],
                ]
            )

            spline_list.append(
                _Bezier(degrees=[3, 1], control_points=spline_9)
            )

            spline_10 = _np.array(
                [
                    [cl_2, v_one_half],
                    [cl_2 + r, v_one_half],
                    [v_one_half - cl_2 - r, a_inv],
                    [v_one_half - cl_2, a_inv],
                    [cl_2, v_one_half + a],
                    [cl_2 + r, v_one_half + a],
                    [v_one_half - cl_2 - r, v_one],
                    [v_one_half - cl_2, v_one],
                ]
            ) + [v_one_half, v_zero]

            spline_list.append(
                _Bezier(degrees=[3, 1], control_points=spline_10)
            )
            return spline_list
        elif closure == "y_max":
            spline_1 = _np.array(
                [
                    [v_zero, v_zero],
                    [cl_2, v_zero],
                    [v_zero, v_one],
                    [cl_2, v_one],
                ]
            )
            spline_list.append(
                _Bezier(degrees=[1, 1], control_points=spline_1)
            )
            spline_2 = _np.array(
                [
                    [cl_2_inv, v_zero],
                    [v_one, v_zero],
                    [cl_2_inv, v_one],
                    [v_one, v_one],
                ]
            )
            spline_list.append(
                _Bezier(degrees=[1, 1], control_points=spline_2)
            )
            spline_3 = _np.array(
                [
                    [v_one_half - cl_2, v_one_half - b],
                    [v_one_half + cl_2, v_one_half - b],
                    [v_one_half - cl_2, v_one],
                    [v_one_half + cl_2, v_one],
                ]
            )
            spline_list.append(
                _Bezier(degrees=[1, 1], control_points=spline_3)
            )
            spline_4 = _np.array(
                [
                    [cl_2, v_zero],
                    [cl_2 + r, v_zero],
                    [v_one_half - cl_2 - r, v_one_half - b],
                    [v_one_half - cl_2, v_one_half - b],
                    [cl_2, v_one],
                    [cl_2 + r, v_one],
                    [v_one_half - cl_2 - r, v_one],
                    [v_one_half - cl_2, v_one],
                ]
            )
            spline_list.append(
                _Bezier(degrees=[3, 1], control_points=spline_4)
            )
            spline_5 = _np.array(
                [
                    [cl_2, v_one_half - b],
                    [cl_2 + r, v_one_half - b],
                    [v_one_half - cl_2 - r, v_zero],
                    [v_one_half - cl_2, v_zero],
                    [cl_2, v_one],
                    [cl_2 + r, v_one],
                    [v_one_half - cl_2 - r, v_one],
                    [v_one_half - cl_2, v_one],
                ]
            ) + [v_one_half, v_zero]

            spline_list.append(
                _Bezier(degrees=[3, 1], control_points=spline_5)
            )
            return (spline_list, None)
        else:
            raise ValueError(
                "Closing tile is only implemented for y-enclosure"
            )

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,  # TODO
        contact_length=0.1,
        a=0.1,
        b=0.2,
        c=0.3,
        r=0.15,
        closure=None,
        **kwargs,  # noqa ARG002
    ):
        """Create a microtile based on the parameters that describe the wall
        thicknesses.

        Thickness parameters are used to describe the inner radius of the
        outward facing branches

        Parameters
        ----------
        parameters : np.array
          Currently, no parameter is used, (First test)
        parameter_sensitivities: list(np.ndarray)
          Describes the parameter sensitivities with respect to some design
          variable. In case the design variables directly apply to the
          parameter itself, they evaluate as delta_ij, currently unused
        contact_length : float
          the length of the wall that contacts neighboring microstructures
        a : float
          height/ thickness of the thinner/upper beam
        b : float
          height/ thickness of the lower/thicker beam
        c : float
          offset to the upper beam (for consistent snap-through must fulfill
          2*c<1-b)
        r : float
          'radius' of the cubic bezier
        closure : str
          string specifying the closing dimensions (e.g., x_min)

        Returns
        -------
        microtile_list : list(splines)
        derivatives: list<list<splines>> / None
        """

        for param in [a, b, c, r, contact_length]:
            if not isinstance(param, float):
                raise ValueError(f"Invalid Type, {param} is not float")
            if param < 0:
                raise ValueError("Invalid parameter, must be > 0.")

        if not ((contact_length > 0) and (contact_length < 0.49)):
            raise ValueError("The length of a side must be in (0.01, 0.49)")

        # Check horizontal parameters
        if not ((r + contact_length) < 0.5):
            raise ValueError(
                "Inconsistent parameters, must fulfill : 2*r + contact_length"
                " < 0.5"
            )

        # Check vertical parameters
        if not ((2 * c + b) < 1.0) or a > c:
            raise ValueError(
                "Inconsistent parameters, must be 2*c<1-c and a<c"
            )

        if parameters is not None:
            raise NotImplementedError(
                "Parametriazation is not implemented for this tile yet"
            )

        if parameter_sensitivities is not None:
            raise NotImplementedError(
                "Derivatives are not implemented for this tile yet"
            )

        if closure is not None:
            return self._closing_tile(
                parameters=None,
                parameter_sensitivities=parameter_sensitivities,
                closure=closure,
                contact_length=contact_length,
                a=a,
                b=b,
                r=r,
                **kwargs,
            )

        v_zero = 0.0
        v_one_half = 0.5
        v_one = 1.0
        cl_2 = contact_length * 0.5
        cl_2_inv = 1 - contact_length * 0.5
        a_inv = v_one - a

        spline_list = []

        # set points:
        spline_1 = _np.array(
            [
                [v_zero, v_zero],
                [cl_2, v_zero],
                [v_zero, b],
                [cl_2, b],
            ]
        )

        spline_list.append(_Bezier(degrees=[1, 1], control_points=spline_1))
        spline_2 = _np.array(
            [
                [cl_2_inv, v_zero],
                [v_one, v_zero],
                [cl_2_inv, b],
                [v_one, b],
            ]
        )

        spline_list.append(_Bezier(degrees=[1, 1], control_points=spline_2))
        spline_3 = _np.array(
            [
                [v_zero, a_inv],
                [cl_2, a_inv],
                [v_zero, v_one],
                [cl_2, v_one],
            ]
        )

        spline_list.append(_Bezier(degrees=[1, 1], control_points=spline_3))
        spline_4 = _np.array(
            [
                [cl_2_inv, a_inv],
                [v_one, a_inv],
                [cl_2_inv, v_one],
                [v_one, v_one],
            ]
        )

        spline_list.append(_Bezier(degrees=[1, 1], control_points=spline_4))

        spline_5 = _np.array(
            [
                [v_one_half - cl_2, v_one_half - b],
                [v_one_half + cl_2, v_one_half - b],
                [v_one_half - cl_2, v_one_half],
                [v_one_half + cl_2, v_one_half],
            ]
        )

        spline_list.append(_Bezier(degrees=[1, 1], control_points=spline_5))

        spline_6 = _np.array(
            [
                [v_one_half - cl_2, v_one_half],
                [v_one_half + cl_2, v_one_half],
                [v_one_half - cl_2, v_one_half + a],
                [v_one_half + cl_2, v_one_half + a],
            ]
        )

        spline_list.append(_Bezier(degrees=[1, 1], control_points=spline_6))

        spline_7 = _np.array(
            [
                [cl_2, v_zero],
                [cl_2 + r, v_zero],
                [v_one_half - cl_2 - r, v_one_half - b],
                [v_one_half - cl_2, v_one_half - b],
                [cl_2, b],
                [cl_2 + r, b],
                [v_one_half - cl_2 - r, v_one_half],
                [v_one_half - cl_2, v_one_half],
            ]
        )

        spline_list.append(_Bezier(degrees=[3, 1], control_points=spline_7))

        spline_8 = _np.array(
            [
                [cl_2, v_one_half - b],
                [cl_2 + r, v_one_half - b],
                [v_one_half - cl_2 - r, v_zero],
                [v_one_half - cl_2, v_zero],
                [cl_2, v_one_half],
                [cl_2 + r, v_one_half],
                [v_one_half - cl_2 - r, b],
                [v_one_half - cl_2, b],
            ]
        ) + [v_one_half, v_zero]

        spline_list.append(_Bezier(degrees=[3, 1], control_points=spline_8))

        spline_9 = _np.array(
            [
                [cl_2, a_inv],
                [cl_2 + r, a_inv],
                [v_one_half - cl_2 - r, v_one_half],
                [v_one_half - cl_2, v_one_half],
                [cl_2, v_one],
                [cl_2 + r, v_one],
                [v_one_half - cl_2 - r, v_one_half + a],
                [v_one_half - cl_2, v_one_half + a],
            ]
        )

        spline_list.append(_Bezier(degrees=[3, 1], control_points=spline_9))

        spline_10 = _np.array(
            [
                [cl_2, v_one_half],
                [cl_2 + r, v_one_half],
                [v_one_half - cl_2 - r, a_inv],
                [v_one_half - cl_2, a_inv],
                [cl_2, v_one_half + a],
                [cl_2 + r, v_one_half + a],
                [v_one_half - cl_2 - r, v_one],
                [v_one_half - cl_2, v_one],
            ]
        ) + [v_one_half, v_zero]

        spline_list.append(_Bezier(degrees=[3, 1], control_points=spline_10))

        return (spline_list, None)
