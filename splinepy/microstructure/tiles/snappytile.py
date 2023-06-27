import numpy as np

from splinepy.bezier import Bezier
from splinepy.microstructure.tiles.tilebase import TileBase


class SnappyTile(TileBase):
    def __init__(self):
        """Snap-through tile consisting of a thin truss and a thick truss that
        collide into each other"""
        self._dim = 2
        # Dummy point
        self._evaluation_points = np.array([[0.5, 0.5]])
        self._n_info_per_eval_point = 1

    def closing_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        closure=None,
        contact_length=0.1,
        a=0.1,
        b=0.2,
        c=0.3,
        r=0.15,
        **kwargs,
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
        c : float
          offset to the upper beam (for consistent snap-through must fulfill
          2*c<1-b)
        r : float
          'radius' of the cubic bezier

        Results
        -------
        spline_list : list
        """
        if closure is None:
            raise ValueError("No closing direction given")

        self.check_params(parameters)

        spline_list = []
        v_zero = 0.0
        v_one_half = 0.5
        v_one = 1.0
        cl_2 = contact_length * 0.5
        cl_2_inv = 1 - contact_length * 0.5
        a_inv = v_one - a

        if closure == "y_min":
            # set points:
            spline_1 = np.array(
                [
                    [v_zero, v_zero],
                    [cl_2, v_zero],
                    [v_zero, b],
                    [cl_2, b],
                ]
            )

            spline_list.append(Bezier(degrees=[1, 1], control_points=spline_1))
            spline_2 = np.array(
                [
                    [cl_2_inv, v_zero],
                    [v_one, v_zero],
                    [cl_2_inv, b],
                    [v_one, b],
                ]
            )

            spline_list.append(Bezier(degrees=[1, 1], control_points=spline_2))
            spline_3 = np.array(
                [
                    [v_zero, a_inv],
                    [cl_2, a_inv],
                    [v_zero, v_one],
                    [cl_2, v_one],
                ]
            )

            spline_list.append(Bezier(degrees=[1, 1], control_points=spline_3))
            spline_4 = np.array(
                [
                    [cl_2_inv, a_inv],
                    [v_one, a_inv],
                    [cl_2_inv, v_one],
                    [v_one, v_one],
                ]
            )

            spline_list.append(Bezier(degrees=[1, 1], control_points=spline_4))

            spline_5 = np.array(
                [
                    [v_one_half - cl_2, v_zero],
                    [v_one_half + cl_2, v_zero],
                    [v_one_half - cl_2, v_one_half],
                    [v_one_half + cl_2, v_one_half],
                ]
            )

            spline_list.append(Bezier(degrees=[1, 1], control_points=spline_5))

            spline_6 = np.array(
                [
                    [v_one_half - cl_2, v_one_half],
                    [v_one_half + cl_2, v_one_half],
                    [v_one_half - cl_2, v_one_half + a],
                    [v_one_half + cl_2, v_one_half + a],
                ]
            )

            spline_list.append(Bezier(degrees=[1, 1], control_points=spline_6))

            spline_7 = np.array(
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

            spline_list.append(Bezier(degrees=[3, 1], control_points=spline_7))

            spline_8 = np.array(
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

            spline_list.append(Bezier(degrees=[3, 1], control_points=spline_8))

            spline_9 = np.array(
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

            spline_list.append(Bezier(degrees=[3, 1], control_points=spline_9))

            spline_10 = np.array(
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
                Bezier(degrees=[3, 1], control_points=spline_10)
            )
            return spline_list
        elif closure == "y_max":
            spline_1 = np.array(
                [
                    [v_zero, v_zero],
                    [cl_2, v_zero],
                    [v_zero, v_one],
                    [cl_2, v_one],
                ]
            )
            spline_list.append(Bezier(degrees=[1, 1], control_points=spline_1))
            spline_2 = np.array(
                [
                    [cl_2_inv, v_zero],
                    [v_one, v_zero],
                    [cl_2_inv, v_one],
                    [v_one, v_one],
                ]
            )
            spline_list.append(Bezier(degrees=[1, 1], control_points=spline_2))
            spline_3 = np.array(
                [
                    [v_one_half - cl_2, v_one_half - b],
                    [v_one_half + cl_2, v_one_half - b],
                    [v_one_half - cl_2, v_one],
                    [v_one_half + cl_2, v_one],
                ]
            )
            spline_list.append(Bezier(degrees=[1, 1], control_points=spline_3))
            spline_4 = np.array(
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
            spline_list.append(Bezier(degrees=[3, 1], control_points=spline_4))
            spline_5 = np.array(
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

            spline_list.append(Bezier(degrees=[3, 1], control_points=spline_5))
            return spline_list
        else:
            raise ValueError(
                "Closing tile is only implemented for y-enclosure"
            )

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        contact_length=0.1,
        a=0.1,
        b=0.2,
        c=0.3,
        r=0.15,
        **kwargs,
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

        Returns
        -------
        microtile_list : list(splines)
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
                "Inconsistent parameters, must fullfil : 2*r + contact_length"
                " < 0.5"
            )

        # Check vertical parameters
        if not ((2 * c + b) < 1.0) or a > c:
            raise ValueError(
                "Inconsistent parameters, must be 2*c<1-c and a<c"
            )

        if parameters is None:
            self._logd("Setting parameters to default values (0.2)")
            parameters = np.array(
                np.ones(
                    (len(self._evaluation_points), self._n_info_per_eval_point)
                )
                * 0.2
            )

        self.check_params(parameters)

        v_zero = 0.0
        v_one_half = 0.5
        v_one = 1.0
        cl_2 = contact_length * 0.5
        cl_2_inv = 1 - contact_length * 0.5
        a_inv = v_one - a

        spline_list = []

        # set points:
        spline_1 = np.array(
            [
                [v_zero, v_zero],
                [cl_2, v_zero],
                [v_zero, b],
                [cl_2, b],
            ]
        )

        spline_list.append(Bezier(degrees=[1, 1], control_points=spline_1))
        spline_2 = np.array(
            [
                [cl_2_inv, v_zero],
                [v_one, v_zero],
                [cl_2_inv, b],
                [v_one, b],
            ]
        )

        spline_list.append(Bezier(degrees=[1, 1], control_points=spline_2))
        spline_3 = np.array(
            [
                [v_zero, a_inv],
                [cl_2, a_inv],
                [v_zero, v_one],
                [cl_2, v_one],
            ]
        )

        spline_list.append(Bezier(degrees=[1, 1], control_points=spline_3))
        spline_4 = np.array(
            [
                [cl_2_inv, a_inv],
                [v_one, a_inv],
                [cl_2_inv, v_one],
                [v_one, v_one],
            ]
        )

        spline_list.append(Bezier(degrees=[1, 1], control_points=spline_4))

        spline_5 = np.array(
            [
                [v_one_half - cl_2, v_one_half - b],
                [v_one_half + cl_2, v_one_half - b],
                [v_one_half - cl_2, v_one_half],
                [v_one_half + cl_2, v_one_half],
            ]
        )

        spline_list.append(Bezier(degrees=[1, 1], control_points=spline_5))

        spline_6 = np.array(
            [
                [v_one_half - cl_2, v_one_half],
                [v_one_half + cl_2, v_one_half],
                [v_one_half - cl_2, v_one_half + a],
                [v_one_half + cl_2, v_one_half + a],
            ]
        )

        spline_list.append(Bezier(degrees=[1, 1], control_points=spline_6))

        spline_7 = np.array(
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

        spline_list.append(Bezier(degrees=[3, 1], control_points=spline_7))

        spline_8 = np.array(
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

        spline_list.append(Bezier(degrees=[3, 1], control_points=spline_8))

        spline_9 = np.array(
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

        spline_list.append(Bezier(degrees=[3, 1], control_points=spline_9))

        spline_10 = np.array(
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

        spline_list.append(Bezier(degrees=[3, 1], control_points=spline_10))

        return spline_list
