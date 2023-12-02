import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class Armadillo(_TileBase):
    def __init__(self):
        """
        Tile in the in a shape of a multisided dice, where every side connects
        to the neighbor tile in the center of the surface
        """
        self._dim = 3
        self._para_dim = 3
        self._evaluation_points = _np.array(
            [
                [0.5, 0.5, 0.5],
            ]
        )
        self._n_info_per_eval_point = 1

    def closing_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,  # TODO
        contact_length=0.3,
        closure=None,
        **kwargs,  # noqa ARG002
    ):
        """Create a closing tile to match with closed surface.

        Parameters
        ----------
        parameters: np.ndarray
            An evaluation point with one parameter is used. This parameter
            describes the thickness of the wall. The parameters must be a
            two-dimensional np.array, where the value must be between 0.01
            and 0.49
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
        derivative_list : list / None
        """
        if closure is None:
            raise ValueError("No closing direction given")

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

        self.check_params(parameters)

        if parameter_sensitivities is not None:
            raise NotImplementedError(
                "Derivatives are not implemented for this tile yet"
            )

        if not (_np.all(parameters > 0) and _np.all(parameters < 0.5)):
            raise ValueError(
                "The thickness of the wall must be in (0.01 and 0.49)"
            )

        v_wall_thickness = parameters[0, 0]
        spline_list = []
        v_zero = 0.0
        v_one_half = 0.5
        v_one = 1.0
        v_half_contact_length = contact_length * 0.5
        v_inner_half_contact_length = contact_length * parameters[0, 0]

        if closure == "x_min":
            # set points:
            right = _np.array(
                [
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            connection_front_right = _np.array(
                [
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            front = _np.array(
                [
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            connection_back_left = _np.array(
                [
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_one,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_zero,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            left = _np.array(
                [
                    [
                        v_zero,
                        v_one,
                        v_one,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_one,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one,
                        v_zero,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_zero,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_front_left = _np.array(
                [
                    [
                        v_zero,
                        v_one,
                        v_one,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one,
                        v_zero,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            back = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_back_right = _np.array(
                [
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            bottom = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            top = _np.array(
                [
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                ]
            )

            connection_front_bottom = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            connection_front_top = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                ]
            )

            connection_back_bottom = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_back_top = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_top_right = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one,
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_top_left = _np.array(
                [
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_one,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_one,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_bottom_left = _np.array(
                [
                    [
                        v_zero,
                        v_one,
                        v_zero,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            connection_bottom_right = _np.array(
                [
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

        elif closure == "x_max":
            # set points:
            right = _np.array(
                [
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_one,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_zero,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one,
                        v_one,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one,
                        v_zero,
                    ],
                ]
            )

            connection_front_right = _np.array(
                [
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one,
                        v_one,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one,
                        v_zero,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            front = _np.array(
                [
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            connection_back_left = _np.array(
                [
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            left = _np.array(
                [
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_front_left = _np.array(
                [
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            back = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_back_right = _np.array(
                [
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_one,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_zero,
                    ],
                ]
            )

            bottom = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            top = _np.array(
                [
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                ]
            )

            connection_front_bottom = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            connection_front_top = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                ]
            )

            connection_back_bottom = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_back_top = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_top_right = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one,
                        v_one,
                        v_one,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_one,
                    ],
                    [
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_top_left = _np.array(
                [
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_bottom_left = _np.array(
                [
                    [
                        v_zero,
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            connection_bottom_right = _np.array(
                [
                    [
                        v_one,
                        v_one,
                        v_zero,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_zero,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

        elif closure == "y_min":
            # set points:
            right = _np.array(
                [
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            connection_front_right = _np.array(
                [
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            front = _np.array(
                [
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            connection_back_left = _np.array(
                [
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_one,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_zero,
                    ],
                ]
            )

            left = _np.array(
                [
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_front_left = _np.array(
                [
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            back = _np.array(
                [
                    [
                        v_one,
                        v_zero,
                        v_one,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_zero,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_one,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_zero,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_back_right = _np.array(
                [
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_one,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_zero,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            bottom = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            top = _np.array(
                [
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                ]
            )

            connection_front_bottom = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            connection_front_top = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                ]
            )

            connection_back_bottom = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_zero,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_back_top = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_one,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_top_right = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one,
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_top_left = _np.array(
                [
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_bottom_left = _np.array(
                [
                    [
                        v_zero,
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            connection_bottom_right = _np.array(
                [
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

        elif closure == "y_max":
            # set points:
            # set points:
            right = _np.array(
                [
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            connection_front_right = _np.array(
                [
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one,
                        v_one,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one,
                        v_zero,
                    ],
                ]
            )

            front = _np.array(
                [
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one,
                        v_one,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one,
                        v_zero,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one,
                        v_one,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one,
                        v_zero,
                    ],
                ]
            )

            connection_back_left = _np.array(
                [
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            left = _np.array(
                [
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_top_left = _np.array(
                [
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one,
                        v_one,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one,
                        v_zero,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            back = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_back_right = _np.array(
                [
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            bottom = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            top = _np.array(
                [
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                ]
            )

            connection_front_bottom = _np.array(
                [
                    [
                        v_one,
                        v_one,
                        v_zero,
                    ],
                    [
                        v_zero,
                        v_one,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            connection_front_top = _np.array(
                [
                    [
                        v_one,
                        v_one,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_one,
                        v_one,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                ]
            )

            connection_back_bottom = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_back_top = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_top_right = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one,
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_front_left = _np.array(
                [
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_bottom_left = _np.array(
                [
                    [
                        v_zero,
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            connection_bottom_right = _np.array(
                [
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

        elif closure == "z_max":
            # set points:
            right = _np.array(
                [
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            connection_front_right = _np.array(
                [
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            front = _np.array(
                [
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            connection_back_left = _np.array(
                [
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            left = _np.array(
                [
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_front_left = _np.array(
                [
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            back = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_back_right = _np.array(
                [
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            bottom = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            top = _np.array(
                [
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one,
                        v_one,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_one,
                        v_one,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_one,
                    ],
                ]
            )

            connection_front_bottom = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            connection_front_top = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_one,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                ]
            )

            connection_back_bottom = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_back_top = _np.array(
                [
                    [
                        v_one,
                        v_zero,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_top_right = _np.array(
                [
                    [
                        v_one,
                        v_one,
                        v_one,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one,
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_top_left = _np.array(
                [
                    [
                        v_zero,
                        v_one,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_bottom_left = _np.array(
                [
                    [
                        v_zero,
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            connection_bottom_right = _np.array(
                [
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

        elif closure == "z_min":
            # set points:
            # set points:
            right = _np.array(
                [
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            connection_front_right = _np.array(
                [
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            front = _np.array(
                [
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            connection_back_left = _np.array(
                [
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            left = _np.array(
                [
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_front_left = _np.array(
                [
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            back = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_back_right = _np.array(
                [
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        -v_wall_thickness + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                ]
            )

            bottom = _np.array(
                [
                    [
                        v_one,
                        v_one,
                        v_zero,
                    ],
                    [
                        v_zero,
                        v_one,
                        v_zero,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_zero,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            top = _np.array(
                [
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                ]
            )

            connection_front_bottom = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one,
                        v_zero,
                    ],
                    [
                        v_zero,
                        v_one,
                        v_zero,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            connection_front_top = _np.array(
                [
                    [
                        v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_one,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                ]
            )

            connection_back_bottom = _np.array(
                [
                    [
                        v_one,
                        v_zero,
                        v_zero,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_zero,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_zero,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                    ],
                ]
            )

            connection_back_top = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        -v_half_contact_length + v_one_half,
                        v_zero,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_top_right = _np.array(
                [
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one,
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one_half + v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half + v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_top_left = _np.array(
                [
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                        v_one,
                    ],
                    [
                        v_zero,
                        v_one_half + v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one_half - v_half_contact_length,
                        v_one_half + v_half_contact_length,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_wall_thickness,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                    ],
                ]
            )

            connection_bottom_left = _np.array(
                [
                    [
                        v_zero,
                        v_one_half + v_half_contact_length,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one_half - v_half_contact_length,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_zero,
                        v_one,
                        v_zero,
                    ],
                    [
                        v_zero,
                        v_zero,
                        v_zero,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_wall_thickness,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

            connection_bottom_right = _np.array(
                [
                    [
                        v_one,
                        v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_one,
                        -v_half_contact_length + v_one_half,
                        v_one_half - v_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_wall_thickness + v_one_half,
                        -v_inner_half_contact_length + v_one_half,
                        v_one_half - v_inner_half_contact_length,
                    ],
                    [
                        v_one,
                        v_one,
                        v_zero,
                    ],
                    [
                        v_one,
                        v_zero,
                        v_zero,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                    [
                        v_one_half + v_inner_half_contact_length,
                        v_one_half - v_inner_half_contact_length,
                        v_one_half - v_wall_thickness,
                    ],
                ]
            )

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_bottom_left)
        )

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_bottom_right)
        )

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_back_bottom)
        )

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_front_bottom)
        )
        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_front_top)
        )

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_back_top)
        )

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_top_right)
        )

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_front_left)
        )

        spline_list.append(_Bezier(degrees=[1, 1, 1], control_points=right))

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_front_right)
        )

        spline_list.append(_Bezier(degrees=[1, 1, 1], control_points=back))

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_back_left)
        )

        spline_list.append(_Bezier(degrees=[1, 1, 1], control_points=left))

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_top_left)
        )

        spline_list.append(_Bezier(degrees=[1, 1, 1], control_points=front))

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_back_right)
        )

        spline_list.append(_Bezier(degrees=[1, 1, 1], control_points=bottom))

        spline_list.append(_Bezier(degrees=[1, 1, 1], control_points=top))

        return (spline_list, None)

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,  # TODO
        contact_length=0.3,
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
            One evaluation point with one parameter is used. This parameter
            describes the thickness of the wall. The parameters must be a
            two-dimensional np.array, where the value must be between 0.01
            and 0.49
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
        derivative_list : list / None
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

        self.check_params(parameters)

        if parameter_sensitivities is not None:
            raise NotImplementedError(
                "Derivatives are not implemented for this tile yet"
            )

        if not (_np.all(parameters > 0) and _np.all(parameters < 0.5)):
            raise ValueError(
                "The thickness of the wall must be in (0.01 and 0.49)"
            )

        if closure is not None:
            return self._closing_tile(
                parameters=parameters,
                parameter_sensitivities=parameter_sensitivities,
                contact_length=contact_length,
                closure=closure,
                **kwargs,
            )

        v_wall_thickness = parameters[0, 0]
        v_zero = 0.0
        v_one_half = 0.5
        v_one = 1.0
        v_half_contact_length = contact_length * 0.5
        v_half_contact_length = contact_length * 0.5
        v_inner_half_contact_length = contact_length * parameters[0, 0]

        spline_list = []

        # set points:
        right = _np.array(
            [
                [
                    v_wall_thickness + v_one_half,
                    -v_inner_half_contact_length + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_one,
                    -v_half_contact_length + v_one_half,
                    v_one_half + v_half_contact_length,
                ],
                [
                    v_wall_thickness + v_one_half,
                    -v_inner_half_contact_length + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    v_one,
                    -v_half_contact_length + v_one_half,
                    v_one_half - v_half_contact_length,
                ],
                [
                    v_wall_thickness + v_one_half,
                    v_inner_half_contact_length + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_one,
                    v_half_contact_length + v_one_half,
                    v_one_half + v_half_contact_length,
                ],
                [
                    v_wall_thickness + v_one_half,
                    v_inner_half_contact_length + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    v_one,
                    v_half_contact_length + v_one_half,
                    v_one_half - v_half_contact_length,
                ],
            ]
        )

        connection_front_right = _np.array(
            [
                [
                    v_wall_thickness + v_one_half,
                    v_inner_half_contact_length + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_one,
                    v_half_contact_length + v_one_half,
                    v_one_half + v_half_contact_length,
                ],
                [
                    v_wall_thickness + v_one_half,
                    v_inner_half_contact_length + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    v_one,
                    v_half_contact_length + v_one_half,
                    v_one_half - v_half_contact_length,
                ],
                [
                    v_inner_half_contact_length + v_one_half,
                    v_wall_thickness + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_half_contact_length + v_one_half,
                    v_one,
                    v_one_half + v_half_contact_length,
                ],
                [
                    v_inner_half_contact_length + v_one_half,
                    v_wall_thickness + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    v_half_contact_length + v_one_half,
                    v_one,
                    v_one_half - v_half_contact_length,
                ],
            ]
        )

        front = _np.array(
            [
                [
                    v_inner_half_contact_length + v_one_half,
                    v_wall_thickness + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_half_contact_length + v_one_half,
                    v_one,
                    v_one_half + v_half_contact_length,
                ],
                [
                    v_inner_half_contact_length + v_one_half,
                    v_wall_thickness + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    v_half_contact_length + v_one_half,
                    v_one,
                    v_one_half - v_half_contact_length,
                ],
                [
                    -v_inner_half_contact_length + v_one_half,
                    v_wall_thickness + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    -v_half_contact_length + v_one_half,
                    v_one,
                    v_one_half + v_half_contact_length,
                ],
                [
                    -v_inner_half_contact_length + v_one_half,
                    v_wall_thickness + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    -v_half_contact_length + v_one_half,
                    v_one,
                    v_one_half - v_half_contact_length,
                ],
            ]
        )

        connection_back_left = _np.array(
            [
                [
                    -v_wall_thickness + v_one_half,
                    -v_inner_half_contact_length + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_zero,
                    -v_half_contact_length + v_one_half,
                    v_one_half + v_half_contact_length,
                ],
                [
                    -v_wall_thickness + v_one_half,
                    -v_inner_half_contact_length + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    v_zero,
                    -v_half_contact_length + v_one_half,
                    v_one_half - v_half_contact_length,
                ],
                [
                    -v_inner_half_contact_length + v_one_half,
                    -v_wall_thickness + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    -v_half_contact_length + v_one_half,
                    v_zero,
                    v_one_half + v_half_contact_length,
                ],
                [
                    -v_inner_half_contact_length + v_one_half,
                    -v_wall_thickness + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    -v_half_contact_length + v_one_half,
                    v_zero,
                    v_one_half - v_half_contact_length,
                ],
            ]
        )

        left = _np.array(
            [
                [
                    v_zero,
                    -v_half_contact_length + v_one_half,
                    v_one_half + v_half_contact_length,
                ],
                [
                    -v_wall_thickness + v_one_half,
                    -v_inner_half_contact_length + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_zero,
                    -v_half_contact_length + v_one_half,
                    v_one_half - v_half_contact_length,
                ],
                [
                    -v_wall_thickness + v_one_half,
                    -v_inner_half_contact_length + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    v_zero,
                    v_half_contact_length + v_one_half,
                    v_one_half + v_half_contact_length,
                ],
                [
                    -v_wall_thickness + v_one_half,
                    v_inner_half_contact_length + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_zero,
                    v_half_contact_length + v_one_half,
                    v_one_half - v_half_contact_length,
                ],
                [
                    -v_wall_thickness + v_one_half,
                    v_inner_half_contact_length + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
            ]
        )

        connection_front_left = _np.array(
            [
                [
                    v_zero,
                    v_half_contact_length + v_one_half,
                    v_one_half + v_half_contact_length,
                ],
                [
                    -v_wall_thickness + v_one_half,
                    v_inner_half_contact_length + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_zero,
                    v_half_contact_length + v_one_half,
                    v_one_half - v_half_contact_length,
                ],
                [
                    -v_wall_thickness + v_one_half,
                    v_inner_half_contact_length + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    -v_half_contact_length + v_one_half,
                    v_one,
                    v_one_half + v_half_contact_length,
                ],
                [
                    -v_inner_half_contact_length + v_one_half,
                    v_wall_thickness + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    -v_half_contact_length + v_one_half,
                    v_one,
                    v_one_half - v_half_contact_length,
                ],
                [
                    -v_inner_half_contact_length + v_one_half,
                    v_wall_thickness + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
            ]
        )

        back = _np.array(
            [
                [
                    v_half_contact_length + v_one_half,
                    v_zero,
                    v_one_half + v_half_contact_length,
                ],
                [
                    v_inner_half_contact_length + v_one_half,
                    -v_wall_thickness + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_half_contact_length + v_one_half,
                    v_zero,
                    v_one_half - v_half_contact_length,
                ],
                [
                    v_inner_half_contact_length + v_one_half,
                    -v_wall_thickness + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    -v_half_contact_length + v_one_half,
                    v_zero,
                    v_one_half + v_half_contact_length,
                ],
                [
                    -v_inner_half_contact_length + v_one_half,
                    -v_wall_thickness + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    -v_half_contact_length + v_one_half,
                    v_zero,
                    v_one_half - v_half_contact_length,
                ],
                [
                    -v_inner_half_contact_length + v_one_half,
                    -v_wall_thickness + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
            ]
        )

        connection_back_right = _np.array(
            [
                [
                    v_inner_half_contact_length + v_one_half,
                    -v_wall_thickness + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_half_contact_length + v_one_half,
                    v_zero,
                    v_one_half + v_half_contact_length,
                ],
                [
                    v_inner_half_contact_length + v_one_half,
                    -v_wall_thickness + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    v_half_contact_length + v_one_half,
                    v_zero,
                    v_one_half - v_half_contact_length,
                ],
                [
                    v_wall_thickness + v_one_half,
                    -v_inner_half_contact_length + v_one_half,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_one,
                    -v_half_contact_length + v_one_half,
                    v_one_half + v_half_contact_length,
                ],
                [
                    v_wall_thickness + v_one_half,
                    -v_inner_half_contact_length + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    v_one,
                    -v_half_contact_length + v_one_half,
                    v_one_half - v_half_contact_length,
                ],
            ]
        )

        bottom = _np.array(
            [
                [
                    v_one_half + v_half_contact_length,
                    v_one_half + v_half_contact_length,
                    v_zero,
                ],
                [
                    v_one_half - v_half_contact_length,
                    v_one_half + v_half_contact_length,
                    v_zero,
                ],
                [
                    v_one_half + v_half_contact_length,
                    v_one_half - v_half_contact_length,
                    v_zero,
                ],
                [
                    v_one_half - v_half_contact_length,
                    v_one_half - v_half_contact_length,
                    v_zero,
                ],
                [
                    v_one_half + v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half - v_wall_thickness,
                ],
                [
                    v_one_half - v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half - v_wall_thickness,
                ],
                [
                    v_one_half + v_inner_half_contact_length,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half - v_wall_thickness,
                ],
                [
                    v_one_half - v_inner_half_contact_length,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half - v_wall_thickness,
                ],
            ]
        )

        top = _np.array(
            [
                [
                    v_one_half + v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half + v_wall_thickness,
                ],
                [
                    v_one_half - v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half + v_wall_thickness,
                ],
                [
                    v_one_half + v_inner_half_contact_length,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half + v_wall_thickness,
                ],
                [
                    v_one_half - v_inner_half_contact_length,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half + v_wall_thickness,
                ],
                [
                    v_one_half + v_half_contact_length,
                    v_one_half + v_half_contact_length,
                    v_one,
                ],
                [
                    v_one_half - v_half_contact_length,
                    v_one_half + v_half_contact_length,
                    v_one,
                ],
                [
                    v_one_half + v_half_contact_length,
                    v_one_half - v_half_contact_length,
                    v_one,
                ],
                [
                    v_one_half - v_half_contact_length,
                    v_one_half - v_half_contact_length,
                    v_one,
                ],
            ]
        )

        connection_front_bottom = _np.array(
            [
                [
                    v_half_contact_length + v_one_half,
                    v_one,
                    v_one_half - v_half_contact_length,
                ],
                [
                    -v_half_contact_length + v_one_half,
                    v_one,
                    v_one_half - v_half_contact_length,
                ],
                [
                    v_one_half + v_half_contact_length,
                    v_one_half + v_half_contact_length,
                    v_zero,
                ],
                [
                    v_one_half - v_half_contact_length,
                    v_one_half + v_half_contact_length,
                    v_zero,
                ],
                [
                    v_inner_half_contact_length + v_one_half,
                    v_one_half + v_wall_thickness,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    -v_inner_half_contact_length + v_one_half,
                    v_one_half + v_wall_thickness,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    v_one_half + v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half - v_wall_thickness,
                ],
                [
                    v_one_half - v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half - v_wall_thickness,
                ],
            ]
        )

        connection_front_top = _np.array(
            [
                [
                    v_half_contact_length + v_one_half,
                    v_one,
                    v_one_half + v_half_contact_length,
                ],
                [
                    -v_half_contact_length + v_one_half,
                    v_one,
                    v_one_half + v_half_contact_length,
                ],
                [
                    v_inner_half_contact_length + v_one_half,
                    v_one_half + v_wall_thickness,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    -v_inner_half_contact_length + v_one_half,
                    v_one_half + v_wall_thickness,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_one_half + v_half_contact_length,
                    v_one_half + v_half_contact_length,
                    v_one,
                ],
                [
                    v_one_half - v_half_contact_length,
                    v_one_half + v_half_contact_length,
                    v_one,
                ],
                [
                    v_one_half + v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half + v_wall_thickness,
                ],
                [
                    v_one_half - v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half + v_wall_thickness,
                ],
            ]
        )

        connection_back_bottom = _np.array(
            [
                [
                    v_one_half + v_half_contact_length,
                    v_one_half - v_half_contact_length,
                    v_zero,
                ],
                [
                    v_one_half - v_half_contact_length,
                    v_one_half - v_half_contact_length,
                    v_zero,
                ],
                [
                    v_one_half + v_half_contact_length,
                    v_zero,
                    v_one_half - v_half_contact_length,
                ],
                [
                    v_one_half - v_half_contact_length,
                    v_zero,
                    v_one_half - v_half_contact_length,
                ],
                [
                    v_one_half + v_inner_half_contact_length,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half - v_wall_thickness,
                ],
                [
                    v_one_half - v_inner_half_contact_length,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half - v_wall_thickness,
                ],
                [
                    v_inner_half_contact_length + v_one_half,
                    v_one_half - v_wall_thickness,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    -v_inner_half_contact_length + v_one_half,
                    v_one_half - v_wall_thickness,
                    v_one_half - v_inner_half_contact_length,
                ],
            ]
        )

        connection_back_top = _np.array(
            [
                [
                    v_one_half + v_half_contact_length,
                    v_one_half - v_half_contact_length,
                    v_one,
                ],
                [
                    v_one_half - v_half_contact_length,
                    v_one_half - v_half_contact_length,
                    v_one,
                ],
                [
                    v_one_half + v_inner_half_contact_length,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half + v_wall_thickness,
                ],
                [
                    v_one_half - v_inner_half_contact_length,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half + v_wall_thickness,
                ],
                [
                    v_half_contact_length + v_one_half,
                    v_zero,
                    v_one_half + v_half_contact_length,
                ],
                [
                    -v_half_contact_length + v_one_half,
                    v_zero,
                    v_one_half + v_half_contact_length,
                ],
                [
                    v_inner_half_contact_length + v_one_half,
                    v_one_half - v_wall_thickness,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    -v_inner_half_contact_length + v_one_half,
                    v_one_half - v_wall_thickness,
                    v_one_half + v_inner_half_contact_length,
                ],
            ]
        )

        connection_top_right = _np.array(
            [
                [
                    v_one_half + v_half_contact_length,
                    v_one_half + v_half_contact_length,
                    v_one,
                ],
                [
                    v_one_half + v_half_contact_length,
                    v_one_half - v_half_contact_length,
                    v_one,
                ],
                [
                    v_one_half + v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half + v_wall_thickness,
                ],
                [
                    v_one_half + v_inner_half_contact_length,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half + v_wall_thickness,
                ],
                [
                    v_one,
                    v_one_half + v_half_contact_length,
                    v_one_half + v_half_contact_length,
                ],
                [
                    v_one,
                    v_one_half - v_half_contact_length,
                    v_one_half + v_half_contact_length,
                ],
                [
                    v_one_half + v_wall_thickness,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_one_half + v_wall_thickness,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                ],
            ]
        )

        connection_top_left = _np.array(
            [
                [
                    v_one_half - v_half_contact_length,
                    v_one_half + v_half_contact_length,
                    v_one,
                ],
                [
                    v_one_half - v_half_contact_length,
                    v_one_half - v_half_contact_length,
                    v_one,
                ],
                [
                    v_zero,
                    v_one_half + v_half_contact_length,
                    v_one_half + v_half_contact_length,
                ],
                [
                    v_zero,
                    v_one_half - v_half_contact_length,
                    v_one_half + v_half_contact_length,
                ],
                [
                    v_one_half - v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half + v_wall_thickness,
                ],
                [
                    v_one_half - v_inner_half_contact_length,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half + v_wall_thickness,
                ],
                [
                    v_one_half - v_wall_thickness,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                ],
                [
                    v_one_half - v_wall_thickness,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                ],
            ]
        )

        connection_bottom_left = _np.array(
            [
                [
                    v_zero,
                    v_one_half + v_half_contact_length,
                    v_one_half - v_half_contact_length,
                ],
                [
                    v_zero,
                    v_one_half - v_half_contact_length,
                    v_one_half - v_half_contact_length,
                ],
                [
                    v_one_half - v_half_contact_length,
                    v_one_half + v_half_contact_length,
                    v_zero,
                ],
                [
                    v_one_half - v_half_contact_length,
                    v_one_half - v_half_contact_length,
                    v_zero,
                ],
                [
                    v_one_half - v_wall_thickness,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    v_one_half - v_wall_thickness,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    v_one_half - v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half - v_wall_thickness,
                ],
                [
                    v_one_half - v_inner_half_contact_length,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half - v_wall_thickness,
                ],
            ]
        )

        connection_bottom_right = _np.array(
            [
                [
                    v_one,
                    v_half_contact_length + v_one_half,
                    v_one_half - v_half_contact_length,
                ],
                [
                    v_one,
                    -v_half_contact_length + v_one_half,
                    v_one_half - v_half_contact_length,
                ],
                [
                    v_wall_thickness + v_one_half,
                    v_inner_half_contact_length + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    v_wall_thickness + v_one_half,
                    -v_inner_half_contact_length + v_one_half,
                    v_one_half - v_inner_half_contact_length,
                ],
                [
                    v_one_half + v_half_contact_length,
                    v_one_half + v_half_contact_length,
                    v_zero,
                ],
                [
                    v_one_half + v_half_contact_length,
                    v_one_half - v_half_contact_length,
                    v_zero,
                ],
                [
                    v_one_half + v_inner_half_contact_length,
                    v_one_half + v_inner_half_contact_length,
                    v_one_half - v_wall_thickness,
                ],
                [
                    v_one_half + v_inner_half_contact_length,
                    v_one_half - v_inner_half_contact_length,
                    v_one_half - v_wall_thickness,
                ],
            ]
        )

        spline_list.append(_Bezier(degrees=[1, 1, 1], control_points=right))

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_front_right)
        )

        spline_list.append(_Bezier(degrees=[1, 1, 1], control_points=back))

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_back_left)
        )

        spline_list.append(_Bezier(degrees=[1, 1, 1], control_points=left))

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_front_left)
        )

        spline_list.append(_Bezier(degrees=[1, 1, 1], control_points=front))

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_back_right)
        )

        spline_list.append(_Bezier(degrees=[1, 1, 1], control_points=bottom))

        spline_list.append(_Bezier(degrees=[1, 1, 1], control_points=top))

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_front_top)
        )

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_front_bottom)
        )

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_back_bottom)
        )

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_back_top)
        )

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_top_right)
        )

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_top_left)
        )

        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_bottom_left)
        )
        spline_list.append(
            _Bezier(degrees=[1, 1, 1], control_points=connection_bottom_right)
        )

        return (spline_list, None)
