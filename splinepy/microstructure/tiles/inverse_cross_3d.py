import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class InverseCross3D(_TileBase):
    def __init__(self):
        """Simple inverse crosstile to tile with linear-quadratic branches and
        a trilinear center spline.

        Class that provides necessary functions to create inverse microtile,
        that can be used to describe the domain surrounding cross-tile
        microstructure.
        """
        self._dim = 3
        self._para_dim = 3
        self._evaluation_points = _np.array(
            [
                [0.0, 0.5, 0.5],
                [1.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 1.0, 0.5],
                [0.5, 0.5, 0.0],
                [0.5, 0.5, 1.0],
            ]
        )
        self._n_info_per_eval_point = 1

    def _closing_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        closure=None,
        boundary_width=0.1,
        filling_height=0.5,
        seperator_distance=None,
        **kwargs,  # noqa ARG002
    ):
        """Create a closing tile to match with closed surface.

        Parameters
        ----------
        parameters : tuple(np.ndarray)
          radii of fitting cylinder at evaluation points
        boundary_width : float
          with of the boundary surrounding branch
        filling_height : float
          portion of the height that is filled in parametric domain
        seperator_distance : float
          Describes the position of the separator layer in the control point
          domain
        closure : str
          parametric dimension that needs to be closed.
          Must be one of {"z_min", "z_max"}

        Returns
        -------
        list_of_splines : list
        derivatives: list<list<splines>> / None
        """
        # Check parameters
        if closure is None:
            raise ValueError("No closing direction given")

        # Set default values
        if seperator_distance is None:
            seperator_distance = 0.3

        if parameters is None:
            self._logd("Tile request is not parametrized, setting default 0.2")
            parameters = (
                _np.ones(
                    (
                        self._evaluation_points.shape[0],
                        self._n_info_per_eval_point,
                    )
                )
                * 0.2
            )

        self.check_params(parameters)

        if parameter_sensitivities is not None:
            raise NotImplementedError(
                "Derivatives are not implemented for this tile yet"
            )

        if not (_np.all(parameters > 0) and _np.all(parameters < 0.5)):
            raise ValueError("Thickness out of range (0, .5)")

        if not (0.0 < float(boundary_width) < 0.5):
            raise ValueError("Boundary Width is out of range")

        if not (0.0 < float(filling_height) < 1.0):
            raise ValueError("Filling must  be in (0,1)")

        # Precompute auxiliary values
        inv_filling_height = 1.0 - filling_height
        ctps_mid_height_top = (1 + filling_height) * 0.5
        ctps_mid_height_bottom = 1.0 - ctps_mid_height_top
        center_width = 1.0 - 2 * boundary_width
        r_center = center_width * 0.5
        half_r_center = (r_center + 0.5) * 0.5
        aux_column_width = 0.5 - 2 * (0.5 - seperator_distance)

        spline_list = []
        if closure == "z_min":
            branch_thickness = parameters.flatten()[5]
            branch_neighbor_x_min_ctps = _np.array(
                [
                    [-0.5, -r_center, filling_height],
                    [-half_r_center, -r_center, filling_height],
                    [-r_center, -r_center, filling_height],
                    [-0.5, r_center, filling_height],
                    [-half_r_center, r_center, filling_height],
                    [-r_center, r_center, filling_height],
                    [-0.5, -aux_column_width, ctps_mid_height_top],
                    [
                        -seperator_distance,
                        -aux_column_width,
                        ctps_mid_height_top,
                    ],
                    [
                        -branch_thickness,
                        -branch_thickness,
                        ctps_mid_height_top,
                    ],
                    [-0.5, aux_column_width, ctps_mid_height_top],
                    [
                        -seperator_distance,
                        aux_column_width,
                        ctps_mid_height_top,
                    ],
                    [-branch_thickness, branch_thickness, ctps_mid_height_top],
                    [-0.5, -aux_column_width, 1.0],
                    [-seperator_distance, -aux_column_width, 1.0],
                    [-branch_thickness, -branch_thickness, 1.0],
                    [-0.5, aux_column_width, 1.0],
                    [-seperator_distance, aux_column_width, 1.0],
                    [-branch_thickness, branch_thickness, 1.0],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[2, 1, 2],
                    control_points=branch_neighbor_x_min_ctps,
                )
            )

            branch_neighbor_x_max_ctps = _np.array(
                [
                    [r_center, -r_center, filling_height],
                    [half_r_center, -r_center, filling_height],
                    [0.5, -r_center, filling_height],
                    [r_center, r_center, filling_height],
                    [half_r_center, r_center, filling_height],
                    [0.5, r_center, filling_height],
                    [branch_thickness, -branch_thickness, ctps_mid_height_top],
                    [
                        seperator_distance,
                        -aux_column_width,
                        ctps_mid_height_top,
                    ],
                    [0.5, -aux_column_width, ctps_mid_height_top],
                    [branch_thickness, branch_thickness, ctps_mid_height_top],
                    [
                        seperator_distance,
                        aux_column_width,
                        ctps_mid_height_top,
                    ],
                    [0.5, aux_column_width, ctps_mid_height_top],
                    [branch_thickness, -branch_thickness, 1.0],
                    [seperator_distance, -aux_column_width, 1.0],
                    [0.5, -aux_column_width, 1.0],
                    [branch_thickness, branch_thickness, 1.0],
                    [seperator_distance, aux_column_width, 1.0],
                    [0.5, aux_column_width, 1.0],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[2, 1, 2],
                    control_points=branch_neighbor_x_max_ctps,
                )
            )

            branch_neighbor_y_min_ctps = _np.array(
                [
                    [-r_center, -0.5, filling_height],
                    [r_center, -0.5, filling_height],
                    [-r_center, -half_r_center, filling_height],
                    [r_center, -half_r_center, filling_height],
                    [-r_center, -r_center, filling_height],
                    [r_center, -r_center, filling_height],
                    [-aux_column_width, -0.5, ctps_mid_height_top],
                    [aux_column_width, -0.5, ctps_mid_height_top],
                    [
                        -aux_column_width,
                        -seperator_distance,
                        ctps_mid_height_top,
                    ],
                    [
                        aux_column_width,
                        -seperator_distance,
                        ctps_mid_height_top,
                    ],
                    [
                        -branch_thickness,
                        -branch_thickness,
                        ctps_mid_height_top,
                    ],
                    [branch_thickness, -branch_thickness, ctps_mid_height_top],
                    [-aux_column_width, -0.5, 1.0],
                    [aux_column_width, -0.5, 1.0],
                    [-aux_column_width, -seperator_distance, 1.0],
                    [aux_column_width, -seperator_distance, 1.0],
                    [-branch_thickness, -branch_thickness, 1.0],
                    [branch_thickness, -branch_thickness, 1.0],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[1, 2, 2],
                    control_points=branch_neighbor_y_min_ctps,
                )
            )

            branch_neighbor_y_max_ctps = _np.array(
                [
                    [-r_center, r_center, filling_height],
                    [r_center, r_center, filling_height],
                    [-r_center, half_r_center, filling_height],
                    [r_center, half_r_center, filling_height],
                    [-r_center, 0.5, filling_height],
                    [r_center, 0.5, filling_height],
                    [-branch_thickness, branch_thickness, ctps_mid_height_top],
                    [branch_thickness, branch_thickness, ctps_mid_height_top],
                    [
                        -aux_column_width,
                        seperator_distance,
                        ctps_mid_height_top,
                    ],
                    [
                        aux_column_width,
                        seperator_distance,
                        ctps_mid_height_top,
                    ],
                    [-aux_column_width, 0.5, ctps_mid_height_top],
                    [aux_column_width, 0.5, ctps_mid_height_top],
                    [-branch_thickness, branch_thickness, 1.0],
                    [branch_thickness, branch_thickness, 1.0],
                    [-aux_column_width, seperator_distance, 1.0],
                    [aux_column_width, seperator_distance, 1.0],
                    [-aux_column_width, 0.5, 1.0],
                    [aux_column_width, 0.5, 1.0],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[1, 2, 2],
                    control_points=branch_neighbor_y_max_ctps,
                )
            )

            branch_x_min_y_min_ctps = _np.array(
                [
                    [-0.5, -0.5, filling_height],
                    [-half_r_center, -0.5, filling_height],
                    [-r_center, -0.5, filling_height],
                    [-0.5, -half_r_center, filling_height],
                    [-half_r_center, -half_r_center, filling_height],
                    [-r_center, -half_r_center, filling_height],
                    [-0.5, -r_center, filling_height],
                    [-half_r_center, -r_center, filling_height],
                    [-r_center, -r_center, filling_height],
                    [-0.5, -0.5, ctps_mid_height_top],
                    [-seperator_distance, -0.5, ctps_mid_height_top],
                    [-aux_column_width, -0.5, ctps_mid_height_top],
                    [-0.5, -seperator_distance, ctps_mid_height_top],
                    [
                        -seperator_distance,
                        -seperator_distance,
                        ctps_mid_height_top,
                    ],
                    [
                        -aux_column_width,
                        -seperator_distance,
                        ctps_mid_height_top,
                    ],
                    [-0.5, -aux_column_width, ctps_mid_height_top],
                    [
                        -seperator_distance,
                        -aux_column_width,
                        ctps_mid_height_top,
                    ],
                    [
                        -branch_thickness,
                        -branch_thickness,
                        ctps_mid_height_top,
                    ],
                    [-0.5, -0.5, 1.0],
                    [-seperator_distance, -0.5, 1.0],
                    [-aux_column_width, -0.5, 1.0],
                    [-0.5, -seperator_distance, 1.0],
                    [-seperator_distance, -seperator_distance, 1.0],
                    [-aux_column_width, -seperator_distance, 1.0],
                    [-0.5, -aux_column_width, 1.0],
                    [-seperator_distance, -aux_column_width, 1.0],
                    [-branch_thickness, -branch_thickness, 1.0],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[2, 2, 2], control_points=branch_x_min_y_min_ctps
                )
            )

            branch_x_min_y_max_ctps = _np.array(
                [
                    [-0.5, r_center, filling_height],
                    [-half_r_center, r_center, filling_height],
                    [-r_center, r_center, filling_height],
                    [-0.5, half_r_center, filling_height],
                    [-half_r_center, half_r_center, filling_height],
                    [-r_center, half_r_center, filling_height],
                    [-0.5, 0.5, filling_height],
                    [-half_r_center, 0.5, filling_height],
                    [-r_center, 0.5, filling_height],
                    [-0.5, aux_column_width, ctps_mid_height_top],
                    [
                        -seperator_distance,
                        aux_column_width,
                        ctps_mid_height_top,
                    ],
                    [-branch_thickness, branch_thickness, ctps_mid_height_top],
                    [-0.5, seperator_distance, ctps_mid_height_top],
                    [
                        -seperator_distance,
                        seperator_distance,
                        ctps_mid_height_top,
                    ],
                    [
                        -aux_column_width,
                        seperator_distance,
                        ctps_mid_height_top,
                    ],
                    [-0.5, 0.5, ctps_mid_height_top],
                    [-seperator_distance, 0.5, ctps_mid_height_top],
                    [-aux_column_width, 0.5, ctps_mid_height_top],
                    [-0.5, aux_column_width, 1.0],
                    [-seperator_distance, aux_column_width, 1.0],
                    [-branch_thickness, branch_thickness, 1.0],
                    [-0.5, seperator_distance, 1.0],
                    [-seperator_distance, seperator_distance, 1.0],
                    [-aux_column_width, seperator_distance, 1.0],
                    [-0.5, 0.5, 1.0],
                    [-seperator_distance, 0.5, 1.0],
                    [-aux_column_width, 0.5, 1.0],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[2, 2, 2], control_points=branch_x_min_y_max_ctps
                )
            )

            branch_x_max_y_min_ctps = _np.array(
                [
                    [r_center, -0.5, filling_height],
                    [half_r_center, -0.5, filling_height],
                    [0.5, -0.5, filling_height],
                    [r_center, -half_r_center, filling_height],
                    [half_r_center, -half_r_center, filling_height],
                    [0.5, -half_r_center, filling_height],
                    [r_center, -r_center, filling_height],
                    [half_r_center, -r_center, filling_height],
                    [0.5, -r_center, filling_height],
                    [aux_column_width, -0.5, ctps_mid_height_top],
                    [seperator_distance, -0.5, ctps_mid_height_top],
                    [0.5, -0.5, ctps_mid_height_top],
                    [
                        aux_column_width,
                        -seperator_distance,
                        ctps_mid_height_top,
                    ],
                    [
                        seperator_distance,
                        -seperator_distance,
                        ctps_mid_height_top,
                    ],
                    [0.5, -seperator_distance, ctps_mid_height_top],
                    [branch_thickness, -branch_thickness, ctps_mid_height_top],
                    [
                        seperator_distance,
                        -aux_column_width,
                        ctps_mid_height_top,
                    ],
                    [0.5, -aux_column_width, ctps_mid_height_top],
                    [aux_column_width, -0.5, 1.0],
                    [seperator_distance, -0.5, 1.0],
                    [0.5, -0.5, 1.0],
                    [aux_column_width, -seperator_distance, 1.0],
                    [seperator_distance, -seperator_distance, 1.0],
                    [0.5, -seperator_distance, 1.0],
                    [branch_thickness, -branch_thickness, 1.0],
                    [seperator_distance, -aux_column_width, 1.0],
                    [0.5, -aux_column_width, 1.0],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[2, 2, 2], control_points=branch_x_max_y_min_ctps
                )
            )

            branch_x_max_y_max_ctps = _np.array(
                [
                    [r_center, r_center, filling_height],
                    [half_r_center, r_center, filling_height],
                    [0.5, r_center, filling_height],
                    [r_center, half_r_center, filling_height],
                    [half_r_center, half_r_center, filling_height],
                    [0.5, half_r_center, filling_height],
                    [r_center, 0.5, filling_height],
                    [half_r_center, 0.5, filling_height],
                    [0.5, 0.5, filling_height],
                    [branch_thickness, branch_thickness, ctps_mid_height_top],
                    [
                        seperator_distance,
                        aux_column_width,
                        ctps_mid_height_top,
                    ],
                    [0.5, aux_column_width, ctps_mid_height_top],
                    [
                        aux_column_width,
                        seperator_distance,
                        ctps_mid_height_top,
                    ],
                    [
                        seperator_distance,
                        seperator_distance,
                        ctps_mid_height_top,
                    ],
                    [0.5, seperator_distance, ctps_mid_height_top],
                    [aux_column_width, 0.5, ctps_mid_height_top],
                    [seperator_distance, 0.5, ctps_mid_height_top],
                    [0.5, 0.5, ctps_mid_height_top],
                    [branch_thickness, branch_thickness, 1.0],
                    [seperator_distance, aux_column_width, 1.0],
                    [0.5, aux_column_width, 1.0],
                    [aux_column_width, seperator_distance, 1.0],
                    [seperator_distance, seperator_distance, 1.0],
                    [0.5, seperator_distance, 1.0],
                    [aux_column_width, 0.5, 1.0],
                    [seperator_distance, 0.5, 1.0],
                    [0.5, 0.5, 1.0],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[2, 2, 2], control_points=branch_x_max_y_max_ctps
                )
            )

            return (spline_list, None)

        elif closure == "z_max":
            branch_thickness = parameters.flatten()[4]
            branch_neighbor_x_min_ctps = _np.array(
                [
                    [-0.5, -aux_column_width, 0.0],
                    [-seperator_distance, -aux_column_width, 0.0],
                    [-branch_thickness, -branch_thickness, 0.0],
                    [-0.5, aux_column_width, 0.0],
                    [-seperator_distance, aux_column_width, 0.0],
                    [-branch_thickness, branch_thickness, 0.0],
                    [-0.5, -aux_column_width, ctps_mid_height_bottom],
                    [
                        -seperator_distance,
                        -aux_column_width,
                        ctps_mid_height_bottom,
                    ],
                    [
                        -branch_thickness,
                        -branch_thickness,
                        ctps_mid_height_bottom,
                    ],
                    [-0.5, aux_column_width, ctps_mid_height_bottom],
                    [
                        -seperator_distance,
                        aux_column_width,
                        ctps_mid_height_bottom,
                    ],
                    [
                        -branch_thickness,
                        branch_thickness,
                        ctps_mid_height_bottom,
                    ],
                    [-0.5, -r_center, inv_filling_height],
                    [-half_r_center, -r_center, inv_filling_height],
                    [-r_center, -r_center, inv_filling_height],
                    [-0.5, r_center, inv_filling_height],
                    [-half_r_center, r_center, inv_filling_height],
                    [-r_center, r_center, inv_filling_height],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[2, 1, 2],
                    control_points=branch_neighbor_x_min_ctps,
                )
            )

            branch_neighbor_x_max_ctps = _np.array(
                [
                    [branch_thickness, -branch_thickness, 0.0],
                    [seperator_distance, -aux_column_width, 0.0],
                    [0.5, -aux_column_width, 0.0],
                    [branch_thickness, branch_thickness, 0.0],
                    [seperator_distance, aux_column_width, 0.0],
                    [0.5, aux_column_width, 0.0],
                    [
                        branch_thickness,
                        -branch_thickness,
                        ctps_mid_height_bottom,
                    ],
                    [
                        seperator_distance,
                        -aux_column_width,
                        ctps_mid_height_bottom,
                    ],
                    [0.5, -aux_column_width, ctps_mid_height_bottom],
                    [
                        branch_thickness,
                        branch_thickness,
                        ctps_mid_height_bottom,
                    ],
                    [
                        seperator_distance,
                        aux_column_width,
                        ctps_mid_height_bottom,
                    ],
                    [0.5, aux_column_width, ctps_mid_height_bottom],
                    [r_center, -r_center, inv_filling_height],
                    [half_r_center, -r_center, inv_filling_height],
                    [0.5, -r_center, inv_filling_height],
                    [r_center, r_center, inv_filling_height],
                    [half_r_center, r_center, inv_filling_height],
                    [0.5, r_center, inv_filling_height],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[2, 1, 2],
                    control_points=branch_neighbor_x_max_ctps,
                )
            )

            branch_neighbor_y_min_ctps = _np.array(
                [
                    [-aux_column_width, -0.5, 0.0],
                    [aux_column_width, -0.5, 0.0],
                    [-aux_column_width, -seperator_distance, 0.0],
                    [aux_column_width, -seperator_distance, 0.0],
                    [-branch_thickness, -branch_thickness, 0.0],
                    [branch_thickness, -branch_thickness, 0.0],
                    [-aux_column_width, -0.5, ctps_mid_height_bottom],
                    [aux_column_width, -0.5, ctps_mid_height_bottom],
                    [
                        -aux_column_width,
                        -seperator_distance,
                        ctps_mid_height_bottom,
                    ],
                    [
                        aux_column_width,
                        -seperator_distance,
                        ctps_mid_height_bottom,
                    ],
                    [
                        -branch_thickness,
                        -branch_thickness,
                        ctps_mid_height_bottom,
                    ],
                    [
                        branch_thickness,
                        -branch_thickness,
                        ctps_mid_height_bottom,
                    ],
                    [-r_center, -0.5, inv_filling_height],
                    [r_center, -0.5, inv_filling_height],
                    [-r_center, -half_r_center, inv_filling_height],
                    [r_center, -half_r_center, inv_filling_height],
                    [-r_center, -r_center, inv_filling_height],
                    [r_center, -r_center, inv_filling_height],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[1, 2, 2],
                    control_points=branch_neighbor_y_min_ctps,
                )
            )

            branch_neighbor_y_max_ctps = _np.array(
                [
                    [-branch_thickness, branch_thickness, 0.0],
                    [branch_thickness, branch_thickness, 0.0],
                    [-aux_column_width, seperator_distance, 0.0],
                    [aux_column_width, seperator_distance, 0.0],
                    [-aux_column_width, 0.5, 0.0],
                    [aux_column_width, 0.5, 0.0],
                    [
                        -branch_thickness,
                        branch_thickness,
                        ctps_mid_height_bottom,
                    ],
                    [
                        branch_thickness,
                        branch_thickness,
                        ctps_mid_height_bottom,
                    ],
                    [
                        -aux_column_width,
                        seperator_distance,
                        ctps_mid_height_bottom,
                    ],
                    [
                        aux_column_width,
                        seperator_distance,
                        ctps_mid_height_bottom,
                    ],
                    [-aux_column_width, 0.5, ctps_mid_height_bottom],
                    [aux_column_width, 0.5, ctps_mid_height_bottom],
                    [-r_center, r_center, inv_filling_height],
                    [r_center, r_center, inv_filling_height],
                    [-r_center, half_r_center, inv_filling_height],
                    [r_center, half_r_center, inv_filling_height],
                    [-r_center, 0.5, inv_filling_height],
                    [r_center, 0.5, inv_filling_height],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[1, 2, 2],
                    control_points=branch_neighbor_y_max_ctps,
                )
            )

            branch_x_min_y_min_ctps = _np.array(
                [
                    [-0.5, -0.5, 0.0],
                    [-seperator_distance, -0.5, 0.0],
                    [-aux_column_width, -0.5, 0.0],
                    [-0.5, -seperator_distance, 0.0],
                    [-seperator_distance, -seperator_distance, 0.0],
                    [-aux_column_width, -seperator_distance, 0.0],
                    [-0.5, -aux_column_width, 0.0],
                    [-seperator_distance, -aux_column_width, 0.0],
                    [-branch_thickness, -branch_thickness, 0.0],
                    [-0.5, -0.5, ctps_mid_height_bottom],
                    [-seperator_distance, -0.5, ctps_mid_height_bottom],
                    [-aux_column_width, -0.5, ctps_mid_height_bottom],
                    [-0.5, -seperator_distance, ctps_mid_height_bottom],
                    [
                        -seperator_distance,
                        -seperator_distance,
                        ctps_mid_height_bottom,
                    ],
                    [
                        -aux_column_width,
                        -seperator_distance,
                        ctps_mid_height_bottom,
                    ],
                    [-0.5, -aux_column_width, ctps_mid_height_bottom],
                    [
                        -seperator_distance,
                        -aux_column_width,
                        ctps_mid_height_bottom,
                    ],
                    [
                        -branch_thickness,
                        -branch_thickness,
                        ctps_mid_height_bottom,
                    ],
                    [-0.5, -0.5, inv_filling_height],
                    [-half_r_center, -0.5, inv_filling_height],
                    [-r_center, -0.5, inv_filling_height],
                    [-0.5, -half_r_center, inv_filling_height],
                    [-half_r_center, -half_r_center, inv_filling_height],
                    [-r_center, -half_r_center, inv_filling_height],
                    [-0.5, -r_center, inv_filling_height],
                    [-half_r_center, -r_center, inv_filling_height],
                    [-r_center, -r_center, inv_filling_height],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[2, 2, 2], control_points=branch_x_min_y_min_ctps
                )
            )

            branch_x_max_y_max_ctps = _np.array(
                [
                    [branch_thickness, branch_thickness, 0.0],
                    [seperator_distance, aux_column_width, 0.0],
                    [0.5, aux_column_width, 0.0],
                    [aux_column_width, seperator_distance, 0.0],
                    [seperator_distance, seperator_distance, 0.0],
                    [0.5, seperator_distance, 0.0],
                    [aux_column_width, 0.5, 0.0],
                    [seperator_distance, 0.5, 0.0],
                    [0.5, 0.5, 0.0],
                    [
                        branch_thickness,
                        branch_thickness,
                        ctps_mid_height_bottom,
                    ],
                    [
                        seperator_distance,
                        aux_column_width,
                        ctps_mid_height_bottom,
                    ],
                    [0.5, aux_column_width, ctps_mid_height_bottom],
                    [
                        aux_column_width,
                        seperator_distance,
                        ctps_mid_height_bottom,
                    ],
                    [
                        seperator_distance,
                        seperator_distance,
                        ctps_mid_height_bottom,
                    ],
                    [0.5, seperator_distance, ctps_mid_height_bottom],
                    [aux_column_width, 0.5, ctps_mid_height_bottom],
                    [seperator_distance, 0.5, ctps_mid_height_bottom],
                    [0.5, 0.5, ctps_mid_height_bottom],
                    [r_center, r_center, inv_filling_height],
                    [half_r_center, r_center, inv_filling_height],
                    [0.5, r_center, inv_filling_height],
                    [r_center, half_r_center, inv_filling_height],
                    [half_r_center, half_r_center, inv_filling_height],
                    [0.5, half_r_center, inv_filling_height],
                    [r_center, 0.5, inv_filling_height],
                    [half_r_center, 0.5, inv_filling_height],
                    [0.5, 0.5, inv_filling_height],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[2, 2, 2], control_points=branch_x_max_y_max_ctps
                )
            )

            branch_x_max_y_min_ctps = _np.array(
                [
                    [aux_column_width, -0.5, 0.0],
                    [seperator_distance, -0.5, 0.0],
                    [0.5, -0.5, 0.0],
                    [aux_column_width, -seperator_distance, 0.0],
                    [seperator_distance, -seperator_distance, 0.0],
                    [0.5, -seperator_distance, 0.0],
                    [branch_thickness, -branch_thickness, 0.0],
                    [seperator_distance, -aux_column_width, 0.0],
                    [0.5, -aux_column_width, 0.0],
                    [aux_column_width, -0.5, ctps_mid_height_bottom],
                    [seperator_distance, -0.5, ctps_mid_height_bottom],
                    [0.5, -0.5, ctps_mid_height_bottom],
                    [
                        aux_column_width,
                        -seperator_distance,
                        ctps_mid_height_bottom,
                    ],
                    [
                        seperator_distance,
                        -seperator_distance,
                        ctps_mid_height_bottom,
                    ],
                    [0.5, -seperator_distance, ctps_mid_height_bottom],
                    [
                        branch_thickness,
                        -branch_thickness,
                        ctps_mid_height_bottom,
                    ],
                    [
                        seperator_distance,
                        -aux_column_width,
                        ctps_mid_height_bottom,
                    ],
                    [0.5, -aux_column_width, ctps_mid_height_bottom],
                    [r_center, -0.5, inv_filling_height],
                    [half_r_center, -0.5, inv_filling_height],
                    [0.5, -0.5, inv_filling_height],
                    [r_center, -half_r_center, inv_filling_height],
                    [half_r_center, -half_r_center, inv_filling_height],
                    [0.5, -half_r_center, inv_filling_height],
                    [r_center, -r_center, inv_filling_height],
                    [half_r_center, -r_center, inv_filling_height],
                    [0.5, -r_center, inv_filling_height],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[2, 2, 2], control_points=branch_x_max_y_min_ctps
                )
            )

            branch_x_min_y_max_ctps = _np.array(
                [
                    [-0.5, aux_column_width, 0.0],
                    [-seperator_distance, aux_column_width, 0.0],
                    [-branch_thickness, branch_thickness, 0.0],
                    [-0.5, seperator_distance, 0.0],
                    [-seperator_distance, seperator_distance, 0.0],
                    [-aux_column_width, seperator_distance, 0.0],
                    [-0.5, 0.5, 0.0],
                    [-seperator_distance, 0.5, 0.0],
                    [-aux_column_width, 0.5, 0.0],
                    [-0.5, aux_column_width, ctps_mid_height_bottom],
                    [
                        -seperator_distance,
                        aux_column_width,
                        ctps_mid_height_bottom,
                    ],
                    [
                        -branch_thickness,
                        branch_thickness,
                        ctps_mid_height_bottom,
                    ],
                    [-0.5, seperator_distance, ctps_mid_height_bottom],
                    [
                        -seperator_distance,
                        seperator_distance,
                        ctps_mid_height_bottom,
                    ],
                    [
                        -aux_column_width,
                        seperator_distance,
                        ctps_mid_height_bottom,
                    ],
                    [-0.5, 0.5, ctps_mid_height_bottom],
                    [-seperator_distance, 0.5, ctps_mid_height_bottom],
                    [-aux_column_width, 0.5, ctps_mid_height_bottom],
                    [-0.5, r_center, inv_filling_height],
                    [-half_r_center, r_center, inv_filling_height],
                    [-r_center, r_center, inv_filling_height],
                    [-0.5, half_r_center, inv_filling_height],
                    [-half_r_center, half_r_center, inv_filling_height],
                    [-r_center, half_r_center, inv_filling_height],
                    [-0.5, 0.5, inv_filling_height],
                    [-half_r_center, 0.5, inv_filling_height],
                    [-r_center, 0.5, inv_filling_height],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(
                    degrees=[2, 2, 2], control_points=branch_x_min_y_max_ctps
                )
            )

            return (spline_list, None)
        else:
            raise ValueError("Corner Type not supported")

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        seperator_distance=None,
        center_expansion=1.0,
        closure=None,
        **kwargs,  # noqa ARG002
    ):
        """Create an inverse microtile based on the parameters that describe
        the branch thicknesses.

        Thickness parameters are used to describe the inner radius of the
        outward facing branches

        Parameters
        ----------
        parameters : np.array
          only first entry is used, defines the internal radii of the
          branches
        seperator_distance : float
          Control point distance to separation layer of biquadratic degrees,
          determines the minimum branch thickness (defaults to 0.4)
        center_expansion : float
          thickness of center is expanded by a factor (default to 1.0), which
          determines the maximum branch thickness
        closure : str
          parametric dimension that needs to be closed.
          Must be one of {"z_min", "z_max"}

        Returns
        -------
        microtile_list : list(splines)
        """

        if not isinstance(center_expansion, float):
            raise ValueError("Invalid Type")

        if not ((center_expansion > 0.5) and (center_expansion < 1.5)):
            raise ValueError("Center Expansion must be in (.5, 1.5)")

        # Set default values
        if seperator_distance is None:
            seperator_distance = 0.3

        if center_expansion is None:
            center_expansion = 1.0

        # Check if all radii are in allowed range
        max_radius = min(0.5, (0.5 / center_expansion))
        max_radius = min(max_radius, seperator_distance)
        min_radius = max(0.5 - seperator_distance, 0)

        # set to default if nothing is given
        if parameters is None:
            self._logd("Setting branch thickness to default 0.2")
            parameters = (
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

        if _np.any(parameters < min_radius) or _np.any(
            parameters > max_radius
        ):
            raise ValueError(
                f"Radii must be in (0,{max_radius}) for "
                f"center_expansion {center_expansion}"
            )

        if closure is not None:
            return self._closing_tile(
                parameters=parameters,
                parameter_sensitivities=parameter_sensitivities,
                seperator_distance=seperator_distance,
                closure=closure,
                **kwargs,
            )

        [
            x_min_r,
            x_max_r,
            y_min_r,
            y_max_r,
            z_min_r,
            z_max_r,
        ] = parameters.flatten()

        # center radius
        center_r = _np.sum(parameters) / 6.0 * center_expansion

        # Auxiliary values for smooothing (mid-branch thickness)
        [
            aux_x_min,
            aux_x_max,
            aux_y_min,
            aux_y_max,
            aux_z_min,
            aux_z_max,
        ] = _np.minimum(parameters.ravel(), center_r)
        # Branch midlength
        hd_center = 0.5 * (0.5 + center_r)
        aux_column_width = 0.5 - 2 * (0.5 - seperator_distance)

        # Init return type
        spline_list = []

        # Start with branch interconnections
        x_min_y_min = _np.array(
            [
                [-0.5, -0.5, -aux_column_width],
                [-seperator_distance, -0.5, -aux_column_width],
                [-y_min_r, -0.5, -y_min_r],
                [-0.5, -seperator_distance, -aux_column_width],
                [-hd_center, -hd_center, -aux_column_width],
                [-aux_y_min, -hd_center, -aux_y_min],
                [-0.5, -x_min_r, -x_min_r],
                [-hd_center, -aux_x_min, -aux_x_min],
                [-center_r, -center_r, -center_r],
                [-0.5, -0.5, aux_column_width],
                [-seperator_distance, -0.5, aux_column_width],
                [-y_min_r, -0.5, y_min_r],
                [-0.5, -seperator_distance, aux_column_width],
                [-hd_center, -hd_center, aux_column_width],
                [-aux_y_min, -hd_center, aux_y_min],
                [-0.5, -x_min_r, x_min_r],
                [-hd_center, -aux_x_min, aux_x_min],
                [-center_r, -center_r, center_r],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 2, 1], control_points=x_min_y_min)
        )

        x_max_y_min = _np.array(
            [
                [y_min_r, -0.5, -y_min_r],
                [seperator_distance, -0.5, -aux_column_width],
                [0.5, -0.5, -aux_column_width],
                [aux_y_min, -hd_center, -aux_y_min],
                [hd_center, -hd_center, -aux_column_width],
                [0.5, -seperator_distance, -aux_column_width],
                [center_r, -center_r, -center_r],
                [hd_center, -aux_x_max, -aux_x_max],
                [0.5, -x_max_r, -x_max_r],
                [y_min_r, -0.5, y_min_r],
                [seperator_distance, -0.5, aux_column_width],
                [0.5, -0.5, aux_column_width],
                [aux_y_min, -hd_center, aux_y_min],
                [hd_center, -hd_center, aux_column_width],
                [0.5, -seperator_distance, aux_column_width],
                [center_r, -center_r, center_r],
                [hd_center, -aux_x_max, aux_x_max],
                [0.5, -x_max_r, x_max_r],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 2, 1], control_points=x_max_y_min)
        )

        x_min_y_max = _np.array(
            [
                [-0.5, x_min_r, -x_min_r],
                [-hd_center, aux_x_min, -aux_x_min],
                [-center_r, center_r, -center_r],
                [-0.5, seperator_distance, -aux_column_width],
                [-hd_center, hd_center, -aux_column_width],
                [-aux_y_max, hd_center, -aux_y_max],
                [-0.5, 0.5, -aux_column_width],
                [-seperator_distance, 0.5, -aux_column_width],
                [-y_max_r, 0.5, -y_max_r],
                [-0.5, x_min_r, x_min_r],
                [-hd_center, aux_x_min, aux_x_min],
                [-center_r, center_r, center_r],
                [-0.5, seperator_distance, aux_column_width],
                [-hd_center, hd_center, aux_column_width],
                [-aux_y_max, hd_center, aux_y_max],
                [-0.5, 0.5, aux_column_width],
                [-seperator_distance, 0.5, aux_column_width],
                [-y_max_r, 0.5, y_max_r],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 2, 1], control_points=x_min_y_max)
        )

        x_max_y_max = _np.array(
            [
                [center_r, center_r, -center_r],
                [hd_center, aux_x_max, -aux_x_max],
                [0.5, x_max_r, -x_max_r],
                [aux_y_max, hd_center, -aux_y_max],
                [hd_center, hd_center, -aux_column_width],
                [0.5, seperator_distance, -aux_column_width],
                [y_max_r, 0.5, -y_max_r],
                [seperator_distance, 0.5, -aux_column_width],
                [0.5, 0.5, -aux_column_width],
                [center_r, center_r, center_r],
                [hd_center, aux_x_max, aux_x_max],
                [0.5, x_max_r, x_max_r],
                [aux_y_max, hd_center, aux_y_max],
                [hd_center, hd_center, aux_column_width],
                [0.5, seperator_distance, aux_column_width],
                [y_max_r, 0.5, y_max_r],
                [seperator_distance, 0.5, aux_column_width],
                [0.5, 0.5, aux_column_width],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 2, 1], control_points=x_max_y_max)
        )

        x_min_z_min = _np.array(
            [
                [-0.5, -aux_column_width, -0.5],
                [-seperator_distance, -aux_column_width, -0.5],
                [-z_min_r, -z_min_r, -0.5],
                [-0.5, aux_column_width, -0.5],
                [-seperator_distance, aux_column_width, -0.5],
                [-z_min_r, z_min_r, -0.5],
                [-0.5, -aux_column_width, -seperator_distance],
                [-hd_center, -aux_column_width, -hd_center],
                [-aux_z_min, -aux_z_min, -hd_center],
                [-0.5, aux_column_width, -seperator_distance],
                [-hd_center, aux_column_width, -hd_center],
                [-aux_z_min, aux_z_min, -hd_center],
                [-0.5, -x_min_r, -x_min_r],
                [-hd_center, -aux_x_min, -aux_x_min],
                [-center_r, -center_r, -center_r],
                [-0.5, x_min_r, -x_min_r],
                [-hd_center, aux_x_min, -aux_x_min],
                [-center_r, center_r, -center_r],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 1, 2], control_points=x_min_z_min)
        )

        x_max_z_min = _np.array(
            [
                [z_min_r, -z_min_r, -0.5],
                [seperator_distance, -aux_column_width, -0.5],
                [0.5, -aux_column_width, -0.5],
                [z_min_r, z_min_r, -0.5],
                [seperator_distance, aux_column_width, -0.5],
                [0.5, aux_column_width, -0.5],
                [aux_z_min, -aux_z_min, -hd_center],
                [hd_center, -aux_column_width, -hd_center],
                [0.5, -aux_column_width, -seperator_distance],
                [aux_z_min, aux_z_min, -hd_center],
                [hd_center, aux_column_width, -hd_center],
                [0.5, aux_column_width, -seperator_distance],
                [center_r, -center_r, -center_r],
                [hd_center, -aux_x_max, -aux_x_max],
                [0.5, -x_max_r, -x_max_r],
                [center_r, center_r, -center_r],
                [hd_center, aux_x_max, -aux_x_max],
                [0.5, x_max_r, -x_max_r],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 1, 2], control_points=x_max_z_min)
        )

        x_min_z_max = _np.array(
            [
                [-0.5, -x_min_r, x_min_r],
                [-hd_center, -aux_x_min, aux_x_min],
                [-center_r, -center_r, center_r],
                [-0.5, x_min_r, x_min_r],
                [-hd_center, aux_x_min, aux_x_min],
                [-center_r, center_r, center_r],
                [-0.5, -aux_column_width, seperator_distance],
                [-hd_center, -aux_column_width, hd_center],
                [-aux_z_max, -aux_z_max, hd_center],
                [-0.5, aux_column_width, seperator_distance],
                [-hd_center, aux_column_width, hd_center],
                [-aux_z_max, aux_z_max, hd_center],
                [-0.5, -aux_column_width, 0.5],
                [-seperator_distance, -aux_column_width, 0.5],
                [-z_max_r, -z_max_r, 0.5],
                [-0.5, aux_column_width, 0.5],
                [-seperator_distance, aux_column_width, 0.5],
                [-z_max_r, z_max_r, 0.5],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 1, 2], control_points=x_min_z_max)
        )

        x_max_z_max = _np.array(
            [
                [center_r, -center_r, center_r],
                [hd_center, -aux_x_max, aux_x_max],
                [0.5, -x_max_r, x_max_r],
                [center_r, center_r, center_r],
                [hd_center, aux_x_max, aux_x_max],
                [0.5, x_max_r, x_max_r],
                [aux_z_max, -aux_z_max, hd_center],
                [hd_center, -aux_column_width, hd_center],
                [0.5, -aux_column_width, seperator_distance],
                [aux_z_max, aux_z_max, hd_center],
                [hd_center, aux_column_width, hd_center],
                [0.5, aux_column_width, seperator_distance],
                [z_max_r, -z_max_r, 0.5],
                [seperator_distance, -aux_column_width, 0.5],
                [0.5, -aux_column_width, 0.5],
                [z_max_r, z_max_r, 0.5],
                [seperator_distance, aux_column_width, 0.5],
                [0.5, aux_column_width, 0.5],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 1, 2], control_points=x_max_z_max)
        )

        y_min_z_min = _np.array(
            [
                [-aux_column_width, -0.5, -0.5],
                [aux_column_width, -0.5, -0.5],
                [-aux_column_width, -seperator_distance, -0.5],
                [aux_column_width, -seperator_distance, -0.5],
                [-z_min_r, -z_min_r, -0.5],
                [z_min_r, -z_min_r, -0.5],
                [-aux_column_width, -0.5, -seperator_distance],
                [aux_column_width, -0.5, -seperator_distance],
                [-aux_column_width, -hd_center, -hd_center],
                [aux_column_width, -hd_center, -hd_center],
                [-aux_z_min, -aux_z_min, -hd_center],
                [aux_z_min, -aux_z_min, -hd_center],
                [-y_min_r, -0.5, -y_min_r],
                [y_min_r, -0.5, -y_min_r],
                [-aux_y_min, -hd_center, -aux_y_min],
                [aux_y_min, -hd_center, -aux_y_min],
                [-center_r, -center_r, -center_r],
                [center_r, -center_r, -center_r],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[1, 2, 2], control_points=y_min_z_min)
        )

        y_max_z_min = _np.array(
            [
                [-z_min_r, z_min_r, -0.5],
                [z_min_r, z_min_r, -0.5],
                [-aux_column_width, seperator_distance, -0.5],
                [aux_column_width, seperator_distance, -0.5],
                [-aux_column_width, 0.5, -0.5],
                [aux_column_width, 0.5, -0.5],
                [-aux_z_min, aux_z_min, -hd_center],
                [aux_z_min, aux_z_min, -hd_center],
                [-aux_column_width, hd_center, -hd_center],
                [aux_column_width, hd_center, -hd_center],
                [-aux_column_width, 0.5, -seperator_distance],
                [aux_column_width, 0.5, -seperator_distance],
                [-center_r, center_r, -center_r],
                [center_r, center_r, -center_r],
                [-aux_y_max, hd_center, -aux_y_max],
                [aux_y_max, hd_center, -aux_y_max],
                [-y_max_r, 0.5, -y_max_r],
                [y_max_r, 0.5, -y_max_r],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[1, 2, 2], control_points=y_max_z_min)
        )

        y_min_z_max = _np.array(
            [
                [-y_min_r, -0.5, y_min_r],
                [y_min_r, -0.5, y_min_r],
                [-aux_y_min, -hd_center, aux_y_min],
                [aux_y_min, -hd_center, aux_y_min],
                [-center_r, -center_r, center_r],
                [center_r, -center_r, center_r],
                [-aux_column_width, -0.5, seperator_distance],
                [aux_column_width, -0.5, seperator_distance],
                [-aux_column_width, -hd_center, hd_center],
                [aux_column_width, -hd_center, hd_center],
                [-aux_z_max, -aux_z_max, hd_center],
                [aux_z_max, -aux_z_max, hd_center],
                [-aux_column_width, -0.5, 0.5],
                [aux_column_width, -0.5, 0.5],
                [-aux_column_width, -seperator_distance, 0.5],
                [aux_column_width, -seperator_distance, 0.5],
                [-z_max_r, -z_max_r, 0.5],
                [z_max_r, -z_max_r, 0.5],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[1, 2, 2], control_points=y_min_z_max)
        )

        y_max_z_max = _np.array(
            [
                [-center_r, center_r, center_r],
                [center_r, center_r, center_r],
                [-aux_y_max, hd_center, aux_y_max],
                [aux_y_max, hd_center, aux_y_max],
                [-y_max_r, 0.5, y_max_r],
                [y_max_r, 0.5, y_max_r],
                [-aux_z_max, aux_z_max, hd_center],
                [aux_z_max, aux_z_max, hd_center],
                [-aux_column_width, hd_center, hd_center],
                [aux_column_width, hd_center, hd_center],
                [-aux_column_width, 0.5, seperator_distance],
                [aux_column_width, 0.5, seperator_distance],
                [-z_max_r, z_max_r, 0.5],
                [z_max_r, z_max_r, 0.5],
                [-aux_column_width, seperator_distance, 0.5],
                [aux_column_width, seperator_distance, 0.5],
                [-aux_column_width, 0.5, 0.5],
                [aux_column_width, 0.5, 0.5],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[1, 2, 2], control_points=y_max_z_max)
        )

        x_min_y_min_z_min = _np.array(
            [
                [-0.5, -0.5, -0.5],
                [-seperator_distance, -0.5, -0.5],
                [-aux_column_width, -0.5, -0.5],
                [-0.5, -seperator_distance, -0.5],
                [-seperator_distance, -seperator_distance, -0.5],
                [-aux_column_width, -seperator_distance, -0.5],
                [-0.5, -aux_column_width, -0.5],
                [-seperator_distance, -aux_column_width, -0.5],
                [-z_min_r, -z_min_r, -0.5],
                [-0.5, -0.5, -seperator_distance],
                [-seperator_distance, -0.5, -seperator_distance],
                [-aux_column_width, -0.5, -seperator_distance],
                [-0.5, -seperator_distance, -seperator_distance],
                [-hd_center, -hd_center, -hd_center],
                [-aux_column_width, -hd_center, -hd_center],
                [-0.5, -aux_column_width, -seperator_distance],
                [-hd_center, -aux_column_width, -hd_center],
                [-aux_z_min, -aux_z_min, -hd_center],
                [-0.5, -0.5, -aux_column_width],
                [-seperator_distance, -0.5, -aux_column_width],
                [-y_min_r, -0.5, -y_min_r],
                [-0.5, -seperator_distance, -aux_column_width],
                [-hd_center, -hd_center, -aux_column_width],
                [-aux_y_min, -hd_center, -aux_y_min],
                [-0.5, -x_min_r, -x_min_r],
                [-hd_center, -aux_x_min, -aux_x_min],
                [-center_r, -center_r, -center_r],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 2, 2], control_points=x_min_y_min_z_min)
        )

        x_max_y_min_z_min = _np.array(
            [
                [aux_column_width, -0.5, -0.5],
                [seperator_distance, -0.5, -0.5],
                [0.5, -0.5, -0.5],
                [aux_column_width, -seperator_distance, -0.5],
                [seperator_distance, -seperator_distance, -0.5],
                [0.5, -seperator_distance, -0.5],
                [z_min_r, -z_min_r, -0.5],
                [seperator_distance, -aux_column_width, -0.5],
                [0.5, -aux_column_width, -0.5],
                [aux_column_width, -0.5, -seperator_distance],
                [seperator_distance, -0.5, -seperator_distance],
                [0.5, -0.5, -seperator_distance],
                [aux_column_width, -hd_center, -hd_center],
                [hd_center, -hd_center, -hd_center],
                [0.5, -seperator_distance, -seperator_distance],
                [aux_z_min, -aux_z_min, -hd_center],
                [hd_center, -aux_column_width, -hd_center],
                [0.5, -aux_column_width, -seperator_distance],
                [y_min_r, -0.5, -y_min_r],
                [seperator_distance, -0.5, -aux_column_width],
                [0.5, -0.5, -aux_column_width],
                [aux_y_min, -hd_center, -aux_y_min],
                [hd_center, -hd_center, -aux_column_width],
                [0.5, -seperator_distance, -aux_column_width],
                [center_r, -center_r, -center_r],
                [hd_center, -aux_x_max, -aux_x_max],
                [0.5, -x_max_r, -x_max_r],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 2, 2], control_points=x_max_y_min_z_min)
        )

        x_min_y_max_z_min = _np.array(
            [
                [-0.5, aux_column_width, -0.5],
                [-seperator_distance, aux_column_width, -0.5],
                [-z_min_r, z_min_r, -0.5],
                [-0.5, seperator_distance, -0.5],
                [-seperator_distance, seperator_distance, -0.5],
                [-aux_column_width, seperator_distance, -0.5],
                [-0.5, 0.5, -0.5],
                [-seperator_distance, 0.5, -0.5],
                [-aux_column_width, 0.5, -0.5],
                [-0.5, aux_column_width, -seperator_distance],
                [-hd_center, aux_column_width, -hd_center],
                [-aux_z_min, aux_z_min, -hd_center],
                [-0.5, seperator_distance, -seperator_distance],
                [-hd_center, hd_center, -hd_center],
                [-aux_column_width, hd_center, -hd_center],
                [-0.5, 0.5, -seperator_distance],
                [-seperator_distance, 0.5, -seperator_distance],
                [-aux_column_width, 0.5, -seperator_distance],
                [-0.5, x_min_r, -x_min_r],
                [-hd_center, aux_x_min, -aux_x_min],
                [-center_r, center_r, -center_r],
                [-0.5, seperator_distance, -aux_column_width],
                [-hd_center, hd_center, -aux_column_width],
                [-aux_y_max, hd_center, -aux_y_max],
                [-0.5, 0.5, -aux_column_width],
                [-seperator_distance, 0.5, -aux_column_width],
                [-y_max_r, 0.5, -y_max_r],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 2, 2], control_points=x_min_y_max_z_min)
        )

        x_max_y_max_z_min = _np.array(
            [
                [z_min_r, z_min_r, -0.5],
                [seperator_distance, aux_column_width, -0.5],
                [0.5, aux_column_width, -0.5],
                [aux_column_width, seperator_distance, -0.5],
                [seperator_distance, seperator_distance, -0.5],
                [0.5, seperator_distance, -0.5],
                [aux_column_width, 0.5, -0.5],
                [seperator_distance, 0.5, -0.5],
                [0.5, 0.5, -0.5],
                [aux_z_min, aux_z_min, -hd_center],
                [hd_center, aux_column_width, -hd_center],
                [0.5, aux_column_width, -seperator_distance],
                [aux_column_width, hd_center, -hd_center],
                [hd_center, hd_center, -hd_center],
                [0.5, seperator_distance, -seperator_distance],
                [aux_column_width, 0.5, -seperator_distance],
                [seperator_distance, 0.5, -seperator_distance],
                [0.5, 0.5, -seperator_distance],
                [center_r, center_r, -center_r],
                [hd_center, aux_x_max, -aux_x_max],
                [0.5, x_max_r, -x_max_r],
                [aux_y_max, hd_center, -aux_y_max],
                [hd_center, hd_center, -aux_column_width],
                [0.5, seperator_distance, -aux_column_width],
                [y_max_r, 0.5, -y_max_r],
                [seperator_distance, 0.5, -aux_column_width],
                [0.5, 0.5, -aux_column_width],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 2, 2], control_points=x_max_y_max_z_min)
        )

        x_min_y_min_z_max = _np.array(
            [
                [-0.5, -0.5, aux_column_width],
                [-seperator_distance, -0.5, aux_column_width],
                [-y_min_r, -0.5, y_min_r],
                [-0.5, -seperator_distance, aux_column_width],
                [-hd_center, -hd_center, aux_column_width],
                [-aux_y_min, -hd_center, aux_y_min],
                [-0.5, -x_min_r, x_min_r],
                [-hd_center, -aux_x_min, aux_x_min],
                [-center_r, -center_r, center_r],
                [-0.5, -0.5, seperator_distance],
                [-seperator_distance, -0.5, seperator_distance],
                [-aux_column_width, -0.5, seperator_distance],
                [-0.5, -seperator_distance, seperator_distance],
                [-hd_center, -hd_center, hd_center],
                [-aux_column_width, -hd_center, hd_center],
                [-0.5, -aux_column_width, seperator_distance],
                [-hd_center, -aux_column_width, hd_center],
                [-aux_z_max, -aux_z_max, hd_center],
                [-0.5, -0.5, 0.5],
                [-seperator_distance, -0.5, 0.5],
                [-aux_column_width, -0.5, 0.5],
                [-0.5, -seperator_distance, 0.5],
                [-seperator_distance, -seperator_distance, 0.5],
                [-aux_column_width, -seperator_distance, 0.5],
                [-0.5, -aux_column_width, 0.5],
                [-seperator_distance, -aux_column_width, 0.5],
                [-z_max_r, -z_max_r, 0.5],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 2, 2], control_points=x_min_y_min_z_max)
        )

        x_max_y_min_z_max = _np.array(
            [
                [y_min_r, -0.5, y_min_r],
                [seperator_distance, -0.5, aux_column_width],
                [0.5, -0.5, aux_column_width],
                [aux_y_min, -hd_center, aux_y_min],
                [hd_center, -hd_center, aux_column_width],
                [0.5, -seperator_distance, aux_column_width],
                [center_r, -center_r, center_r],
                [hd_center, -aux_x_max, aux_x_max],
                [0.5, -x_max_r, x_max_r],
                [aux_column_width, -0.5, seperator_distance],
                [seperator_distance, -0.5, seperator_distance],
                [0.5, -0.5, seperator_distance],
                [aux_column_width, -hd_center, hd_center],
                [hd_center, -hd_center, hd_center],
                [0.5, -seperator_distance, seperator_distance],
                [aux_z_max, -aux_z_max, hd_center],
                [hd_center, -aux_column_width, hd_center],
                [0.5, -aux_column_width, seperator_distance],
                [aux_column_width, -0.5, 0.5],
                [seperator_distance, -0.5, 0.5],
                [0.5, -0.5, 0.5],
                [aux_column_width, -seperator_distance, 0.5],
                [seperator_distance, -seperator_distance, 0.5],
                [0.5, -seperator_distance, 0.5],
                [z_max_r, -z_max_r, 0.5],
                [seperator_distance, -aux_column_width, 0.5],
                [0.5, -aux_column_width, 0.5],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 2, 2], control_points=x_max_y_min_z_max)
        )

        x_min_y_max_z_max = _np.array(
            [
                [-0.5, x_min_r, x_min_r],
                [-hd_center, aux_x_min, aux_x_min],
                [-center_r, center_r, center_r],
                [-0.5, seperator_distance, aux_column_width],
                [-hd_center, hd_center, aux_column_width],
                [-aux_y_max, hd_center, aux_y_max],
                [-0.5, 0.5, aux_column_width],
                [-seperator_distance, 0.5, aux_column_width],
                [-y_max_r, 0.5, y_max_r],
                [-0.5, aux_column_width, seperator_distance],
                [-hd_center, aux_column_width, hd_center],
                [-aux_z_max, aux_z_max, hd_center],
                [-0.5, seperator_distance, seperator_distance],
                [-hd_center, hd_center, hd_center],
                [-aux_column_width, hd_center, hd_center],
                [-0.5, 0.5, seperator_distance],
                [-seperator_distance, 0.5, seperator_distance],
                [-aux_column_width, 0.5, seperator_distance],
                [-0.5, aux_column_width, 0.5],
                [-seperator_distance, aux_column_width, 0.5],
                [-z_max_r, z_max_r, 0.5],
                [-0.5, seperator_distance, 0.5],
                [-seperator_distance, seperator_distance, 0.5],
                [-aux_column_width, seperator_distance, 0.5],
                [-0.5, 0.5, 0.5],
                [-seperator_distance, 0.5, 0.5],
                [-aux_column_width, 0.5, 0.5],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 2, 2], control_points=x_min_y_max_z_max)
        )

        x_max_y_max_z_max = _np.array(
            [
                [center_r, center_r, center_r],
                [hd_center, aux_x_max, aux_x_max],
                [0.5, x_max_r, x_max_r],
                [aux_y_max, hd_center, aux_y_max],
                [hd_center, hd_center, aux_column_width],
                [0.5, seperator_distance, aux_column_width],
                [y_max_r, 0.5, y_max_r],
                [seperator_distance, 0.5, aux_column_width],
                [0.5, 0.5, aux_column_width],
                [aux_z_max, aux_z_max, hd_center],
                [hd_center, aux_column_width, hd_center],
                [0.5, aux_column_width, seperator_distance],
                [aux_column_width, hd_center, hd_center],
                [hd_center, hd_center, hd_center],
                [0.5, seperator_distance, seperator_distance],
                [aux_column_width, 0.5, seperator_distance],
                [seperator_distance, 0.5, seperator_distance],
                [0.5, 0.5, seperator_distance],
                [z_max_r, z_max_r, 0.5],
                [seperator_distance, aux_column_width, 0.5],
                [0.5, aux_column_width, 0.5],
                [aux_column_width, seperator_distance, 0.5],
                [seperator_distance, seperator_distance, 0.5],
                [0.5, seperator_distance, 0.5],
                [aux_column_width, 0.5, 0.5],
                [seperator_distance, 0.5, 0.5],
                [0.5, 0.5, 0.5],
            ]
        ) + _np.array([0.5, 0.5, 0.5])

        spline_list.append(
            _Bezier(degrees=[2, 2, 2], control_points=x_max_y_max_z_max)
        )
        return (spline_list, None)
