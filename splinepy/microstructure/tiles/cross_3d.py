import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class Cross3D(_TileBase):
    def __init__(self):
        """Simple crosstile with linear-quadratic branches and a trilinear
        center spline."""
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
        **kwargs,  # noqa ARG002
    ):
        """Create a closing tile to match with closed surface.

        Parameters
        ----------
        parameters : np.ndarray(6, 1)
          Six evaluation points with one parameter is used. This parameter
          describes the radius of the cylinder at the evaluation point.
          The parameters must be a two-dimensional np.array, where the
          value must be between 0.01 and 0.49
        boundary_width : float
          with of the boundary surrounding branch
        filling_height : float
          portion of the height that is filled in parametric domain
        closure : str
          parametric dimension that needs to be closed.
          Must be {"z_min", "z_max"}

        Returns
        -------
        list_of_splines : list
        derivative_list : list / None
        """
        # Check parameters
        if closure is None:
            raise ValueError("No closing direction given")

        if parameters is None:
            self._logd("Tile request is not parametrized, setting default 0.2")
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

        if not (_np.all(parameters > 0) and _np.all(parameters < 0.5)):
            raise ValueError("Thickness out of range (0, .5)")

        if not (0.0 < float(boundary_width) < 0.5):
            raise ValueError("Boundary Width is out of range")

        if not (0.0 < float(filling_height) < 1.0):
            raise ValueError("Filling must  be in (0,1)")

        inv_boundary_width = 1.0 - boundary_width
        inv_filling_height = 1.0 - filling_height
        center_width = 1.0 - 2 * boundary_width
        ctps_mid_height_top = (1 + filling_height) * 0.5
        ctps_mid_height_bottom = 1.0 - ctps_mid_height_top
        r_center = center_width * 0.5

        spline_list = []
        if closure == "z_min":
            # The branch is located at zmin of current tile
            branch_thickness = parameters[5, 0]
            ctps_corner = _np.array(
                [
                    [0.0, 0.0, 0.0],
                    [boundary_width, 0.0, 0.0],
                    [0.0, boundary_width, 0.0],
                    [boundary_width, boundary_width, 0.0],
                    [0.0, 0.0, filling_height],
                    [boundary_width, 0.0, filling_height],
                    [0.0, boundary_width, filling_height],
                    [boundary_width, boundary_width, filling_height],
                ]
            )

            spline_list.append(
                _Bezier(degrees=[1, 1, 1], control_points=ctps_corner)
            )

            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1],
                    control_points=(
                        ctps_corner + _np.array([0.0, inv_boundary_width, 0.0])
                    ),
                )
            )

            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1],
                    control_points=(
                        ctps_corner + _np.array([inv_boundary_width, 0.0, 0.0])
                    ),
                )
            )

            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1],
                    control_points=(
                        ctps_corner
                        + _np.array(
                            [inv_boundary_width, inv_boundary_width, 0.0]
                        )
                    ),
                )
            )

            center_ctps = _np.array(
                [
                    [boundary_width, boundary_width, 0.0],
                    [inv_boundary_width, boundary_width, 0.0],
                    [boundary_width, inv_boundary_width, 0.0],
                    [inv_boundary_width, inv_boundary_width, 0.0],
                    [boundary_width, boundary_width, filling_height],
                    [inv_boundary_width, boundary_width, filling_height],
                    [boundary_width, inv_boundary_width, filling_height],
                    [inv_boundary_width, inv_boundary_width, filling_height],
                ]
            )

            spline_list.append(
                _Bezier(degrees=[1, 1, 1], control_points=center_ctps)
            )

            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1],
                    control_points=_np.maximum(
                        center_ctps - _np.array([center_width, 0, 0]), 0
                    ),
                )
            )

            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1],
                    control_points=_np.maximum(
                        center_ctps - _np.array([0, center_width, 0]), 0
                    ),
                )
            )

            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1],
                    control_points=_np.minimum(
                        center_ctps + _np.array([center_width, 0, 0]), 1.0
                    ),
                )
            )

            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1],
                    control_points=_np.minimum(
                        center_ctps + _np.array([0, center_width, 0]), 1.0
                    ),
                )
            )
            branch_ctps = _np.array(
                [
                    [-r_center, -r_center, filling_height],
                    [r_center, -r_center, filling_height],
                    [-r_center, r_center, filling_height],
                    [r_center, r_center, filling_height],
                    [
                        -branch_thickness,
                        -branch_thickness,
                        ctps_mid_height_top,
                    ],
                    [branch_thickness, -branch_thickness, ctps_mid_height_top],
                    [-branch_thickness, branch_thickness, ctps_mid_height_top],
                    [branch_thickness, branch_thickness, ctps_mid_height_top],
                    [-branch_thickness, -branch_thickness, 1.0],
                    [branch_thickness, -branch_thickness, 1.0],
                    [-branch_thickness, branch_thickness, 1.0],
                    [branch_thickness, branch_thickness, 1.0],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(degrees=[1, 1, 2], control_points=branch_ctps)
            )

            return (spline_list, None)
        elif closure == "z_max":
            # The branch is located at zmax of current tile
            branch_thickness = parameters[4, 0]
            ctps_corner = _np.array(
                [
                    [0.0, 0.0, inv_filling_height],
                    [boundary_width, 0.0, inv_filling_height],
                    [0.0, boundary_width, inv_filling_height],
                    [boundary_width, boundary_width, inv_filling_height],
                    [0.0, 0.0, 1.0],
                    [boundary_width, 0.0, 1.0],
                    [0.0, boundary_width, 1.0],
                    [boundary_width, boundary_width, 1.0],
                ]
            )

            spline_list.append(
                _Bezier(degrees=[1, 1, 1], control_points=ctps_corner)
            )

            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1],
                    control_points=(
                        ctps_corner + _np.array([0.0, inv_boundary_width, 0.0])
                    ),
                )
            )

            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1],
                    control_points=(
                        ctps_corner + _np.array([inv_boundary_width, 0.0, 0.0])
                    ),
                )
            )

            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1],
                    control_points=(
                        ctps_corner
                        + _np.array(
                            [inv_boundary_width, inv_boundary_width, 0.0]
                        )
                    ),
                )
            )

            center_ctps = _np.array(
                [
                    [boundary_width, boundary_width, inv_filling_height],
                    [inv_boundary_width, boundary_width, inv_filling_height],
                    [boundary_width, inv_boundary_width, inv_filling_height],
                    [
                        inv_boundary_width,
                        inv_boundary_width,
                        inv_filling_height,
                    ],
                    [boundary_width, boundary_width, 1.0],
                    [inv_boundary_width, boundary_width, 1.0],
                    [boundary_width, inv_boundary_width, 1.0],
                    [inv_boundary_width, inv_boundary_width, 1.0],
                ]
            )

            spline_list.append(
                _Bezier(degrees=[1, 1, 1], control_points=center_ctps)
            )

            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1],
                    control_points=_np.maximum(
                        center_ctps - _np.array([center_width, 0, 0]), 0
                    ),
                )
            )

            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1],
                    control_points=_np.maximum(
                        center_ctps - _np.array([0, center_width, 0]), 0
                    ),
                )
            )

            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1],
                    control_points=_np.minimum(
                        center_ctps + _np.array([center_width, 0, 0]), 1.0
                    ),
                )
            )

            spline_list.append(
                _Bezier(
                    degrees=[1, 1, 1],
                    control_points=_np.minimum(
                        center_ctps + _np.array([0, center_width, 0]), 1.0
                    ),
                )
            )

            branch_ctps = _np.array(
                [
                    [-branch_thickness, -branch_thickness, 0.0],
                    [branch_thickness, -branch_thickness, 0.0],
                    [-branch_thickness, branch_thickness, 0.0],
                    [branch_thickness, branch_thickness, 0.0],
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
                    [-r_center, -r_center, inv_filling_height],
                    [r_center, -r_center, inv_filling_height],
                    [-r_center, r_center, inv_filling_height],
                    [r_center, r_center, inv_filling_height],
                ]
            ) + _np.array([0.5, 0.5, 0.0])

            spline_list.append(
                _Bezier(degrees=[1, 1, 2], control_points=branch_ctps)
            )

            return (spline_list, None)
        else:
            raise NotImplementedError(
                "Requested closing dimension is not supported"
            )

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        center_expansion=1.0,
        closure=None,
        **kwargs,  # noqa ARG002
    ):
        """Create a microtile based on the parameters that describe the branch
        thicknesses.

        Thickness parameters are used to describe the inner radius of the
        outward facing branches

        Parameters
        ----------
        parameters : np.array(6, 1)
          Six evaluation points with one parameter is used. This parameter
          describes the radius of the cylinder at the evaluation point.
          The parameters must be a two-dimensional np.array, where the
          value must be between 0.01 and 0.49
        center_expansion : float
          thickness of center is expanded by a factor
        closure : str
          parametric dimension that needs to be closed.
          Must be {"z_min", "z_max"}
        **kwargs
          Will be passed to _closing_tile if given

        Returns
        -------
        microtile_list : list(splines)
        derivative_list : list / None
        """

        if not isinstance(center_expansion, float):
            raise ValueError("Invalid Type")

        if not ((center_expansion > 0.5) and (center_expansion < 1.5)):
            raise ValueError("Center Expansion must be in (.5,1.5)")

        # Max radius, so there is no tanglement in the crosstile
        max_radius = min(0.5, (0.5 / center_expansion))

        # set to default if nothing is given
        if parameters is None:
            self._logd("Setting branch thickness to default 0.2")
            parameters = (
                _np.ones(
                    (
                        self._evaluation_points.shape[0],
                        self._n_info_per_eval_point,
                    )
                )
                * 0.2
            )

        if parameter_sensitivities is not None:
            raise NotImplementedError(
                "Derivatives are not implemented for this tile yet"
            )

        if closure is not None:
            return self._closing_tile(
                parameters=parameters,
                parameter_sensitivities=parameter_sensitivities,  # TODO
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
        ] = parameters.ravel()

        # Check for type and consistency
        self.check_params(parameters)
        if _np.any(parameters <= 0) or _np.any(parameters > max_radius):
            raise ValueError(
                f"Radii must be in (0,{max_radius}) for "
                f"center_expansion {center_expansion}"
            )

        # center radius - mean of radii
        center_r = _np.sum(parameters) / 6.0 * center_expansion
        hd_center = 0.5 * (0.5 + center_r)

        # Initialize return list
        spline_list = []

        # Create the center-tile
        center_points = _np.array(
            [
                [-center_r, -center_r, -center_r],
                [center_r, -center_r, -center_r],
                [-center_r, center_r, -center_r],
                [center_r, center_r, -center_r],
                [-center_r, -center_r, center_r],
                [center_r, -center_r, center_r],
                [-center_r, center_r, center_r],
                [center_r, center_r, center_r],
            ]
        )

        spline_list.append(
            _Bezier(
                degrees=[1, 1, 1],
                control_points=center_points + [0.5, 0.5, 0.5],
            )
        )

        # X-Axis branches
        # X-Min-Branch
        aux_x_min = min(x_min_r, center_r)
        x_min_ctps = _np.array(
            [
                [-0.5, -x_min_r, -x_min_r],
                [-hd_center, -aux_x_min, -aux_x_min],
                center_points[0, :],
                [-0.5, x_min_r, -x_min_r],
                [-hd_center, aux_x_min, -aux_x_min],
                center_points[2, :],
                [-0.5, -x_min_r, x_min_r],
                [-hd_center, -aux_x_min, aux_x_min],
                center_points[4, :],
                [-0.5, x_min_r, x_min_r],
                [-hd_center, aux_x_min, aux_x_min],
                center_points[6, :],
            ]
        )
        spline_list.append(
            _Bezier(
                degrees=[2, 1, 1], control_points=x_min_ctps + [0.5, 0.5, 0.5]
            )
        )

        # X-Min-Branch
        aux_x_max = min(x_max_r, center_r)
        x_max_ctps = _np.array(
            [
                center_points[1, :],
                [hd_center, -aux_x_max, -aux_x_max],
                [0.5, -x_max_r, -x_max_r],
                center_points[3, :],
                [hd_center, aux_x_max, -aux_x_max],
                [0.5, x_max_r, -x_max_r],
                center_points[5, :],
                [hd_center, -aux_x_max, aux_x_max],
                [0.5, -x_max_r, x_max_r],
                center_points[7, :],
                [hd_center, aux_x_max, aux_x_max],
                [0.5, x_max_r, x_max_r],
            ]
        )
        spline_list.append(
            _Bezier(
                degrees=[2, 1, 1], control_points=x_max_ctps + [0.5, 0.5, 0.5]
            )
        )

        # Y-Axis branches
        # Y-Min-Branch
        aux_y_min = min(y_min_r, center_r)
        y_min_ctps = _np.array(
            [
                [-y_min_r, -0.5, -y_min_r],
                [y_min_r, -0.5, -y_min_r],
                [-aux_y_min, -hd_center, -aux_y_min],
                [aux_y_min, -hd_center, -aux_y_min],
                center_points[0, :],
                center_points[1, :],
                [-y_min_r, -0.5, y_min_r],
                [y_min_r, -0.5, y_min_r],
                [-aux_y_min, -hd_center, aux_y_min],
                [aux_y_min, -hd_center, aux_y_min],
                center_points[4, :],
                center_points[5, :],
            ]
        )
        spline_list.append(
            _Bezier(
                degrees=[1, 2, 1], control_points=y_min_ctps + [0.5, 0.5, 0.5]
            )
        )

        # Y-Min-Branch
        aux_y_max = min(y_max_r, center_r)
        y_max_ctps = _np.array(
            [
                center_points[2, :],
                center_points[3, :],
                [-aux_y_max, hd_center, -aux_y_max],
                [aux_y_max, hd_center, -aux_y_max],
                [-y_max_r, 0.5, -y_max_r],
                [y_max_r, 0.5, -y_max_r],
                center_points[6, :],
                center_points[7, :],
                [-aux_y_max, hd_center, aux_y_max],
                [aux_y_max, hd_center, aux_y_max],
                [-y_max_r, 0.5, y_max_r],
                [y_max_r, 0.5, y_max_r],
            ]
        )
        spline_list.append(
            _Bezier(
                degrees=[1, 2, 1], control_points=y_max_ctps + [0.5, 0.5, 0.5]
            )
        )

        # Y-Axis branches
        # Y-Min-Branch
        aux_z_min = min(z_min_r, center_r)
        z_min_ctps = _np.array(
            [
                [-z_min_r, -z_min_r, -0.5],
                [z_min_r, -z_min_r, -0.5],
                [-z_min_r, z_min_r, -0.5],
                [z_min_r, z_min_r, -0.5],
                [-aux_z_min, -aux_z_min, -hd_center],
                [aux_z_min, -aux_z_min, -hd_center],
                [-aux_z_min, aux_z_min, -hd_center],
                [aux_z_min, aux_z_min, -hd_center],
                center_points[0, :],
                center_points[1, :],
                center_points[2, :],
                center_points[3, :],
            ]
        )
        spline_list.append(
            _Bezier(
                degrees=[1, 1, 2], control_points=z_min_ctps + [0.5, 0.5, 0.5]
            )
        )

        # Y-Min-Branch
        aux_z_max = min(z_max_r, center_r)
        z_max_ctps = _np.array(
            [
                center_points[4, :],
                center_points[5, :],
                center_points[6, :],
                center_points[7, :],
                [-aux_z_max, -aux_z_max, hd_center],
                [aux_z_max, -aux_z_max, hd_center],
                [-aux_z_max, aux_z_max, hd_center],
                [aux_z_max, aux_z_max, hd_center],
                [-z_max_r, -z_max_r, 0.5],
                [z_max_r, -z_max_r, 0.5],
                [-z_max_r, z_max_r, 0.5],
                [z_max_r, z_max_r, 0.5],
            ]
        )
        spline_list.append(
            _Bezier(
                degrees=[1, 1, 2], control_points=z_max_ctps + [0.5, 0.5, 0.5]
            )
        )

        return (spline_list, None)
